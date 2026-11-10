# How the pipeline works

This describes what each stage computes, the statistical model, and where the
outputs go. For running it, see [GETTING_STARTED.md](GETTING_STARTED.md); for each
setting, [CONFIG.md](CONFIG.md).

## The idea

A TWAS asks whether people whose genotypes predict higher expression of a gene also
have a higher trait value. Its power depends on how well expression can be predicted
from cis-SNPs, how much trait variance acts through the gene, how big the GWAS is,
and whether the trait actually works through the gene as a whole or through one of
its isoforms. isoQTL_sim makes all of this known by construction, then measures how
often each kind of test detects it.

The one thing not simulated is the genotypes: they are real, so linkage
disequilibrium, allele frequencies and the number of SNPs near each gene are real.

## Stages

`isoqtl-sim plan` lists the stages for a config. In order:

| Stage | Jobs | What it does |
|---|---|---|
| `prepare` (run separately) | 1 | Reads the GTF and transcript FASTA, keeps protein-coding autosomal genes, samples `n_genes` of them, draws per-gene parameters, lists the samples from the VCF header. Writes `params/`. |
| `simulate` | genes ÷ `genes_per_job` | Per gene: reads genotypes in the training window, simulates isoform expression and one GWAS per scenario. Writes `simulated/<gene>.rds` and `simulated_status/`. |
| `truth_matrix` | 1 | Only with reads. Joins all genes' isoform expression into one samples × isoforms matrix. |
| `reads` | samples | Only with reads. Simulates a FASTQ pair per sample with `Rsubread::simReads` from the true isoform abundances. |
| `salmon_index_<annot>` | 1 | Only with reads and no `salmon_index` given. Builds the Salmon index for that annotation. |
| `salmon_<annot>` | samples | Quantifies each sample (`salmon quant -l A --validateMappings --seqBias`). |
| `aggregate_salmon_<annot>` | 1 | `tximport` to isoform and gene counts, filter by `min_prop_expressed`, TMM-normalize, log2. Writes `quants/salmon_<annot>_{gene,isoform}.rds` and a BED file. |
| `star`, `featurecounts_<annot>` | samples, 1 | Only if `quant.star` is set. Genome alignment and gene counting. **Never run yet.** |
| `fastqc`, `multiqc` | samples, 1 | Only if `quant.fastqc: true`. **Never run yet.** |
| `train_<source>` | genes ÷ `genes_per_job` | Per gene: trains the models on that source's expression, tests each against every scenario's GWAS. One stage per source: `truth`, then `salmon_<annot>` and `star_<annot>` for each quantification. Writes `twas/<source>/chunk_*.tsv.gz`. |
| `summarize` | 1 | Combines isoform tests per gene, computes power and calibration. Writes `summary/`. |

Each job writes a marker in `.done/<stage>/<index>` when it finishes. `run` and
`submit` skip jobs with markers, so an interrupted run continues where it stopped,
and a failed job can be rerun alone. On LSF, each stage's array starts when the
previous one has finished successfully (`-w done(...)`).

## The simulation model

Everything below is in `R/sim_core.R`, which has no dependencies on the pipeline and
is covered by `tests/test_sim_core.R`.

**Genes and isoforms.** A gene has M isoforms, those in the transcript FASTA for its
gene ID. About 14% of protein-coding genes have one isoform; the mean is 8 and the
maximum 192 (GENCODE v38).

**Causal SNPs.** Each isoform has `n_causal` causal SNPs within `windows.causal` of
the gene. A share `prop_shared` of them (drawn per gene from the configured list,
rounded up) is common to every isoform; the rest are specific to each isoform. With
`n_causal: 4` and `prop_shared: 0.75`, a gene with 3 isoforms has 3 shared and 3
specific causal SNPs, 6 in all. A gene whose window has fewer SNPs than it needs is
recorded as `too_few_snps` and left out.

**Effects.** Shared SNPs' effects on the M isoforms are drawn from a multivariate
normal with a random correlation matrix whose entries lie in `cov_B` (0.7–0.9 by
default: the shared SNPs push isoforms in the same direction). Specific SNPs' effects
are independent normals. Effects on each isoform are scaled so that the genetic
component has variance `h2_g` in the sample; `h2_g` is drawn per gene from a Beta
distribution (default Beta(2, 20), mean 0.09) or fixed. Residual noise is
multivariate normal with variance `1 − h2_g` and correlations in `cov_E`.
Expression is on a log scale: `Y = X B + E`, samples × isoforms. Gene expression is
`log(sum(exp(Y)))`, so the gene's heritability is an emergent quantity, reported as
`gene_h2`.

**Traits and GWAS.** A scenario names a trait mediator (`gene`: the trait acts
through the sum of the isoforms' genetic values; `isoform`: through one isoform,
chosen at random per gene), an `h2ge` (share of trait variance explained through the
gene) and a GWAS size `n_gwas`. The trait's per-SNP effects are the gene's causal
effects, weighted accordingly and scaled to `h2ge`. GWAS marginal Z-scores for every
SNP in the training window are then drawn as

    Z ~ N( sqrt(n_gwas / (1 − h2ge)) · R β,  R )

where R is the LD matrix of the real panel. Noise is generated as `Xᵀe / sqrt(n−1)`
with standardized X, whose covariance is exactly R even though there are more SNPs
than samples. With `h2ge: 0` the scenario is a null trait and measures the
false-positive rate. A randomly signed effect avoids a directional bias.

**Prediction models.** For each gene (total expression) and each isoform, the
pipeline fits every method in `models.methods` on the standardized genotypes in the
training window (`windows.train`, 1 Mb by default) with 5-fold cross-validation:
`enet` (elastic net, glmnet, α = 0.5), `blup` (rrBLUP) and `susie` (susieR). The
method with the best CV R² is kept. If `models.isotwas` is set, the gene's isoforms
are also modeled jointly with `isotwas::compute_isotwas` using the listed methods,
and the best per isoform is kept as a separate set of models. A single-isoform gene's
isoform model is its gene model.

**Tests.** Each model is tested against each scenario's GWAS with the standard
summary-statistic TWAS statistic, `Z = wᵀz / sqrt(wᵀ R w)`. A model with all-zero
weights (a null model: nothing was predictable) is recorded but not tested. The
expected TWAS Z for a true effect is about `sqrt(n_gwas · h2ge)` times the
correlation between predicted and true expression, which is why `h2ge` values of
0.0001–0.001 with `n_gwas` of 100,000 span the interesting range.

**Combining isoforms.** For the `isoform_acat` and `isotwas_acat` levels, the
p-values of a gene's isoform tests are combined with ACAT (the Cauchy combination
test), which is valid under arbitrary dependence between the isoforms. Null models
are excluded first; a p-value of exactly 1 would otherwise dominate the combination.

**Summary.** For each source, level and scenario, `rejection_rate` is the share of
*simulated* genes with p below `alpha`, and `rejection_rate_bonferroni` uses
`alpha / n_genes_simulated`. A gene with no testable model at a level counts as a
miss. Genes that could not be simulated are reported separately and excluded.

## Read-level comparison

With `run_reads: true`, the true isoform abundances (`exp(Y)`, scaled to a million
per sample) drive `Rsubread::simReads` to produce a FASTQ pair per sample with the
configured read length and library size. Only the simulated genes are expressed, so
`library_size` should be scaled to the number of genes (see CONFIG.md). Salmon
quantifies each sample against each annotation in `quant.salmon`; `tximport`
summarizes to genes; features with TPM > 0.1 in fewer than `min_prop_expressed` of
samples are dropped; counts are TMM-normalized and log2-transformed. The training
stage then runs again on these estimates, giving `salmon_<annot>` rows beside the
`truth` rows. The difference is what quantification (and, with a different
annotation, annotation mismatch) costs in prediction accuracy and power.

Measured on 100 genes: Salmon cost 2–5 points of power at the gene level and up to
11 at the isoform level; isoform CV R² fell by 28%, gene CV R² by 8%.

## Output directory

```
<output_dir>/
  params/
    config.lock.yaml        the settings this run was prepared with (checked on every command)
    genes.tsv               chosen genes: coordinates, n_isoforms, prop_shared, h2_g, n_causal
    transcripts.tsv         their transcripts
    samples.txt             sample IDs, in the order used everywhere
  simulated/<gene>.rds      Y, B, causal SNP IDs, per-scenario GWAS z, isoform and gene h2
  simulated_status/         one TSV per chunk: gene_id, status (ok / no_snps / too_few_snps)
  simulated_all/            isoform_log_expression.rds (with reads)
  reads/<sample>_R{1,2}.fastq.gz          (with reads)
  salmon_index/<annot>/, salmon/<annot>/<sample>/quant.sf
  quants/<source>_{gene,isoform}.rds      normalized expression matrices
  bed/<source>_gene.bed                   gene expression in QTLtools BED layout
  twas/<source>/chunk_*.tsv.gz            one row per model × scenario: method, cv_r2, twas_z, twas_p, ...
  summary/power_summary.tsv               the result table
  summary/gene_tests.tsv.gz               one row per source × level × scenario × gene
  logs/<stage>/<index>.log                job output
  .done/<stage>/<index>                   completion markers
  submit_lsf.sh                           written by `submit`
```

## Reproducibility

Every random draw is seeded from `seed` plus the gene chunk or sample index, so the
same config gives the same genes, samples, expression, GWAS and models, and jobs can
run in any order or be rerun alone. `params/config.lock.yaml` records the settings
at `prepare` time; any later command with a different config for the same
`output_dir` is refused, except for changes to `threads`, `resources`, `lsf` and
`output_dir` itself, which cannot change results. Software versions are pinned by
`environment.yml` and `tools/install_r_source_packages.sh`.

## Run time

Measured on one LSF cluster core (R 4.3, OpenBLAS single-threaded):

| Work | Time |
|---|---|
| One univariate model (enet + blup + susie, 5-fold CV, ~5,000 SNPs) | ~22 s |
| One gene without isoTWAS (M + 1 models) | ~3.5 min for the average gene (M ≈ 8) |
| isoTWAS, all eight methods | ~83 s per isoform, on top of the above |
| `simulate` | seconds per gene |
| All 19,022 genes, nine scenarios | ~1,100 CPU-hours without isoTWAS; ~4,600 with |

Run time is linear in the number of isoforms. With `genes_per_job: 20` and isoTWAS,
jobs took a median of 5 hours and at most about 11.
