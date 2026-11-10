# Configuration reference

A run is described by one YAML file. `config/example.yaml` is a complete template
with the defaults; copy it and change what you need. This page explains every
setting. "Locked" means the setting is fixed once `isoqtl-sim prepare` has run for
an `output_dir`; to change it, use a new `output_dir`. Everything else can change
between commands.

Paths must be absolute and written out in full (`~` and `$HOME` are not expanded).

## Top level

| Setting | Locked | Meaning |
|---|---|---|
| `output_dir` | no | Where everything for this run is written. One directory per run. |
| `seed` | yes | Integer. Same seed and settings give the same genes, samples, expression, GWAS and models. |
| `genes_per_job` | yes | Genes handled by one `simulate` or `train` job. More genes per job means fewer, longer jobs. 15 gives 7 jobs for 100 genes; use 20–25 for all genes so arrays stay under LSF's limit (usually 1,000). |
| `threads` | no | Cores given to Salmon and STAR jobs. |
| `run_reads` | yes | `false`: simulate expression and test models on it. `true`: also simulate reads, quantify them, and test models on the estimates. Much more compute and disk. |
| `truth_annotation` | yes | Name of the entry under `reference` whose transcripts are the true isoforms. |

## `reference`

A named list of annotations. The one named by `truth_annotation` defines the
simulation; others can be added for quantification (`quant.salmon`,
`quant.star.annotations`) to measure annotation mismatch.

```yaml
reference:
  gencode_v38:
    gtf: /data/gencode_v38/gencode.v38.annotation.gtf
    fasta: /data/gencode_v38/gencode.v38.transcripts.fa
    salmon_index: /data/salmon_index_v38     # optional
```

| Setting | Meaning |
|---|---|
| `gtf` | GENCODE-style GTF. Genes are taken from its `gene` records with `gene_type "protein_coding"` on chromosomes 1–22. |
| `fasta` | The matching transcript FASTA (GENCODE `transcripts.fa`). Its headers give each transcript's gene; this defines the isoforms of each gene. Also used by read simulation and Salmon. |
| `salmon_index` | Optional. An existing Salmon index for this FASTA. Without it, the pipeline builds one in `output_dir/salmon_index/<name>`. |

## `genotypes`

```yaml
genotypes:
  vcf_template: /data/genotypes/eur_chr{chr}.vcf.gz
  chr_prefix: ""
  sample_file: /data/samples.txt     # optional
  n_samples: 500                     # optional
```

| Setting | Locked | Meaning |
|---|---|---|
| `vcf_template` | yes | Path with `{chr}`, replaced by 1…22. Each file must be bgzipped and have a `.tbi` index (`tabix -p vcf`). Only biallelic SNPs with single-base alleles are used; `GT` must be the first FORMAT field; missing genotypes become the mean. |
| `chr_prefix` | yes | How the VCF names chromosomes: `""` for `22`, `"chr"` for `chr22`. Check with `tabix -l file.vcf.gz`. Wrong value → every gene gets status `no_snps`. |
| `sample_file` | yes | Optional. One sample ID per line; only these (and present in the VCF) are used. Default: all samples in the VCF header. |
| `n_samples` | yes | Optional. A random subset of that size. |

The samples are the reference panel for everything: expression is simulated for
them, models are trained on them, and their LD generates the GWAS statistics.
500 is the usual size; fewer samples make prediction harder (lower CV R²), as in a
real small reference panel.

## `windows`

| Setting | Default | Meaning |
|---|---|---|
| `causal` | 250000 | Causal SNPs are drawn within this many bp of the gene body. |
| `train` | 1000000 | Models are trained on, and GWAS statistics generated for, SNPs within this distance. Must be ≥ `causal`. Larger windows mean more SNPs per model and slower training. |

## `simulation`

| Setting | Default | Meaning |
|---|---|---|
| `n_genes` | 100 | Number of genes to sample at random from the annotation, or `all` (about 19,000 for GENCODE v38). |
| `min_isoforms` | 1 | Optional. Keep only genes with at least this many isoforms (e.g. 2 to study multi-isoform genes only). |
| `n_causal` | 4 | Causal SNPs per isoform. |
| `prop_shared` | [0.5, 0.75, 1] | Share of each isoform's causal SNPs that are shared by all isoforms of the gene; drawn per gene from this list. `1` means every isoform has the same causal SNPs; `0.5` with `n_causal: 4` means 2 shared and 2 isoform-specific. |
| `h2_g` | {beta: [2, 20]} | Cis-heritability of each isoform's expression. Either `{beta: [a, b]}` to draw per gene from a Beta(a, b) distribution (Beta(2, 20) has mean 0.09 and 95% of genes between 0.01 and 0.25), or a single number such as `0.1`. |
| `cov_B` | [0.7, 0.9] | Range for the correlation of shared SNPs' effects across isoforms. High values mean shared SNPs move all isoforms together. |
| `cov_E` | [0, 0.5] | Range for the correlation of residual (non-genetic) noise across isoforms. |

## `trait.scenarios`

A list; each entry simulates one GWAS per gene and adds a column block to the
summary.

```yaml
trait:
  scenarios:
    - {name: "null",       h2ge: 0,      n_gwas: 100000, mediator: gene}
    - {name: gene_5e-4,    h2ge: 0.0005, n_gwas: 100000, mediator: gene}
    - {name: isoform_5e-4, h2ge: 0.0005, n_gwas: 100000, mediator: isoform}
```

| Field | Meaning |
|---|---|
| `name` | Label in the results. Must be unique. Quote `"null"`: unquoted, YAML reads it as "no value". |
| `h2ge` | Share of trait variance explained through this gene. `0` is a null trait: its `rejection_rate` is the false-positive rate. The expected TWAS Z is about `sqrt(n_gwas × h2ge)` times the model's correlation with true expression, so with `n_gwas: 100000`, `h2ge` from 0.0001 (Z ≈ 3 for a perfect model) to 0.001 (Z ≈ 10) spans low to high power. |
| `n_gwas` | GWAS sample size. Only the product `n_gwas × h2ge` matters for power, so vary one of them. |
| `mediator` | `gene`: the trait acts through total gene expression (the sum of the isoforms' genetic values). `isoform`: through one isoform, chosen at random per gene and recorded as `mediator_isoform` in the simulated data. |

Scenarios are cheap: each adds one GWAS draw per gene and one test per trained
model, no extra training.

## `reads` (only with `run_reads: true`)

| Setting | Default | Meaning |
|---|---|---|
| `read_length` | 75 | Bases per read. |
| `library_size` | 30000000 | Read pairs per sample **over the simulated genes only**, since nothing else is expressed. 30 million is the depth of a typical whole-transcriptome sample; for 100 genes out of 19,000, use about 150,000 to get the same depth per gene. |
| `paired_end` | true | Paired or single-end reads. |

## `quant` (only with `run_reads: true`)

| Setting | Default | Meaning |
|---|---|---|
| `salmon` | [gencode_v38] | Annotations (names from `reference`) to quantify with Salmon. Each gives a `salmon_<name>` source in the results. |
| `min_prop_expressed` | 0.75 | Keep genes and isoforms with TPM > 0.1 in at least this share of samples. A gene filtered out at a level cannot be tested there and counts as a miss. |
| `fastqc` | false | Run FastQC per sample and MultiQC. **Not yet tested.** |
| `star` | (absent) | `{index: /path/to/star_index, annotations: [gencode_v38]}` adds STAR alignment and featureCounts gene counts as `star_<name>` sources. **Not yet tested anywhere.** |

## `models`

```yaml
models:
  methods: [enet, blup, susie]
  isotwas:
    methods: [multi_enet, mrce_lasso, joinet, spls, sgl, mtlasso, stacking, univariate]
    omega_est: mean
```

| Setting | Meaning |
|---|---|
| `methods` | Univariate methods fitted to each gene and each isoform, with 5-fold cross-validation; the best CV R² wins. `enet` (elastic net, glmnet with α = 0.5), `lasso` (glmnet with α = 1), `blup` (rrBLUP), `susie` (susieR). Fewer methods means proportionally faster training. |
| `isotwas` | Optional. If present, the isoforms of each multi-isoform gene are also modeled jointly with `isotwas::compute_isotwas`, and the results appear as the `isotwas_acat` level. `methods` are isotwas method names; `univariate` is isotwas's own per-isoform fit. `omega_est` (`mean` or `replicates`) is how isotwas estimates the residual correlation; without replicate measurements use `mean`. Adds about 80 s per isoform per gene with all eight methods; `joinet` and `mrce_lasso` are the slowest. |

## `summary`

| Setting | Default | Meaning |
|---|---|---|
| `alpha` | 0.05 | Significance level for `rejection_rate`; `rejection_rate_bonferroni` uses `alpha / n_genes_simulated`. |

## `resources` (optional; LSF only)

Overrides the memory and time requested per stage, keyed by the stage name as shown
by `isoqtl-sim plan`, or by its first word to cover several stages
(`train` matches `train_truth` and `train_salmon_gencode_v38`).

```yaml
resources:
  train: {hours: 24, mem_gb: 32}
  salmon_gencode_v38: {hours: 8}
```

Only `hours` and `mem_gb` are accepted. The defaults are in the `plan` output.

## `lsf` (only for `isoqtl-sim submit`)

| Setting | Default | Meaning |
|---|---|---|
| `job_prefix` | isoqtl | Job names are `<prefix>_<stage>`. Use a different prefix per run so `bkill -J "<prefix>_*"` targets one run. |
| `queue_short`, `queue_long`, `short_hours` | short, medium, 3 | Stages needing up to `short_hours` go to `queue_short`, others to `queue_long`. `bqueues` lists your cluster's queues and their limits. |
| `max_running` | 200 | Cap on simultaneously running jobs per array (`%N` in the array spec). |
| `max_array_size` | 1000 | Your cluster's limit on jobs per array (`bparams -a | grep MAX_JOB_ARRAY_SIZE`). `submit` refuses stages above it and says how to fix the config. |
| `mem_per_gb` | 1 | Set to `1024` on clusters where `bsub -M` counts megabytes. |
| `container_exec` | "" | Command that runs each job inside the software environment, e.g. `conda run --no-capture-output -p /path/to/env` or `apptainer exec --bind /data /path/to/isoqtl-sim.sif`. Leave empty if jobs inherit an activated environment. |
