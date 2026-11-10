# isoQTL_sim

**Power simulations for TWAS and isoform-level TWAS, from real genotypes.**

isoQTL_sim answers a design question: *if a trait is affected by a gene's expression,
or by one of its isoforms, how often will a gene-level TWAS, an isoform-level TWAS or
an isoTWAS find it, at a given GWAS size and effect size?* It simulates isoform
expression for real genes from real 1000 Genomes genotypes, simulates GWAS summary
statistics for traits that act through those genes, trains prediction models the way
TWAS tools do, tests them, and reports power and calibration in one table.

Optionally it also simulates RNA-seq reads from the expression, quantifies them with
Salmon, and repeats the whole comparison on the estimated expression, so you can see
what quantification costs.

Everything is driven by one configuration file and one command. Nothing in the code
is specific to a lab, a cluster or a path.

## What you get

A table like this, one row per expression source, test level and trait scenario
(from a 100-gene pilot run, GWAS of 100,000 people):

```
source  level         scenario      h2ge    mediator  n_genes_simulated  n_genes_modeled  rejection_rate  rejection_rate_bonferroni  mean_cv_r2
truth   gene          null          0       gene      100                99               0.02            0                          0.095
truth   gene          gene_5e-4     0.0005  gene      100                99               0.89            0.81                       0.095
truth   gene          isoform_5e-4  0.0005  isoform   100                99               0.84            0.72                       0.095
truth   isoform_acat  null          0       gene      100                100              0.06            0                          0.064
truth   isoform_acat  gene_5e-4     0.0005  gene      100                100              0.96            0.89                       0.064
truth   isoform_acat  isoform_5e-4  0.0005  isoform   100                100              0.96            0.89                       0.064
truth   isotwas_acat  isoform_5e-4  0.0005  isoform   100                100              0.95            0.93                       0.069
```

`rejection_rate` in the `null` row is the false-positive rate (it should be near
0.05); in the other rows it is power. `level` is the kind of test: `gene` (one
model of total expression), `isoform_acat` (one model per isoform, combined with
ACAT) or `isotwas_acat` (isoforms modeled jointly with
[isotwas](https://github.com/bhattacharya-a-bt/isotwas), combined with ACAT).
A per-gene table (`gene_tests.tsv.gz`) sits beside it.

## How it works

1. **prepare**: pick genes from a GENCODE annotation and samples from your VCFs.
2. **simulate**: for each gene, draw causal SNPs in the real genotypes, simulate
   isoform expression with shared and isoform-specific effects, then simulate GWAS
   Z-scores for each trait scenario using the real LD.
3. **reads → salmon → aggregate** (optional): simulate FASTQ files with Rsubread,
   quantify with Salmon, normalize, and build gene and isoform expression matrices.
4. **train**: fit elastic net, BLUP and SuSiE (and optionally isoTWAS) models for
   each gene and isoform, keep the best by cross-validation, and test each against
   every scenario's GWAS.
5. **summarize**: power, calibration and prediction accuracy per level and scenario.

Each stage is a job array that skips work already finished, so a run can be
stopped and resumed, and runs on one machine (`isoqtl-sim run`) or an LSF cluster
(`isoqtl-sim submit`). [PIPELINE_OVERVIEW.md](PIPELINE_OVERVIEW.md) has the details
and the statistical model.

## What you need

- **A Linux x86_64 machine**, or any machine with Docker. A laptop is enough for a
  100-gene run; all 19,000 genes need a cluster (about 1,100 CPU-hours without
  isoTWAS, 4,600 with).
- **GENCODE annotation**: the GTF and the transcript FASTA of one release.
- **Genotypes**: bgzipped, tabix-indexed VCFs, one per autosome, on the same genome
  build as the annotation. The guide shows how to make them from 1000 Genomes.

## Install and run in 15 minutes

[GETTING_STARTED.md](GETTING_STARTED.md) walks through every step with the exact
commands and what you should see after each. In short:

```bash
git clone https://github.com/nemani-aryun/isoQTL_sim.git
cd isoQTL_sim

# Software, option A: Docker
docker build --platform linux/amd64 -t isoqtl-sim .

# Software, option B: conda/micromamba on Linux
micromamba create -y -p ./env -f environment.yml
micromamba run -p ./env bash tools/install_r_source_packages.sh fetch ./r_src
micromamba run -p ./env bash tools/install_r_source_packages.sh install ./r_src

# Check everything works on built-in synthetic data (a few minutes)
bash tests/run_toy_e2e.sh /tmp/isoqtl_toy        # ends with "End-to-end test passed"

# Your own run
cp config/example.yaml my_run.yaml                # edit the five /path/to/... lines
isoqtl-sim prepare --config my_run.yaml
isoqtl-sim plan    --config my_run.yaml           # stages and job counts
isoqtl-sim run     --config my_run.yaml -j 8      # this machine, 8 jobs at once
isoqtl-sim status  --config my_run.yaml
column -t -s $'\t' results/my_run/summary/power_summary.tsv
```

On an LSF cluster, replace `run` with `submit`, which writes one `bsub` job array per
stage and chains them. For a cluster without internet access,
`tools/make_offline_bundle.sh` builds an installable package on a machine that has
it; [docs/CLUSTER_SETUP.md](docs/CLUSTER_SETUP.md) is a worked example of that route.

## Documentation

| Document | Read it when |
|---|---|
| [GETTING_STARTED.md](GETTING_STARTED.md) | You are running isoQTL_sim for the first time |
| [docs/CONFIG.md](docs/CONFIG.md) | You want to know what a setting does, or change the simulation |
| [PIPELINE_OVERVIEW.md](PIPELINE_OVERVIEW.md) | You want to know what each stage computes, the model, and the output files |
| [docs/CLUSTER_SETUP.md](docs/CLUSTER_SETUP.md) | You are installing on an LSF cluster, especially one without internet access |

## What has been tested

- Simulation, univariate TWAS and isoTWAS on true expression: all 19,022
  protein-coding autosomal genes (GENCODE v38, 500 European 1000 Genomes samples)
  on an LSF cluster.
- Read simulation and Salmon quantification: 100 genes on the same cluster.
- Every step on synthetic data: `tests/run_toy_e2e.sh`, run on macOS and Linux.
- **Not yet tested**: STAR/featureCounts quantification, FastQC/MultiQC, and the
  Docker image build. These are wired in but have not run anywhere; please report
  what you find.

## Citation

```bibtex
@article{nemani2026,
  title={A realistic simulation tool for powering eQTL studies and TWAS},
  author={Nemani, Aryun and Head, Taylor and Bhattacharya, Arjun},
  journal={bioRxiv},
  year={2026}
}
```
