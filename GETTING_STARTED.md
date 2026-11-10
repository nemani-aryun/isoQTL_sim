# Getting started with isoQTL_sim

This guide takes you from nothing to your first power table. Every step has the
exact commands to type and what you should see afterwards. Nothing here assumes
you have used the pipeline before.

**Time needed.** About 30 minutes of typing, plus waiting: 15–40 minutes for the
software install, an hour or more for genotype downloads, and 20 minutes to a few
hours for the first 100-gene run depending on your machine.

**Contents**

1. [Before you start](#1-before-you-start)
2. [Get the code](#2-get-the-code)
3. [Install the software](#3-install-the-software) (Docker, conda, or a cluster without internet)
4. [Run the built-in test](#4-run-the-built-in-test)
5. [Get the reference data](#5-get-the-reference-data) (GENCODE and 1000 Genomes)
6. [Write your configuration](#6-write-your-configuration)
7. [Run it](#7-run-it) (on one machine, or on an LSF cluster)
8. [Read the results](#8-read-the-results)
9. [Your next runs](#9-your-next-runs)
10. [When something goes wrong](#10-when-something-goes-wrong)

---

## 1. Before you start

You need **one** of these:

- A Linux machine (x86_64) where you can install software with conda or micromamba.
- Any machine with [Docker](https://docs.docker.com/get-docker/) (Linux, macOS,
  Windows with WSL). macOS without Docker is not supported, because some R packages
  have no conda build for Apple silicon.
- An LSF cluster. If its nodes have no internet access, see
  [docs/CLUSTER_SETUP.md](docs/CLUSTER_SETUP.md) after step 2.

Disk: 5 GB for the software, 2 GB for the annotation, and a few GB for filtered
genotypes (the raw 1000 Genomes download is tens of GB; you can delete it after
filtering). A 100-gene run without reads needs under 1 GB.

**One folder for everything.** Put the code, data, configs and results under a
single directory. This makes the Docker route work without path translation.
Throughout this guide that directory is `$ISOQTL_HOME`:

```bash
export ISOQTL_HOME=$HOME/isoqtl       # choose any location with enough space
mkdir -p "$ISOQTL_HOME"
cd "$ISOQTL_HOME"
```

Put that `export` line in your `~/.bashrc` (or `~/.zshrc`) so it is set in every
new terminal, or type it again each time.

## 2. Get the code

```bash
cd "$ISOQTL_HOME"
git clone https://github.com/nemani-aryun/isoQTL_sim.git
ls isoQTL_sim
```

You should see, among other things, `bin`, `R`, `config`, `tests`, `tools`,
`environment.yml` and `Dockerfile`.

## 3. Install the software

Pick one option. At the end of this step the command `isoqtl-sim` must work in your
terminal; the rest of the guide uses it.

### Option A: Docker (any operating system)

Build the image once (10–30 minutes; it downloads about 2 GB):

```bash
cd "$ISOQTL_HOME/isoQTL_sim"
docker build --platform linux/amd64 -t isoqtl-sim .
```

> The image has not yet been built by the authors (development happened on a
> machine without Docker). If the build fails, please open an issue with the last 30
> lines of output, and use option B or C meanwhile.

Then define `isoqtl-sim` as a shell function that runs the container with your
folder mounted at the same path inside and outside:

```bash
cat >> ~/.bashrc <<'EOF2'
export ISOQTL_HOME=$HOME/isoqtl
isoqtl-sim() {
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$ISOQTL_HOME:$ISOQTL_HOME" -w "$ISOQTL_HOME" isoqtl-sim isoqtl-sim "$@"
}
isoqtl-shell() {
  docker run --rm -it -u "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$ISOQTL_HOME:$ISOQTL_HOME" -w "$ISOQTL_HOME" isoqtl-sim bash "$@"
}
EOF2
source ~/.bashrc
isoqtl-sim --help
```

(Use `~/.zshrc` instead of `~/.bashrc` on macOS.) You should see the usage text
starting with `isoqtl-sim: run the isoQTL_sim pipeline`. `isoqtl-shell` opens a
shell inside the container with all tools (`R`, `salmon`, `tabix`, `plink2`) on
PATH; use it for the data preparation commands in step 5.

### Option B: conda or micromamba (Linux x86_64)

If you don't have micromamba:

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)   # accept the defaults
source ~/.bashrc
```

Create the environment beside the code and add the two packages that have no conda
build (isotwas and its dependencies, compiled from source; 10–20 minutes):

```bash
cd "$ISOQTL_HOME/isoQTL_sim"
micromamba create -y -p ./env -f environment.yml
micromamba run -p ./env bash tools/install_r_source_packages.sh fetch ./r_src
micromamba run -p ./env bash tools/install_r_source_packages.sh install ./r_src
```

With `conda` instead of `micromamba`: `conda env create -p ./env -f environment.yml`
and `conda run --no-capture-output -p ./env bash ...`. The last line printed should
be `isotwas 2.0.0 loads`.

Then, in **every new terminal**, activate it and put the command on your PATH:

```bash
cd "$ISOQTL_HOME/isoQTL_sim"
micromamba activate ./env             # or: conda activate ./env
chmod +x bin/isoqtl-sim
export PATH="$PWD/bin:$PATH"
isoqtl-sim --help
```

Check the tools are there; each prints a version with no error:

```bash
R --version | head -1        # R 4.3.x
salmon --version
tabix --version | head -1
plink2 --version | head -1
```

### Option C: LSF cluster without internet access

Compute nodes often cannot download anything. `tools/make_offline_bundle.sh` runs on
a machine with internet access (a laptop is fine) and produces a folder with every
package; you copy it to the cluster and install from it offline.
[docs/CLUSTER_SETUP.md](docs/CLUSTER_SETUP.md) has the steps, including how to test
and how to submit. Follow its sections 1–3, then come back to step 5 here for the
data.

## 4. Run the built-in test

The test makes a tiny synthetic genome and genotype set, runs every stage of the
pipeline on it, and checks the outputs. It needs no downloads and takes 3–10
minutes.

```bash
cd "$ISOQTL_HOME/isoQTL_sim"
bash tests/run_toy_e2e.sh "$ISOQTL_HOME/toy_test"       # option B or C
isoqtl-shell -c "bash /opt/isoqtl_sim/tests/run_toy_e2e.sh $ISOQTL_HOME/toy_test"   # option A
```

You should see each stage announced (`[simulate] 3 job(s)`, `[reads] 80 job(s)`, …),
a printed summary table, then:

```
Rerun skipped finished work in 1s
Wrote .../toy_test/results/submit_lsf.sh
End-to-end test passed
```

If it ends with `FAIL:` or an R error, go to [section 10](#10-when-something-goes-wrong).
Also run the statistical checks (about a minute, ends with `All checks passed`):

```bash
Rscript tests/test_sim_core.R                                      # option B or C
isoqtl-shell -c "Rscript /opt/isoqtl_sim/tests/test_sim_core.R"    # option A
```

You can delete `$ISOQTL_HOME/toy_test` afterwards.

## 5. Get the reference data

You need an annotation (which genes and isoforms exist) and genotypes (real LD).
Both must be on the same genome build. The commands below use GENCODE v38 and the
1000 Genomes high-coverage data, both on GRCh38, and produce exactly the inputs the
example config expects. With Docker, run them inside `isoqtl-shell` so `plink2` and
`tabix` are available.

### 5a. Annotation (5 minutes, 2 GB)

```bash
mkdir -p "$ISOQTL_HOME/data/gencode_v38" && cd "$ISOQTL_HOME/data/gencode_v38"
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
gunzip gencode.v38.annotation.gtf.gz gencode.v38.transcripts.fa.gz
ls -lh
```

You should see `gencode.v38.annotation.gtf` (about 1.4 GB) and
`gencode.v38.transcripts.fa` (about 380 MB). (Keeping them uncompressed is simplest;
the pipeline itself can read the `.gz` files too.) Any GENCODE release works the same way;
change `38` throughout.

### 5b. Genotypes from 1000 Genomes (1–3 hours, mostly downloading)

The pipeline needs one bgzipped, tabix-indexed VCF per autosome containing
biallelic SNPs for the samples you want. This recipe takes the 3,202-sample
high-coverage release, keeps the unrelated European-ancestry samples (about 525)
and SNPs with minor allele frequency of at least 1%.

```bash
mkdir -p "$ISOQTL_HOME/data/1kg_raw" "$ISOQTL_HOME/data/genotypes"
cd "$ISOQTL_HOME/data/1kg_raw"

# Sample list: European superpopulation, no parents listed (unrelated).
curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
awk 'NR > 1 && $7 == "EUR" && $3 == "0" && $4 == "0" { print $2 }' \
  20130606_g1k_3202_samples_ped_population.txt > eur_unrelated.txt
wc -l eur_unrelated.txt        # about 525

# Download and filter each chromosome. Start with chr22 to check the recipe works.
BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased
for chr in 22; do      # then: for chr in $(seq 1 21); do
  f=CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
  curl -C - -O "$BASE/$f"
  plink2 --vcf "$f" --keep eur_unrelated.txt --snps-only just-acgt --max-alleles 2 \
    --maf 0.01 --export vcf bgz id-paste=iid --out "$ISOQTL_HOME/data/genotypes/eur_chr${chr}"
  tabix -f -p vcf "$ISOQTL_HOME/data/genotypes/eur_chr${chr}.vcf.gz"
done
```

Check the result for chr22:

```bash
tabix -l "$ISOQTL_HOME/data/genotypes/eur_chr22.vcf.gz"                       # prints the chromosome name
tabix -H "$ISOQTL_HOME/data/genotypes/eur_chr22.vcf.gz" | tail -1 | cut -f1-12  # header; sample IDs from column 10
```

Note what `tabix -l` printed: `22` means `chr_prefix: ""` in the config, `chr22`
means `chr_prefix: "chr"`. (plink2 writes `22` by default.) Then run the loop for
chromosomes 1 to 21; the raw files are large, so delete `1kg_raw` when done.

**Using your own genotypes instead.** Any VCFs work if they are bgzipped and
tabix-indexed, one file per autosome with the chromosome number in the name, have
`GT` first in the FORMAT column, and use the same build as the annotation.
Multi-allelic sites and indels are ignored, and missing genotypes are set to the
mean. Name them so that one template with `{chr}` matches them all, e.g.
`/data/geno/cohort_chr{chr}.vcf.gz`.

## 6. Write your configuration

Copy the template and open it in any text editor:

```bash
cd "$ISOQTL_HOME"
cp isoQTL_sim/config/example.yaml my_run.yaml
nano my_run.yaml           # or any editor
```

Change these five lines (and only these, for a first run):

| Line | Set it to |
|---|---|
| `output_dir:` | `/home/you/isoqtl/results/my_run` (your `$ISOQTL_HOME`, spelled out, plus `/results/my_run`) |
| `gtf:` | `/home/you/isoqtl/data/gencode_v38/gencode.v38.annotation.gtf` |
| `fasta:` | `/home/you/isoqtl/data/gencode_v38/gencode.v38.transcripts.fa` |
| `vcf_template:` | `/home/you/isoqtl/data/genotypes/eur_chr{chr}.vcf.gz` (keep `{chr}` literally) |
| `chr_prefix:` | `""` if `tabix -l` printed `22`, `"chr"` if it printed `chr22` |

Write paths out in full; `$ISOQTL_HOME` and `~` are not expanded inside the file.
Leave `n_genes: 100` and `run_reads: false` for the first run. Every other setting is
explained in [docs/CONFIG.md](docs/CONFIG.md).

Now let the pipeline read the annotation, pick 100 genes at random and list the
samples. This also checks that the files can be opened (1–3 minutes):

```bash
isoqtl-sim prepare --config my_run.yaml
```

You should see `Prepared 100 genes, N transcripts, 525 samples` (your numbers may
differ slightly). The chosen genes are in `results/my_run/params/genes.tsv`. Then:

```bash
isoqtl-sim plan --config my_run.yaml
```

```
simulate                           jobs=7      mem=8G hours=4 cores=1
train_truth                        jobs=7      mem=16G hours=12 cores=1
summarize                          jobs=1      mem=16G hours=1 cores=1
```

Seven jobs because 100 genes at 15 per job make 7 chunks. Memory and hours are the
requests used for cluster jobs; on one machine they are informational.

## 7. Run it

### On one machine

```bash
isoqtl-sim run --config my_run.yaml -j 4
```

`-j` is how many jobs run at once; use the number of cores you can spare. For 100
genes, expect 20 minutes to 2 hours depending on the machine (each gene trains one
model per isoform plus one for the gene, at 5–25 seconds per model). The output
looks like:

```
[simulate] 7 job(s)
[train_truth] 7 job(s)
[summarize] 1 job(s)
Pipeline complete.
```

Progress while it runs, from another terminal:

```bash
isoqtl-sim status --config my_run.yaml       # e.g. "train_truth   3/7"
ls results/my_run/logs/train_truth/           # one log per job
```

If a job fails, `run` prints `FAILED train_truth[3]: see .../logs/train_truth/3.log`
and stops after that stage. Fix the cause (the log says what went wrong) and run the
same command again; finished jobs are skipped.

### On an LSF cluster

`submit` writes a script with one `bsub` job array per stage, each waiting for the
previous one:

```bash
isoqtl-sim submit --config my_run.yaml
cat results/my_run/submit_lsf.sh              # read it once
bash results/my_run/submit_lsf.sh             # submit
bjobs -w                                      # your jobs
isoqtl-sim status --config my_run.yaml        # progress per stage
```

Before the first submission, set the `lsf:` block of the config: your queue names
(`bqueues` lists them), and `container_exec`, the command that puts each job inside
the software environment (for option B: `conda run --no-capture-output -p
/full/path/to/isoQTL_sim/env`). Jobs that fail leave later stages pending; read the
log in `results/my_run/logs/<stage>/<index>.log`, kill the pending jobs with
`bkill -J "isoqtl_*"`, fix the cause, and run the `bash` line again. See
[docs/CLUSTER_SETUP.md](docs/CLUSTER_SETUP.md) for the details.

## 8. Read the results

```bash
column -t -s $'\t' results/my_run/summary/power_summary.tsv
```

One row per `source` (where the expression came from: `truth`, or
`salmon_gencode_v38` after a read-level run), `level` (the kind of test) and
`scenario` (the simulated trait). The columns:

| Column | Meaning |
|---|---|
| `h2ge`, `n_gwas`, `mediator` | The scenario: trait variance explained through this gene, GWAS sample size, and whether the trait acts through total gene expression or through one isoform |
| `n_genes_simulated` | Genes with simulated expression (the denominator for the rates) |
| `n_genes_not_simulated` | Genes skipped because their window had too few SNPs |
| `n_genes_modeled` | Genes with at least one usable model at this level; the rest count as misses |
| `rejection_rate` | Share of simulated genes with p < `alpha` (0.05). In the `null` row this is the false-positive rate and should be near 0.05; elsewhere it is power |
| `rejection_rate_bonferroni` | Same, with p < `alpha` / `n_genes_simulated` |
| `mean_cv_r2` | Mean cross-validated R² of the models used at this level |

The `level` values:

- `gene`: one model of total gene expression, one test.
- `isoform_acat`: one model per isoform, trained separately, and the isoform tests
  combined into one gene-level p-value with ACAT.
- `isotwas_acat`: as above, but the isoforms are modeled jointly with isoTWAS (only
  if `models.isotwas` is set in the config).

For single-isoform genes the isoform model is the gene model, so all levels cover
the same genes. What to look for: the `null` rows near 0.05; how `gene` compares with
`isoform_acat` when the mediator is `gene` versus `isoform`; and how power grows with
`h2ge`.

Per-gene detail is in `results/my_run/summary/gene_tests.tsv.gz` (one row per
source, level, scenario and gene, with its p-value, CV R² and number of isoforms),
and every trained model's test is in `results/my_run/twas/<source>/chunk_*.tsv.gz`.

## 9. Your next runs

**Each run gets its own `output_dir`.** The pipeline refuses to reuse a directory
with a changed config (`These settings changed since prepare ran`), because finished
jobs would silently keep the old settings. Copy the config, change `output_dir`,
then edit what you like. Changing `threads`, `resources` or the `lsf` block is
always allowed.

Common next steps, each a few lines in the config
([docs/CONFIG.md](docs/CONFIG.md) has the details):

- **More scenarios**: add rows under `trait.scenarios` with other `h2ge`, `n_gwas`
  or `mediator` values. Each adds one GWAS per gene and one test per model; cheap.
- **isoTWAS models**: uncomment the `isotwas:` block under `models`. Adds about 80
  seconds per isoform per gene on one core.
- **Read-level comparison**: set `run_reads: true`. This simulates a FASTQ pair per
  sample, quantifies with Salmon, and adds `salmon_<annotation>` rows to the summary.
  Set `reads.library_size` to the read depth you want *for the simulated genes
  only*; the pipeline expresses nothing else, so 100 genes at 150,000 read pairs
  matches the per-gene depth of 30 million over the whole transcriptome.
- **Another annotation for quantification**: add a second entry under `reference`
  (e.g. `gencode_v45`) and list it under `quant.salmon`; `truth_annotation` stays
  v38. This measures what an annotation mismatch costs.
- **All genes**: `n_genes: all` gives about 19,000 genes. Raise `genes_per_job` to
  20–25 so the job arrays stay under LSF limits. Measured on one cluster core: about
  22 seconds per univariate model (a gene with M isoforms trains M + 1), plus 83
  seconds per isoform for isoTWAS; roughly 1,100 CPU-hours without isoTWAS and
  4,600 with. Check `du -sh results/my_run` after a pilot and scale by 190 for disk.

To recompute only the summary after changing `summary.alpha`, or after a code update
that changes the summary step: `Rscript R/run_step.R summarize my_run.yaml 1`
(inside the environment or `isoqtl-shell`).

## 10. When something goes wrong

| What you see | Cause | What to do |
|---|---|---|
| `Parameters not found. Run isoqtl-sim prepare` | `plan`, `run` or `submit` before `prepare` | Run `prepare` |
| `tabix: command not found` or `Rscript: command not found` | Environment not active | Option B: `micromamba activate ./env`; option A: use the `isoqtl-sim` function |
| `Prepared 0 genes` | GTF has no `protein_coding` autosomal genes, or `fasta` doesn't match the GTF | Use the GTF and transcript FASTA of the same GENCODE release |
| Every gene has status `no_snps` (see `results/*/simulated_status/`) | Wrong `chr_prefix`, or VCF naming doesn't match `vcf_template` | Compare `tabix -l` with `chr_prefix`; check `ls` of the template with `{chr}` replaced by `22` |
| `Samples missing from VCF` | `sample_file` has IDs not in the VCF header | Fix the sample file, or remove the setting |
| `FAILED simulate[3]: see .../logs/simulate/3.log` | Look in that log | The last lines say why; usually a path or an out-of-memory kill |
| `These settings changed since prepare ran` | Config edited after `prepare` | Use a new `output_dir` (section 9) |
| `isotwas failed for ENSG...` lines in a train log | One gene's joint models failed | The gene keeps its univariate models; a few of these are normal. Many means a package problem: run `Rscript -e 'library(isotwas)'` |
| LSF: `TERM_MEMLIMIT` / `TERM_RUNLIMIT` in a log | Not enough memory or time | Add `resources: {train: {mem_gb: 32, hours: 24}}` (stage name as in `plan`) and submit again |
| LSF: `needs N array jobs, above the limit` | Too many jobs for one array | Follow the message: raise `genes_per_job`, then `prepare` into a new `output_dir` |

If none of these fit, open an issue with the config (paths removed), the command,
and the last 30 lines of the failing log.
