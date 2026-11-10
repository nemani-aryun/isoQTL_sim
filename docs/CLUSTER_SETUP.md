# Running on an LSF cluster

`isoqtl-sim submit` turns a run into one LSF job array per stage, chained so each
stage starts when the previous one has finished. This guide covers what is specific
to a cluster: installing the software where compute nodes may have no internet
access, the `lsf` block of the config, submitting, watching, recovering from failed
jobs, and scaling to every gene. Read [GETTING_STARTED.md](GETTING_STARTED.md)
first for the data and the config; this page assumes you have both.

Replace `CONDA_BIN` below with your cluster's conda or mamba executable
(`which conda`, or `module avail conda` / `module avail miniforge` and the module's
path), and set `ISOQTL_HOME` to a directory on fast shared storage with enough quota
(usually a scratch filesystem, not your home directory):

```bash
export ISOQTL_HOME=/scratch/$USER/isoqtl
export CONDA_BIN=/path/to/conda
mkdir -p "$ISOQTL_HOME" && cd "$ISOQTL_HOME"
git clone https://github.com/nemani-aryun/isoQTL_sim.git
```

Put the two `export` lines in your `~/.bashrc`, or type them in every new login.

## 1. Install the software environment

Do the install inside an interactive job, since login nodes often limit CPU and
memory. Conda caches packages in your home directory by default, which can exceed a
quota, so point the cache at `ISOQTL_HOME` too:

```bash
bsub -Is -q short -n 4 -W 2:00 -M 16 -R "rusage[mem=16]" bash
export ISOQTL_HOME=/scratch/$USER/isoqtl CONDA_BIN=/path/to/conda   # if not carried over
export CONDA_PKGS_DIRS=$ISOQTL_HOME/conda_pkgs
cd "$ISOQTL_HOME/isoQTL_sim"
```

(Queue names and limits vary; `bqueues` lists yours. `-M` is memory in GB on some
clusters and MB on others; `bqueues -l short` shows the units with the limits.)

### If compute nodes can reach the internet

```bash
"$CONDA_BIN" env create -p "$ISOQTL_HOME/env" -f environment.yml
"$CONDA_BIN" run --no-capture-output -p "$ISOQTL_HOME/env" bash tools/install_r_source_packages.sh fetch "$ISOQTL_HOME/r_src"
"$CONDA_BIN" run --no-capture-output -p "$ISOQTL_HOME/env" bash tools/install_r_source_packages.sh install "$ISOQTL_HOME/r_src"
```

The last line printed should be `isotwas 2.0.0 loads`.

### If compute nodes cannot reach the internet

Build an **offline bundle** on any machine with internet access (a laptop is fine,
macOS or Linux), copy it to the cluster, and install from it there. The bundle is
about 1 GB and contains every package at a fixed version, so the environment is the
same on every cluster.

**On the machine with internet access:**

```bash
git clone https://github.com/nemani-aryun/isoQTL_sim.git && cd isoQTL_sim
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)      # micromamba, if you don't have it
bash tools/make_offline_bundle.sh "$(which micromamba)" ./isoqtl_bundle
```

It prints `478 packages downloaded` (or similar) and ends with the size of the
bundle. Copy the `isoqtl_bundle` folder to `$ISOQTL_HOME/isoqtl_bundle` on the
cluster with whatever transfer tool your cluster allows (`scp -r`, `rsync -a`, or an
SFTP client such as FileZilla in binary mode). The installer checks every file, so
an incomplete transfer is caught.

**On the cluster, in the interactive job:**

```bash
bash "$ISOQTL_HOME/isoqtl_bundle/install_offline.sh" "$CONDA_BIN" "$ISOQTL_HOME/env"
```

It prints `Bundle files verified` first (or names files that arrived incomplete:
copy those again and rerun), installs the conda packages, then compiles seven R
packages from source, and ends with `isotwas 2.0.0 loads`. Allow 15–30 minutes.

### Check the environment

```bash
eval "$("$CONDA_BIN" shell.bash hook)"
conda activate "$ISOQTL_HOME/env"
R --version | head -1           # R 4.3.x
salmon --version
tabix --version | head -1
Rscript -e 'for (p in c("data.table","glmnet","rrBLUP","susieR","Rsubread","tximport","isotwas")) cat(p, as.character(packageVersion(p)), "\n")'
```

## 2. Run the tests

Still in the interactive job with the environment active:

```bash
cd "$ISOQTL_HOME/isoQTL_sim"
chmod +x bin/isoqtl-sim
Rscript tests/test_sim_core.R                           # ends with "All checks passed"
bash tests/run_toy_e2e.sh "$ISOQTL_HOME/toy_test"       # ends with "End-to-end test passed"
exit                                                    # leave the interactive job
```

The second test takes a few minutes and runs every stage on synthetic data. Repeat
both after updating the code, giving the second a new directory (it skips work it
finds finished).

## 3. Set up each login session

These are lost when you log out; your jobs and results are not.

```bash
export ISOQTL_HOME=/scratch/$USER/isoqtl CONDA_BIN=/path/to/conda
eval "$("$CONDA_BIN" shell.bash hook)"
conda activate "$ISOQTL_HOME/env"
cd "$ISOQTL_HOME/isoQTL_sim"
export PATH="$PWD/bin:$PATH"
```

## 4. Configure the `lsf` block

Start from `config/example.yaml` (or `config/example_full_lsf.yaml` for all genes,
`config/example_reads_lsf.yaml` for a read-level run), fill in the data paths as in
GETTING_STARTED, and set these ([CONFIG.md](CONFIG.md) has every field):

```yaml
lsf:
  job_prefix: myrun            # job names start with this; one prefix per run
  queue_short: short           # from `bqueues`; for stages of up to short_hours
  queue_long: medium           # must allow the longest stage (see `plan`)
  short_hours: 3
  max_running: 200             # how many array jobs may run at once
  max_array_size: 1000         # bparams -a | grep MAX_JOB_ARRAY_SIZE
  mem_per_gb: 1                # 1024 if your cluster's bsub -M counts MB
  container_exec: /path/to/conda run --no-capture-output -p /scratch/USER/isoqtl/env
```

`container_exec` is what makes each job run inside the environment without an
activation step on the node; write both paths in full. With a container runtime
instead, use `apptainer exec --bind /scratch /path/to/isoqtl-sim.sif` (see section 8).

## 5. Prepare, submit, watch

`prepare` reads the whole GTF, so run it in an interactive job:

```bash
bsub -Is -q short -n 1 -W 1:00 -M 16 -R "rusage[mem=16]" isoqtl-sim prepare --config my_run.yaml
isoqtl-sim plan   --config my_run.yaml      # stages, job counts, memory and hours per job
isoqtl-sim submit --config my_run.yaml      # writes <output_dir>/submit_lsf.sh
cat  <output_dir>/submit_lsf.sh             # read it once: one bsub per stage
bash <output_dir>/submit_lsf.sh
```

Then:

```bash
bjobs -w                                    # your jobs; later stages show PEND until the previous one finishes
isoqtl-sim status --config my_run.yaml      # finished jobs per stage, e.g. "train_truth  120/952"
ls <output_dir>/logs/train_truth/           # one log per array job
grep -lE "Exited with|TERM_" <output_dir>/logs/*/*.log    # logs of failed jobs
```

**When a job fails,** the stages after it stay pending forever (they wait for
`done(...)`). Read the failed log, fix the cause, remove the pending jobs with
`bkill -J "myrun_*"`, and run the `bash <output_dir>/submit_lsf.sh` line again.
Finished jobs have a marker in `<output_dir>/.done/` and are skipped; only the failed
and pending work runs.

## 6. Scaling to every gene

With `n_genes: all` (about 19,000 genes for GENCODE v38):

- **Array size.** `genes_per_job: 20` gives about 950 jobs per stage, under the usual
  limit of 1,000. `submit` refuses a stage that exceeds `max_array_size` and says what
  to change.
- **Time per job.** Measured on one cluster core: about 22 s per univariate model
  (a gene with M isoforms trains M + 1 models) and about 83 s per isoform for isoTWAS
  with all eight methods. At 20 genes per job that is a median of roughly 1.5 hours
  without isoTWAS and 5 hours with, with the longest jobs at about 3 and 11 hours.
  Run a 100-gene pilot first, read `grep -h "Run time" <output_dir>/logs/train_truth/*.log`,
  and scale by `genes_per_job / 15`.
- **Queue limit.** The default request for `train` stages is 12 hours. Check
  `bqueues -l <queue_long> | grep -A2 RUNLIMIT`, and raise or lower the request with
  `resources: {train: {hours: 24}}` in the config.
- **Disk.** Check `du -sh <output_dir>` after the pilot and multiply by
  `n_genes / 100`. Read-level runs add the FASTQ files (two per sample) and Salmon
  outputs; those dominate.
- **Total compute.** About 1,100 CPU-hours for all genes without isoTWAS and about
  4,600 with, before any read simulation.

## 7. Recomputing only the summary

If you change `summary.alpha`, or update the code in a way that affects only the
summary step, keep the finished jobs and rerun the last step alone:

```bash
cp -r <output_dir>/summary <output_dir>/summary_old
bsub -Is -q short -n 1 -W 1:00 -M 16 -R "rusage[mem=16]" Rscript R/run_step.R summarize my_run.yaml 1
```

## 8. Using a container instead of conda

If your cluster has Apptainer or Singularity (`module avail apptainer`) and you
have built the Docker image (`docker build --platform linux/amd64 -t isoqtl-sim .`)
and pushed it to a registry, convert it once:

```bash
apptainer pull isoqtl-sim.sif docker://REGISTRY/isoqtl-sim:latest
```

and set `container_exec: apptainer exec --bind /scratch,/data /path/to/isoqtl-sim.sif`
with every filesystem your config refers to in `--bind`. No conda environment is
needed then. The image build itself has not yet been exercised by the authors.

## Troubleshooting

| What you see | Likely cause | What to do |
|---|---|---|
| `conda: command not found` inside a job log | `container_exec` has a relative or wrong conda path | Use the full path printed by `which conda` on the login node |
| `Disk quota exceeded` during install | Conda cache in your home directory | `export CONDA_PKGS_DIRS=$ISOQTL_HOME/conda_pkgs` and retry |
| `Failed to connect`, `CondaHTTPError`, or `post-link script failed` during install | The node has no internet access | Use the offline bundle (section 1). A bundle built by the current `make_offline_bundle.sh` needs no downloads at all |
| `These bundle files are missing or incomplete` | Transfer dropped or truncated files | Copy the listed files again, rerun the installer |
| Job killed with `TERM_MEMLIMIT` | Not enough memory requested | Add `resources: {<stage>: {mem_gb: 32}}` and submit again |
| Job killed with `TERM_RUNLIMIT` | Not enough time, or the queue's limit is lower than the request | Raise `hours` in `resources`, or use a queue with a longer limit |
| Jobs pend with "memory" or "resource requirement" errors | LSF counts memory in MB on this cluster | `lsf: {mem_per_gb: 1024}` |
| `needs N array jobs, above the limit` | More jobs than one array may hold | Raise `genes_per_job` and run `prepare` into a new `output_dir` |
| `These settings changed since prepare ran` | Config edited after `prepare` | Use a new `output_dir`; `lsf`, `threads` and `resources` changes are always allowed |
| Every gene has status `no_snps` | Wrong `chr_prefix` or `vcf_template` | Compare `tabix -l` output with `chr_prefix` |
| `Samples missing from VCF` | `sample_file` IDs not in the VCF header | Fix the sample file or remove the setting |
| Later stages `PEND` forever | An earlier job failed | Section 5: read the log, `bkill -J "<prefix>_*"`, fix, resubmit |
