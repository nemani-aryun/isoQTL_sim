# Setting up isoQTL_sim on the LSF cluster

This guide creates a conda environment on the cluster, tests it, then runs a
100-gene pilot (steps 1 to 7), the full run on simulated expression (step 9) and a
read-level pilot (step 10). Allow about an hour for steps 1 to 5, most of it waiting
for installs and the queue.

**What's been tested so far.**

- The 100-gene pilot on simulated expression ran on the cluster (steps 1 to 7).
- Read simulation, Salmon and the TWAS on Salmon estimates have run on a Mac only:
  on synthetic data, and on the real GENCODE v38 reference with synthetic genotypes.
  Step 10 is their first run on the cluster.
- isoTWAS model training (step 11) has run on a Mac only, on synthetic data and on the
  real GENCODE v38 reference with synthetic genotypes. Its environment has been solved
  for Linux but never installed there, so step 11 is the first test of both.
- STAR and featureCounts haven't been run anywhere yet, and no config here uses them.

**Paths and names to check.** The conda path, queue names and data paths come from the
lab's existing scripts. If one fails, step 2 shows how to look up the right value.

## 1. Copy the code to the cluster

The new code is on a local branch that hasn't been pushed. Run this **on your Mac**,
replacing `USER` and `CLUSTER_LOGIN_HOST`:

```bash
cd ~/Library/CloudStorage/OneDrive-InsideMDAnderson/isoQTL_sim
ssh USER@CLUSTER_LOGIN_HOST 'mkdir -p /rsrch5/scratch/epi/USER/isoqtl_pilot'
rsync -av --exclude .git ./ USER@CLUSTER_LOGIN_HOST:/rsrch5/scratch/epi/USER/isoqtl_pilot/code/
```

Scratch may be purged on a schedule, so don't leave the only copy of anything there.

**If rsync is blocked,** use an SFTP client such as FileZilla wherever this guide says
`rsync`. Set Transfer > Transfer type to Binary. For the code, upload the folders `R`,
`bin`, `config`, `tests`, `tools` and `docs` into `isoqtl_pilot/code/` and choose
Overwrite. For a bundle, upload the whole folder; its installer checks every file and
names any that arrived incomplete.

## 2. Log in and check the cluster matches

```bash
ssh USER@CLUSTER_LOGIN_HOST
export PILOT=/rsrch5/scratch/epi/$USER/isoqtl_pilot
export CONDA_BIN=/risapps/rhel8/miniforge3/24.5.0-0/bin/conda

ls -l "$CONDA_BIN"                 # conda path from the lab's scripts
bqueues                            # available queues and their time limits
df -h /rsrch5/scratch              # free space
```

If the conda path is missing, run `module avail conda` or `module avail miniforge` and
use what you find.

## 3. Create the environment

Interactive nodes can't reach the web, so the packages are downloaded on a Mac first
and installed offline on the cluster.

**On your Mac,** copy the bundle (about 900 MB). To build a new one, for example after
`environment.yml` changes, run `tools/make_offline_bundle.sh` on a machine with web access.

```bash
rsync -av BUNDLE_DIR/ USER@CLUSTER_LOGIN_HOST:/rsrch5/scratch/epi/USER/isoqtl_pilot/bundle/
```

**On the cluster,** install from it in an interactive job:

```bash
bsub -Is -q short -n 2 -W 2:00 -M 16 -R "rusage[mem=16]" bash

# Inside the interactive job, set these again if they aren't carried over:
export PILOT=/rsrch5/scratch/epi/$USER/isoqtl_pilot
export CONDA_BIN=/risapps/rhel8/miniforge3/24.5.0-0/bin/conda
export CONDA_PKGS_DIRS=$PILOT/conda_pkgs

bash "$PILOT/bundle/install_offline.sh" "$CONDA_BIN" "$PILOT/env"
eval "$("$CONDA_BIN" shell.bash hook)"
conda activate "$PILOT/env"
```

`-p` (used by the installer) keeps the environment beside the project, out of your home
directory quota. Pipeline jobs don't need web access once this is installed.

Check that it installed correctly. Every line should print a version, with no errors:

```bash
R --version | head -1
python --version
salmon --version
STAR --version
plink2 --version | head -1
tabix --version | head -1
Rscript -e 'for (p in c("data.table","R.utils","yaml","corpcor","glmnet","rrBLUP","susieR","Rsubread","edgeR","tximport")) cat(p, as.character(packageVersion(p)), "\n")'
```

You should see R 4.3.x.

## 4. Run the tests

Still in the interactive job with the environment active:

```bash
cd "$PILOT/code"
chmod +x bin/isoqtl-sim
Rscript tests/test_sim_core.R
```

Expect `All checks passed`. This checks the simulation model against known math and
takes about a minute.

Then run the pipeline on synthetic data, which takes a few minutes:

```bash
bash tests/run_toy_e2e.sh "$PILOT/toy_test"
```

Expect `End-to-end test passed`. If it fails, send the end of the output and the log it
names, which will be under `$PILOT/toy_test/results/logs/`.

Type `exit` to leave the interactive job.

## 5. Configure and run the 100-gene pilot

This pilot skips read simulation, so it needs far less compute and disk than a read-level run.

```bash
export PILOT=/rsrch5/scratch/epi/$USER/isoqtl_pilot
export CONDA_BIN=/risapps/rhel8/miniforge3/24.5.0-0/bin/conda
eval "$("$CONDA_BIN" shell.bash hook)"
conda activate "$PILOT/env"
cd "$PILOT/code"

cp config/pilot_lsf.yaml "$PILOT/pilot.yaml"
sed -i "s/CHANGE_ME/$USER/g" "$PILOT/pilot.yaml"
```

Open `$PILOT/pilot.yaml` and check these four things:

1. **`reference`:** the GTF and FASTA exist (`ls -l` each path).
2. **`genotypes.vcf_template`:** this points to another user's directory. Check
   you can read it with the command below. If it has moved or been purged, the VCFs need
   regenerating, so tell me.
3. **`chr_prefix`:** the check below prints chromosome names. Use `""` if you see `22`
   and `"chr"` if you see `chr22`.
4. **`lsf.queue_short` and `lsf.queue_long`:** these exist in `bqueues`, and the long
   queue allows at least 12 hours.

```bash
VCF=/rsrch9/home/epi/bhattacharya_lab/users/sthead/GTEx_gencode_comp/pass1/files_for_analysis/1KG_vcf/genos_1kg_eur_500_snps_maf_0.01_chr22.vcf.gz
ls -l "$VCF" "$VCF.tbi"
tabix -l "$VCF"
```

Select the genes and samples, then generate the job script:

```bash
bsub -Is -q short -n 1 -W 1:00 -M 16 -R "rusage[mem=16]" \
  python bin/isoqtl-sim prepare --config "$PILOT/pilot.yaml"

python bin/isoqtl-sim plan --config "$PILOT/pilot.yaml"
python bin/isoqtl-sim submit --config "$PILOT/pilot.yaml"
```

Read the generated script, then submit it:

```bash
cat "$PILOT/results/submit_lsf.sh"
bash "$PILOT/results/submit_lsf.sh"
```

## 6. Watch the jobs

```bash
bjobs -w                                               # your jobs
python bin/isoqtl-sim status --config "$PILOT/pilot.yaml"   # progress per stage
ls "$PILOT/results/logs/"                              # one log per job
```

**If a job fails:** later stages stay `PEND` because each one waits for the previous
stage to succeed. Read the log, fix the cause, and remove the stuck jobs with
`bkill -J "isoqtl_pilot_*"`. Then rerun the `bash` command from step 5. Finished work
is skipped.

## 7. Send back the results

When everything is done:

```bash
column -t -s $'\t' "$PILOT/results/summary/power_summary.tsv"
python bin/isoqtl-sim status --config "$PILOT/pilot.yaml"
```

Send me that output, plus the end of any failed log if something went wrong. To find
failed logs: `grep -lE "Exited with|TERM_" "$PILOT"/results/logs/*/*.log`.

Scratch may be purged, so copy the summary to your Mac (run this **on your Mac**):

```bash
rsync -av USER@CLUSTER_LOGIN_HOST:/rsrch5/scratch/epi/USER/isoqtl_pilot/results/summary/ ~/isoqtl_pilot_summary/
```

## 8. Update the code and start a new session

Do this whenever the code changes on your Mac, and at the start of each login.

**On your Mac,** copy the code again. Results and the environment are not touched:

```bash
cd ~/Library/CloudStorage/OneDrive-InsideMDAnderson/isoQTL_sim
rsync -av --exclude .git ./ USER@CLUSTER_LOGIN_HOST:/rsrch5/scratch/epi/USER/isoqtl_pilot/code/
```

**On the cluster,** set up the session. These settings are lost when you log out; the
jobs and their results are not:

```bash
export PILOT=/rsrch5/scratch/epi/$USER/isoqtl_pilot
export CONDA_BIN=/risapps/rhel8/miniforge3/24.5.0-0/bin/conda
eval "$("$CONDA_BIN" shell.bash hook)"
conda activate "$PILOT/env"
cd "$PILOT/code"
```

After a code change, repeat the two tests from step 4 in an interactive job. Give the
second one a new directory, because it skips work it finds already finished:

```bash
bsub -Is -q short -n 4 -W 1:00 -M 16 -R "rusage[mem=16]" bash
# then set up the session again as above
Rscript tests/test_sim_core.R
bash tests/run_toy_e2e.sh "$PILOT/toy_test_$(date +%Y%m%d)"
exit
```

## 9. Full run on simulated expression

This simulates every protein-coding autosomal gene (about 19,000) under nine trait
scenarios, without reads. Each run needs its own `output_dir`; the config sets
`results_full_truth`, so the pilot's results are kept.

```bash
cp config/full_truth_lsf.yaml "$PILOT/full_truth.yaml"
sed -i "s/CHANGE_ME/$USER/g" "$PILOT/full_truth.yaml"
grep chr_prefix "$PILOT/pilot.yaml" "$PILOT/full_truth.yaml"    # make these match
```

Check two things before submitting:

1. **Disk.** Run `du -sh "$PILOT/results/simulated"` for the pilot. The full run needs
   about 190 times that, up to double with the extra scenarios. Compare it with
   `df -h /rsrch5/scratch`.
2. **Time.** `grep -h "Run time" "$PILOT"/results/logs/train_truth/*.log` shows what 15
   genes took in the pilot. The full run has 25 genes per job and a 12-hour limit. If
   the pilot's jobs took more than 6 hours, raise the limit with the `resources`
   block in the config.

Then prepare, review and submit as in step 5:

```bash
bsub -Is -q short -n 1 -W 1:00 -M 16 -R "rusage[mem=16]" \
  python bin/isoqtl-sim prepare --config "$PILOT/full_truth.yaml"

python bin/isoqtl-sim plan --config "$PILOT/full_truth.yaml"
python bin/isoqtl-sim submit --config "$PILOT/full_truth.yaml"
cat "$PILOT/results_full_truth/submit_lsf.sh"
bash "$PILOT/results_full_truth/submit_lsf.sh"
```

Expect `plan` to show about 760 jobs each for `simulate` and `train_truth`. Watch
progress as in step 6, with `--config "$PILOT/full_truth.yaml"` and logs under
`$PILOT/results_full_truth/logs/`. To remove stuck jobs, use `bkill -J "isoqtl_full_*"`.

When it finishes, send me:

```bash
column -t -s $'\t' "$PILOT/results_full_truth/summary/power_summary.tsv"
tail -5 "$PILOT/results_full_truth/logs/summarize/1.log"
```

**Reading the summary.** `rejection_rate` is the share of simulated genes with
p < 0.05, and `rejection_rate_bonferroni` uses 0.05 divided by the number of simulated
genes. In the `null` scenario the first should be near 0.05. Genes with too few SNPs
to simulate are counted in `n_genes_not_simulated` and left out of both rates. A
single-isoform gene's isoform test is its gene test, so both levels cover the same
genes. A trained model whose weights are all zero predicts nothing, so it is not
tested; a gene with no other model at that level drops out of `n_genes_modeled` and
counts as a miss. `gene_tests.tsv.gz` has one row per gene, level and scenario, with
`n_isoforms`, for comparing the levels on multi-isoform genes only.

**Recomputing a summary.** After a code update that changes only the summary, there is
no need to rerun the jobs. Keep the old summary and run the summary step directly:

```bash
cp -r "$PILOT/results_full_truth/summary" "$PILOT/results_full_truth/summary_old"
bsub -Is -q short -n 1 -W 1:00 -M 16 -R "rusage[mem=16]" \
  Rscript R/run_step.R summarize "$PILOT/full_truth.yaml" 1
```

## 10. Read-level pilot

This simulates reads for 100 genes in every sample, quantifies them with Salmon, and
repeats the TWAS on Salmon's estimates beside the true expression. It can run at the
same time as step 9.

```bash
cp config/reads_pilot_lsf.yaml "$PILOT/reads_pilot.yaml"
sed -i "s/CHANGE_ME/$USER/g" "$PILOT/reads_pilot.yaml"
grep chr_prefix "$PILOT/pilot.yaml" "$PILOT/reads_pilot.yaml"   # make these match

bsub -Is -q short -n 1 -W 1:00 -M 16 -R "rusage[mem=16]" \
  python bin/isoqtl-sim prepare --config "$PILOT/reads_pilot.yaml"

python bin/isoqtl-sim plan --config "$PILOT/reads_pilot.yaml"
python bin/isoqtl-sim submit --config "$PILOT/reads_pilot.yaml"
cat "$PILOT/results_reads_pilot/submit_lsf.sh"
bash "$PILOT/results_reads_pilot/submit_lsf.sh"
```

Expect these stages, with one job per sample for `reads` and `salmon_gencode_v38`:

```
simulate, truth_matrix, reads, salmon_index_gencode_v38, salmon_gencode_v38,
aggregate_salmon_gencode_v38, train_truth, train_salmon_gencode_v38, summarize
```

The Salmon index is built once from the transcript FASTA and takes the longest of the
single jobs. To remove stuck jobs, use `bkill -J "isoqtl_reads_*"`. When it finishes,
send me:

```bash
column -t -s $'\t' "$PILOT/results_reads_pilot/summary/power_summary.tsv"
python bin/isoqtl-sim status --config "$PILOT/reads_pilot.yaml"
```

The summary has `truth` rows and `salmon_gencode_v38` rows for the same genes. The gap
between them is what quantification costs.

## 11. isoTWAS models

This adds isoTWAS training from the [isotwas](https://github.com/bhattacharya-a-bt/isotwas)
package: the isoforms of a gene are modeled jointly by several methods and the best
is kept per isoform. Results appear as a third level, `isotwas_acat`, beside `gene`
and `isoform_acat` (univariate isoform models).

isotwas needs more R packages, so it gets a new environment, `env_isotwas`. The old
environment and its results are untouched.

**On your Mac,** copy the code and the new bundle (about 1 GB). To rebuild the bundle,
run `tools/make_offline_bundle.sh MICROMAMBA_BIN ~/isoqtl_bundle`.

```bash
cd ~/Library/CloudStorage/OneDrive-InsideMDAnderson/isoQTL_sim
rsync -av --exclude .git ./ USER@CLUSTER_LOGIN_HOST:/rsrch5/scratch/epi/USER/isoqtl_pilot/code/
rsync -av ~/isoqtl_bundle/ USER@CLUSTER_LOGIN_HOST:/rsrch5/scratch/epi/USER/isoqtl_pilot/isoqtl_bundle/
```

**On the cluster,** install and test in an interactive job:

```bash
bsub -Is -q short -n 4 -W 2:00 -M 16 -R "rusage[mem=16]" bash

export PILOT=/rsrch5/scratch/epi/$USER/isoqtl_pilot
export CONDA_BIN=/risapps/rhel8/miniforge3/24.5.0-0/bin/conda
export CONDA_PKGS_DIRS=$PILOT/conda_pkgs
bash "$PILOT/isoqtl_bundle/install_offline.sh" "$CONDA_BIN" "$PILOT/env_isotwas"
```

It first prints `Bundle files verified`, or names the files to copy again. The last
line should be `isotwas 2.0.0 loads`. Seven packages are compiled or installed
from source here; if one fails, the script prints the end of its log. Send me that.

Still in the interactive job, test the new environment:

```bash
eval "$("$CONDA_BIN" shell.bash hook)"
conda activate "$PILOT/env_isotwas"
cd "$PILOT/code"
Rscript tests/test_sim_core.R
bash tests/run_toy_e2e.sh "$PILOT/toy_test_isotwas"
exit
```

Expect `All checks passed` and `End-to-end test passed`.

**Pilot first.** Back on the login node, with `env_isotwas` active (set up the session
as in step 8, but activate `"$PILOT/env_isotwas"`):

```bash
cp config/isotwas_pilot_lsf.yaml "$PILOT/isotwas_pilot.yaml"
sed -i "s/CHANGE_ME/$USER/g" "$PILOT/isotwas_pilot.yaml"
grep chr_prefix "$PILOT/pilot.yaml" "$PILOT/isotwas_pilot.yaml"    # make these match

bsub -Is -q short -n 1 -W 1:00 -M 16 -R "rusage[mem=16]" \
  python bin/isoqtl-sim prepare --config "$PILOT/isotwas_pilot.yaml"
python bin/isoqtl-sim submit --config "$PILOT/isotwas_pilot.yaml"
bash "$PILOT/results_isotwas_pilot/submit_lsf.sh"
```

It uses the same 100 genes as the first pilot. When it finishes, send me:

```bash
column -t -s $'\t' "$PILOT/results_isotwas_pilot/summary/power_summary.tsv"
grep -h "Run time" "$PILOT"/results_isotwas_pilot/logs/train_truth/*.log
grep -h "isotwas" "$PILOT"/results_isotwas_pilot/logs/train_truth/*.log | sort | uniq -c
```

The run times decide the settings for the full run. From a Mac benchmark, isoTWAS
training should take about three times as long as the univariate models, which would
make the pilot's jobs 2 to 8 hours.

**Full run.** `config/full_isotwas_lsf.yaml` runs every gene, 20 per job, with a
24-hour limit for training. Check that the long queue allows 24 hours
(`bqueues -l medium | grep -A2 RUNLIMIT`), then:

```bash
cp config/full_isotwas_lsf.yaml "$PILOT/full_isotwas.yaml"
sed -i "s/CHANGE_ME/$USER/g" "$PILOT/full_isotwas.yaml"
grep chr_prefix "$PILOT/pilot.yaml" "$PILOT/full_isotwas.yaml"     # make these match

bsub -Is -q short -n 1 -W 1:00 -M 16 -R "rusage[mem=16]" \
  python bin/isoqtl-sim prepare --config "$PILOT/full_isotwas.yaml"
python bin/isoqtl-sim plan --config "$PILOT/full_isotwas.yaml"
python bin/isoqtl-sim submit --config "$PILOT/full_isotwas.yaml"
bash "$PILOT/results_full_isotwas/submit_lsf.sh"
```

It repeats the simulation and the univariate models of step 9 with the same seed, so
its `gene` and `isoform_acat` rows should match step 9's closely.

## Troubleshooting

| What you see | Likely cause | What to do |
|---|---|---|
| `conda: command not found` in a job | Wrong conda path in `container_exec` | Set the full path in `pilot.yaml` |
| `Disk quota exceeded` during install | Conda cache in your home directory | Check `CONDA_PKGS_DIRS` is set as in step 3 and retry |
| Connection or `PackagesNotFoundError` errors during install | Conda is trying the web or can't find a bundle file | Confirm `--offline` is in `install_offline.sh` and `bundle/pkgs` copied fully, then send me the error |
| `An error occurred while installing package ... bioconductor-genomeinfodbdata` | A bundle built before `r_data/` existed: that package downloads its content during install | Upload the current bundle's `explicit.txt`, `install_offline.sh` and `r_data` folder, delete the half-built environment folder, and rerun the install |
| Job killed with `TERM_MEMLIMIT` | Not enough memory requested | Raise `mem_gb` for that stage in the config's `resources` block, then run `submit` again |
| Job killed with `TERM_RUNLIMIT` | Not enough time, or the queue limit is too short | Raise `hours` for that stage in `resources`, or change the queue |
| `These settings changed since prepare ran` | The config was edited after `prepare`; finished jobs would keep the old settings | Set a new `output_dir` and run `prepare` again. Queue, `threads` and `resources` changes are always allowed |
| `needs N array jobs, above the limit` | More array jobs than LSF allows in one array | Follow the message: raise `genes_per_job`, then `prepare` again in a new or emptied `output_dir` |
| `Bad resource requirement syntax` or jobs pending on memory | The cluster counts memory in MB, not GB | Set `lsf.mem_per_gb: 1024` |
| `Samples missing from VCF` | Sample file doesn't match the VCF | Check `genotypes.sample_file` |
| All genes report `no_snps` | Wrong `chr_prefix` | Recheck with `tabix -l` |

## Using a container instead (optional)

The `Dockerfile` hasn't been built yet because the Mac used for development has no
Docker. Once it's been built and pushed to a registry, and if the cluster has Apptainer
(check with `module avail apptainer`), you can replace the conda line in `pilot.yaml`:

```yaml
container_exec: apptainer exec --bind /rsrch5,/rsrch9 /path/to/isoqtl-sim.sif
```
