#!/bin/bash
# End-to-end pipeline test on synthetic data. Usage: run_toy_e2e.sh WORK_DIR
# Passes only if every stage finishes, including rerun skipping and required outputs.
set -euo pipefail
work=${1:?usage: run_toy_e2e.sh WORK_DIR}
root=$(cd "$(dirname "$0")/.." && pwd)
cli="python3 $root/bin/isoqtl-sim"
mkdir -p "$work"
python3 "$root/tests/make_toy_data.py" "$work/data"

cat > "$work/config.yaml" <<EOF
output_dir: $work/results
seed: 1
genes_per_job: 3
threads: 2
run_reads: true
truth_annotation: toy
reference:
  toy: {gtf: $work/data/toy.gtf, fasta: $work/data/toy.fa}
genotypes:
  vcf_template: $work/data/toy_chr{chr}.vcf.gz
  chr_prefix: ""
windows: {causal: 100000, train: 150000}
simulation:
  n_genes: all
  n_causal: 4
  prop_shared: [0.5, 0.75, 1]
  h2_g: 0.4
  cov_B: [0.7, 0.9]
  cov_E: [0, 0.5]
trait:
  scenarios:
    - {name: null_100k, h2ge: 0, n_gwas: 100000, mediator: gene}
    - {name: alt_100k, h2ge: 0.005, n_gwas: 100000, mediator: gene}
    - {name: iso_100k, h2ge: 0.005, n_gwas: 100000, mediator: isoform}
reads: {read_length: 75, library_size: 200000, paired_end: true}
quant: {salmon: [toy], min_prop_expressed: 0.5}
models:
  methods: [enet, blup]
  isotwas: {methods: [multi_enet, univariate]}
summary: {alpha: 0.05}
EOF

$cli prepare --config "$work/config.yaml"
$cli plan --config "$work/config.yaml"
$cli run --config "$work/config.yaml" -j 4
$cli status --config "$work/config.yaml"

summary="$work/results/summary/power_summary.tsv"
for required in "$summary" "$work/results/bed/salmon_toy_gene.bed" \
    "$work/results/simulated_all/isoform_log_expression.rds"; do
  test -s "$required" || { echo "FAIL: missing $required"; exit 1; }
done
grep -q "^truth" "$summary" || { echo "FAIL: no truth results"; exit 1; }
grep -q "^salmon_toy" "$summary" || { echo "FAIL: no read-level results"; exit 1; }
for source in truth salmon_toy; do
  grep -q "^$source	isotwas_acat" "$summary" || { echo "FAIL: no isoTWAS results for $source"; exit 1; }
done
# Gene and isoform tests must cover the same genes, single-isoform genes included
# (column 7 is n_genes_simulated, column 9 is n_genes_modeled).
awk -F'\t' '$1 == "truth" && $9 != $7 { bad = 1 } END { exit bad }' "$summary" \
  || { echo "FAIL: a truth-level test does not cover every simulated gene"; exit 1; }

# A changed setting must be refused, because finished jobs would keep the old one.
sed 's/h2ge: 0.005, n_gwas: 100000, mediator: gene/h2ge: 0.002, n_gwas: 100000, mediator: gene/' \
  "$work/config.yaml" > "$work/changed.yaml"
if $cli plan --config "$work/changed.yaml" > /dev/null 2>&1; then
  echo "FAIL: a changed config was accepted for an existing output_dir"; exit 1
fi

before=$(find "$work/results/.done" -type f | wc -l)
start=$SECONDS
$cli run --config "$work/config.yaml" -j 4 > /dev/null
test "$(find "$work/results/.done" -type f | wc -l)" -eq "$before"
echo "Rerun skipped finished work in $((SECONDS - start))s"

$cli submit --config "$work/config.yaml"
bash -n "$work/results/submit_lsf.sh"
echo "End-to-end test passed"
