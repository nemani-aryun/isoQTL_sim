#!/bin/bash

#######################################################################
# Step 2: Salmon Quantification
# Quantifies simulated reads using Salmon with specified annotation
# Usage: Called by master script with necessary parameters
#######################################################################

echo "***** Step 2: Salmon Quantification *****"
echo "Execution start time: $(date)"

# Parse arguments
SAMPLE_INDEX=$1
PARAM_ROW_READS=$2
TXOME=$3
TXOME_NAME=$4
THREADS=$5
DIR_PROJECT=$6
DIR_GENOTYPES=$7

echo "Sample index: ${SAMPLE_INDEX}"
echo "Transcriptome: ${TXOME_NAME}"
echo "Threads: ${THREADS}"

# Salmon executable path
salmon=/rsrch5/home/epi/bhattacharya_lab/software/salmon-latest_linux_x86_64/bin/salmon

cd ${DIR_PROJECT}/files_for_analysis

# Get sample name
sample_file="${DIR_GENOTYPES}/files_for_analysis/1kg_eur_500_sample_ids"
sample_name=$(awk -v row=${SAMPLE_INDEX} 'NR == row {print $1}' "$sample_file")

# Determine paired-end status
psr_file="${DIR_PROJECT}/files_for_analysis/parameter_space_reads.txt"
paired_end_status=$(awk -v row=$((PARAM_ROW_READS + 1)) 'NR == row {print $3}' "$psr_file")

# Create output directory
out_dir="salmon/param_row_reads_${PARAM_ROW_READS}_${TXOME_NAME}"
mkdir -p "$out_dir"

# Run Salmon based on read type
if [[ "$paired_end_status" == "T" ]]; then
    echo "Running Salmon for paired-end reads..."
    fqz_file1="${DIR_PROJECT}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz"
    fqz_file2="${DIR_PROJECT}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R2.fastq.gz"
    
    $salmon quant -i ${TXOME} -l A \
      -1 ${fqz_file1} -2 ${fqz_file2} \
      -o ${out_dir}/${sample_name} \
      -p ${THREADS} \
      --validateMappings \
      --seqBias
else
    echo "Running Salmon for single-end reads..."
    fqz_file="${DIR_PROJECT}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz"
    
    $salmon quant -i ${TXOME} -l A \
      --unmatedReads ${fqz_file} \
      -o ${out_dir}/${sample_name} \
      -p ${THREADS} \
      --validateMappings \
      --seqBias
fi

echo "***** Salmon quantification complete for sample ${sample_name} *****"
date