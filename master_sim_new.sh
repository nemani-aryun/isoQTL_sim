#!/bin/bash

####################################################################################
# Master Script for Step 1: Expression and Read Simulation
# Usage: Customize directories and parameters below, then submit to cluster
####################################################################################

#######################################################################
# USER-DEFINED PARAMETERS
#######################################################################

# Pass identifiers
PASS=1
GENO_PASS=1

# Expression simulation settings
GENES_PER_JOB=15
TOTAL_GENES=12959
TOTAL_JOBS=$(( (TOTAL_GENES + GENES_PER_JOB - 1) / GENES_PER_JOB ))

# Read parameter space row (which RNA-seq settings to use)
PARAM_ROW_READS=1

# Directory paths (CUSTOMIZE THESE)
DIR_PROJECT="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}"
DIR_GENOTYPES="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${GENO_PASS}"
DIR_SCRIPTS="${DIR_PROJECT}/scripts"

#######################################################################
# SETUP
#######################################################################

cd $DIR_PROJECT

# Create log directories
mkdir -p logs/step1_sim

echo "=========================================="
echo "Master Script: Step 1 Simulation"
echo "=========================================="
echo "Project directory: ${DIR_PROJECT}"
echo "Genotype directory: ${DIR_GENOTYPES}"
echo "Scripts directory: ${DIR_SCRIPTS}"
echo "Total genes: ${TOTAL_GENES}"
echo "Genes per job: ${GENES_PER_JOB}"
echo "Total jobs: ${TOTAL_JOBS}"
echo "=========================================="

#######################################################################
# SUBMIT STEP 1: EXPRESSION & READ SIMULATION
#######################################################################

# This job array handles both expression simulation AND read generation
# - Job indices 1-864: Simulate expression (batches of 15 genes each)
# - Job indices 1-500: Also simulate reads (one sample per index)
# - Job indices 501-864: Only simulate expression (no read simulation)

job_id_step1=$(bsub -n 1 -W 03:00 -J "step1_sim[1-$TOTAL_JOBS]" \
    -q short -M 30 -R rusage[mem=30] \
    -o logs/step1_sim/%I.log \
    ${DIR_SCRIPTS}/step1_sim_expr_and_reads.sh \
        ${PASS} \
        ${GENO_PASS} \
        ${GENES_PER_JOB} \
        ${TOTAL_GENES} \
        ${PARAM_ROW_READS} \
        ${DIR_PROJECT} \
        ${DIR_GENOTYPES} \
        ${DIR_SCRIPTS} | awk '{print $2}' | tr -d '<>')

echo "Submitted Step 1 job array with ID: $job_id_step1"
echo "Job array indices: 1-${TOTAL_JOBS}"
echo "Expression simulation: All ${TOTAL_JOBS} jobs"
echo "Read simulation: Jobs 1-500 only"
echo "=========================================="
echo "Monitor with: bjobs -A"
echo "Check logs in: ${DIR_PROJECT}/logs/step1_sim/"
echo "=========================================="

# Return job ID for potential downstream dependencies
echo $job_id_step1