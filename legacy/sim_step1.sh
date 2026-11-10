#!/bin/bash

#######################################################################
# Step 1: Simulate Expression and Generate Reads
# This script runs both expression simulation and read generation
# Usage: Called by master script with all necessary parameters
#######################################################################

echo "***** Step 1: Expression & Read Simulation *****"
echo "Job ID: $LSB_JOBID"
echo "Job index: $LSB_JOBINDEX"
echo "Node: $(hostname)"
echo "Execution start time: $(date)"

# Parse arguments from master script
PASS=$1
GENO_PASS=$2
GENES_PER_JOB=$3
TOTAL_GENES=$4
PARAM_ROW_READS=$5
DIR_PROJECT=$6
DIR_GENOTYPES=$7
DIR_SCRIPTS=$8

# Load required modules
module load R/4.3.1
module load tabix/0.2.6

#######################################################################
# Part A: Expression Simulation (runs once per job array task)
#######################################################################

echo ">>> Starting expression simulation for job array index $LSB_JOBINDEX"

# Calculate gene range for this job
JOB_INDEX=${LSB_JOBINDEX}
START_INDEX=$(( (JOB_INDEX - 1) * GENES_PER_JOB + 1 ))
END_INDEX=$(( JOB_INDEX * GENES_PER_JOB ))

# Handle boundary case
if [ $END_INDEX -gt $TOTAL_GENES ]; then
    END_INDEX=$TOTAL_GENES
fi

# Loop through assigned genes
for ((i=START_INDEX; i<=END_INDEX; i++))
do
    echo "Simulating expression for gene $i"
    Rscript ${DIR_SCRIPTS}/s01_sim_expr.R \
        $i \
        ${PASS} \
        ${GENO_PASS} \
        ${DIR_PROJECT} \
        ${DIR_GENOTYPES}
done

echo ">>> Expression simulation complete"

#######################################################################
# Part B: Read Simulation (runs for sample corresponding to job index)
#######################################################################

# Only run read simulation if this is a valid sample index (1-500)
if [ $LSB_JOBINDEX -le 500 ]; then
    echo ">>> Starting read simulation for sample $LSB_JOBINDEX"
    
    Rscript ${DIR_SCRIPTS}/s02_sim_reads.R \
        ${LSB_JOBINDEX} \
        ${PASS} \
        ${PARAM_ROW_READS} \
        ${DIR_PROJECT}
    
    echo ">>> Read simulation complete for sample $LSB_JOBINDEX"
else
    echo ">>> Skipping read simulation (job index > 500)"
fi

# Unload modules
module unload R
module unload tabix

echo "***** Step 1 complete *****"
date