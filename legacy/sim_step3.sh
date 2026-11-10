#!/bin/bash

#######################################################################
# Step 3: eQTL Mapping
# Map cis-SNPs to gene expression using QTLtools
#######################################################################

echo "***** STEP 3: eQTL Mapping *****"
echo "Execution start time: $(date)"

# Parse arguments
ANNOT=$1
PARAM_ROW_READS=$2
DIR_PROJECT=$3
DIR_GENOTYPES=$4

echo "Annotation: ${ANNOT}"
echo "Project directory: ${DIR_PROJECT}"

# Load modules
module load R/4.3.1
module load qtltools
module load tabix
module load samtools

#######################################################################
# Run eQTL mapping
#######################################################################

Rscript ${DIR_PROJECT}/scripts/step3_eqtl_mapping.R \
    ${ANNOT} \
    ${PARAM_ROW_READS} \
    ${DIR_PROJECT} \
    ${DIR_GENOTYPES}

# Unload modules
module unload R
module unload qtltools
module unload tabix
module unload samtools

echo "***** STEP 3 complete *****"
date