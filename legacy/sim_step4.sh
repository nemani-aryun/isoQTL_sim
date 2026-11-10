#!/bin/bash

#######################################################################
# Step 4a: Train Gene-Level TWAS Models
# Univariate methods: Elastic Net, BLUP, SuSiE
#######################################################################

echo "***** STEP 4a: Gene-Level TWAS Training *****"
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

#######################################################################
# Run TWAS training
#######################################################################

Rscript ${DIR_PROJECT}/scripts/step4a_train_gene_twas.R \
    ${ANNOT} \
    ${PARAM_ROW_READS} \
    ${DIR_PROJECT} \
    ${DIR_GENOTYPES}

# Unload modules
module unload R

echo "***** STEP 4a complete *****"
date