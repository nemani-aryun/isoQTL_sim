#!/bin/bash

#######################################################################
# Step 4b: Train Isoform-Level TWAS Models
# Multivariate methods using isotwas package
#######################################################################

echo "***** STEP 4b: Isoform-Level TWAS Training *****"
echo "Execution start time: $(date)"

# Parse arguments
PARAM_ROW_READS=$1
DIR_PROJECT=$2
DIR_GENOTYPES=$3

echo "Project directory: ${DIR_PROJECT}"

# Load modules
module load R/4.3.1

#######################################################################
# Run isoform TWAS training
#######################################################################

Rscript ${DIR_PROJECT}/scripts/step4b_train_isotwas.R \
    ${PARAM_ROW_READS} \
    ${DIR_PROJECT} \
    ${DIR_GENOTYPES}

# Unload modules
module unload R

echo "***** STEP 4b complete *****"
date