# TWAS & isoTWAS Model Training Guide

## Overview

After simulating and quantifying expression (Steps 1-2), the next phase is training **predictive models** that map SNPs to expression. These models enable:
1. **Power analysis** for eQTL and TWAS studies
2. **Comparison** of gene-level vs isoform-level prediction accuracy
3. **Understanding** which isoforms drive associations

---

## Why Train Models?

### Gene-Level TWAS (Traditional)
- Fast, robust, standard approach
- **Limitation**: Loses isoform-level information
- Cannot distinguish which transcript(s) are causal

### Isoform-Level TWAS (Novel)
- Captures which specific isoforms respond to genetic variation
- Can model isoform co-expression and collinearity
- Better functional interpretation

---

## Implementation Plan

### Step 4a: Gene-Level TWAS Model Training

**Input**: 
- `bed_files/salmon_gencode_vXX_gene_counts.bed` (gene expression × 500 samples)
- Genotype dosages (SNPs × 500 samples)

**Process** (pseudocode):
```r
for each gene in {simulated genes}:
  1. Extract SNPs within ±1Mb of gene
  2. Standardize SNP dosages and expression
  3. Train three methods:
     a) Elastic Net (cv.glmnet, alpha=0.5)
     b) BLUP (rrBLUP package)
     c) SuSiE (fine-mapping with sparsity)
  4. Perform 5-fold cross-validation
  5. Select best model by R² (explained variance)
  6. Save: SNP weights, R², effect sizes
```

**Output**:
- `TWAS_models/geneID_TWAS.tsv.gz`

**Script**: `step4a_train_gene_twas.R` (to be created)

---

### Step 4b: Isoform-Level TWAS Model Training

**Input**:
- Isoform-level expression (transcript quantifications, TPM)
- Genotype dosages

**Novel Components**:

1. **Isoform Panel Definition**
   - Filter to isoforms expressed in ≥75% samples (>0.1 TPM)
   - For genes with >20 isoforms, keep only top 95% by cumulative expression
   - Learn isoform correlation structure (Omega matrix)

2. **Multivariate Training Methods**
```r
   Methods:
   - multi_enet: Multivariate elastic net (shared penalties across isoforms)
   - mvsusie: Multivariate SuSiE (joint fine-mapping)
   - univariate: Train each isoform separately (baseline)
   - mrce_lasso: Multivariate with robust covariance estimation
```

3. **Covariance Estimation**
   - Estimate isoform co-expression (Omega) from the data
   - Use replicated measures or bootstrap if needed

4. **Cross-Validation**
   - Nested CV: inner CV for parameter tuning, outer for R² estimation
   - Stratified by expression quartiles to avoid bias

5. **Output**:
   - Per isoform: SNP weights × isoforms matrix
   - Isoform correlation structure (Omega)
   - Per-isoform R² values
   - Multivariate R² (explained variance across isoforms)

**Script**: `step4b_train_isoform_twas.R` (to be created, based on isotwas package)

---

## Example Code (Adapted from isotwas)

```r
# This is pseudocode/template - modify for your needs

library(isotwas)
library(glmnet)
library(susieR)
library(rrBLUP)

# === STEP 4a: GENE-LEVEL TWAS ===

# Load expression and genotypes
gene_expr <- fread("bed_files/salmon_gencode_v38_gene_counts.bed")
genotypes <- load_genotypes(sample_ids, snp_ranges)

train_gene_twas <- function(gene_id, genotypes, expression) {
  
  # Extract cis-SNPs (±1Mb from gene)
  snps <- get_cis_snps(gene_id, genotypes, cis_window=1e6)
  
  # Prepare matrices
  X <- as.matrix(genotypes[, snps])
  y <- expression[, gene_id]
  
  # Scale
  X <- scale(X)
  y <- scale(y)
  
  # Three competing methods
  enet_model <- cv.glmnet(X, y, alpha=0.5, nfolds=5, family="gaussian")
  blup_model <- mixed.solve(y, X=X)  # rrBLUP
  susie_model <- susie(X, y, L=10)   # SuSiE
  
  # Select best by R²
  enet_r2 <- calculate_r2(X, y, enet_model)
  blup_r2 <- calculate_r2(X, y, blup_model)
  susie_r2 <- calculate_r2(X, y, susie_model)
  
  best_r2 <- max(c(enet_r2, blup_r2, susie_r2))
  
  # Save weights
  weights_out <- data.frame(
    SNP = snps,
    Weight = get_weights_best_model(...),
    R2 = best_r2,
    Method = c("enet", "blup", "susie")[which.max(c(enet_r2, blup_r2, susie_r2))]
  )
  
  return(weights_out)
}

# === STEP 4b: ISOFORM-LEVEL TWAS ===

train_isoform_twas <- function(gene_id, genotypes, isoform_expr) {
  
  # Get isoforms for this gene
  tx_ids <- colnames(isoform_expr)[grepl(gene_id, colnames(isoform_expr))]
  Y <- as.matrix(isoform_expr[, tx_ids])
  
  # Filter isoforms: keep those expressed >0.1 TPM in ≥75% samples
  expressed <- colMeans(Y > 0.1) >= 0.75
  Y <- Y[, expressed]
  
  # Get cis-SNPs
  snps <- get_cis_snps(gene_id, genotypes)
  X <- as.matrix(genotypes[, snps])
  X <- scale(X)
  Y <- scale(Y)
  
  # Estimate isoform covariance (Omega)
  Omega <- cov(Y)
  
  # Train multivariate models
  m_enet <- compute_isotwas(
    X = X, 
    Y = Y,
    method = c('multi_enet', 'mvsusie', 'univariate'),
    omega = Omega,
    nfolds = 5
  )
  
  # Return: weights per isoform + Omega
  return(list(
    weights = m_enet$weights,
    omega = Omega,
    r2 = m_enet$r2,
    univariate_r2 = m_enet$univariate_r2
  ))
}
```

---

## Next Steps

1. **Adapt isotwas code** to your data structure
   - Modify for single dataset vs cross-validation setup
   - Test on 1000 genes first (power calculation design)
   
2. **Implement gene-level TWAS** (Step 4a)
   - Simpler baseline, no dependencies on isoform methods
   
3. **Implement isoform-level TWAS** (Step 4b)
   - Build on Step 4a infrastructure
   - Use isotwas R package functions

4. **Validation**
   - Train on 90% samples, test on 10%
   - Compare predictions vs held-out data
   - Check cross-validation curves for overfitting

5. **Comparison**
   - isoform-level R² vs gene-level R²
   - Which isoforms drive associations?
   - Power to detect isoform-specific effects
