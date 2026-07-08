#!/usr/bin/env Rscript

####################################################################################
# Step 4a: Train Gene-Level TWAS Models
# Using Elastic Net, BLUP, and SuSiE
####################################################################################

library(data.table)
library(dplyr)
library(glmnet)
library(rrBLUP)
library(susieR)
library(VariantAnnotation)

args <- commandArgs(trailingOnly = TRUE)
annot <- as.character(args[1])
param_row_reads <- as.numeric(args[2])
DIR_PROJECT <- as.character(args[3])
DIR_GENOTYPES <- as.character(args[4])

cat("========================================\n")
cat("STEP 4a: Gene-Level TWAS Training\n")
cat("Annotation:", annot, "\n")
cat("========================================\n")

#######################################################################
# Load Data
#######################################################################

# Load BED file (phenotypes)
bed_file <- paste0(DIR_PROJECT, "/files_for_analysis/bed_files/param_row_reads_",
                   param_row_reads, "/salmon_", annot, "_gene_counts.bed")

cat("Loading phenotypes from:", bed_file, "\n")
phenotypes <- fread(bed_file)
phenotypes_matrix <- as.matrix(phenotypes[, 7:ncol(phenotypes)])
rownames(phenotypes_matrix) <- phenotypes$pid
gene_coords <- phenotypes[, 1:6]

# Load genotypes
vcf_dir <- paste0(DIR_GENOTYPES, "/files_for_analysis/1KG_vcf")

#######################################################################
# Helper Functions
#######################################################################

# Extract genotypes for cis window
get_cis_genotypes <- function(gene_row, vcf_dir, cis_window = 1e6) {
  
  chr <- as.numeric(gsub("chr", "", gene_row["#Chr"]))
  gene_start <- as.numeric(gene_row["start"])
  gene_end <- as.numeric(gene_row["end"])
  
  snp_start <- max(gene_start - cis_window, 1)
  snp_end <- gene_end + cis_window
  
  vcf_file <- paste0(vcf_dir, "/genos_1kg_eur_500_snps_maf_0.01_chr", chr, ".vcf.gz")
  
  if (!file.exists(vcf_file)) return(NULL)
  
  query <- paste0(chr, ":", snp_start, "-", snp_end)
  system(paste("tabix -h", vcf_file, query, ">", "/tmp/temp_region.vcf"))
  
  vcf <- readVcf("/tmp/temp_region.vcf", "hg38")
  
  if (nrow(vcf) == 0) return(NULL)
  
  geno_mat <- geno(vcf)$GT
  geno_numeric <- apply(geno_mat, c(1, 2), function(x) {
    as.numeric(substr(x, 1, 1)) + as.numeric(substr(x, 3, 3))
  })
  geno_numeric <- t(geno_numeric)
  
  vcf_gr <- rowRanges(vcf)
  snp_ids <- paste0(seqnames(vcf_gr), "_", start(vcf_gr), "_", 
                    mcols(vcf_gr)$REF, "_", 
                    as.character(unlist(mcols(vcf_gr)$ALT)))
  rownames(geno_numeric) <- snp_ids
  colnames(geno_numeric) <- colnames(phenotypes_matrix)
  
  system("rm /tmp/temp_region.vcf")
  
  return(geno_numeric)
}

# Train TWAS model
train_gene_twas <- function(gene_row, gene_idx, phenotypes_matrix, vcf_dir) {
  
  gene_id <- gene_row["pid"]
  expr <- phenotypes_matrix[gene_id, ]
  
  # Get cis-SNPs
  X <- get_cis_genotypes(gene_row, vcf_dir)
  
  if (is.null(X) || nrow(X) == 0) {
    return(NULL)
  }
  
  # Standardize
  X <- scale(X)
  expr <- scale(expr)
  
  # Skip if no SNPs or low variance
  if (nrow(X) < 10 || var(expr) < 0.01) {
    return(NULL)
  }
  
  best_result <- NULL
  best_r2 <- 0
  
  #######################################################################
  # Method 1: Elastic Net
  #######################################################################
  
  tryCatch({
    enet_model <- cv.glmnet(X, expr, alpha = 0.5, nfolds = 5, 
                            family = "gaussian", standardize = FALSE)
    
    # Calculate R²
    pred <- predict(enet_model, X, s = "lambda.min")
    ss_res <- sum((expr - pred)^2)
    ss_tot <- sum((expr - mean(expr))^2)
    r2_enet <- 1 - (ss_res / ss_tot)
    
    if (r2_enet > best_r2) {
      best_r2 <- r2_enet
      coef_enet <- coef(enet_model, s = "lambda.min")
      best_result <- data.frame(
        SNP = rownames(coef_enet)[-1],
        Weight = as.numeric(coef_enet[-1, 1]),
        Gene = gene_id,
        R2 = r2_enet,
        Method = "ElasticNet",
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    cat("ElasticNet error for gene", gene_idx, "\n")
  })
  
  #######################################################################
  # Method 2: BLUP
  #######################################################################
  
  tryCatch({
    blup_model <- mixed.solve(expr, X = X)
    
    # Calculate R²
    pred <- X %*% blup_model$u
    ss_res <- sum((expr - pred)^2)
    ss_tot <- sum((expr - mean(expr))^2)
    r2_blup <- 1 - (ss_res / ss_tot)
    
    if (r2_blup > best_r2) {
      best_r2 <- r2_blup
      best_result <- data.frame(
        SNP = rownames(X),
        Weight = as.numeric(blup_model$u),
        Gene = gene_id,
        R2 = r2_blup,
        Method = "BLUP",
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    cat("BLUP error for gene", gene_idx, "\n")
  })
  
  #######################################################################
  # Method 3: SuSiE
  #######################################################################
  
  tryCatch({
    susie_model <- susie(X, expr, L = 10, verbose = FALSE)
    
    # Calculate R²
    pred <- X %*% colSums(susie_model$alpha * susie_model$mu)
    ss_res <- sum((expr - pred)^2)
    ss_tot <- sum((expr - mean(expr))^2)
    r2_susie <- 1 - (ss_res / ss_tot)
    
    if (r2_susie > best_r2) {
      best_r2 <- r2_susie
      weights_susie <- colSums(susie_model$alpha * susie_model$mu)
      best_result <- data.frame(
        SNP = rownames(X),
        Weight = as.numeric(weights_susie),
        Gene = gene_id,
        R2 = r2_susie,
        Method = "SuSiE",
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    cat("SuSiE error for gene", gene_idx, "\n")
  })
  
  return(best_result)
}

#######################################################################
# Main Loop
#######################################################################

output_dir <- paste0(DIR_PROJECT, "/results/twas")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

all_models <- data.frame()

for (i in 1:nrow(gene_coords)) {
  
  if (i %% 100 == 0) {
    cat("Training TWAS for gene", i, "of", nrow(gene_coords), "\n")
  }
  
  gene_model <- tryCatch(
    train_gene_twas(gene_coords[i, ], i, phenotypes_matrix, vcf_dir),
    error = function(e) {
      cat("Error for gene", i, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (!is.null(gene_model)) {
    all_models <- rbind(all_models, gene_model)
  }
}

# Save results
output_file <- paste0(output_dir, "/gene_twas_", annot, ".tsv.gz")
fwrite(all_models, file = output_file, sep = "\t", compress = "gzip")

cat("\nGene-level TWAS results saved to:", output_file, "\n")
cat("Total genes with models:", length(unique(all_models$Gene)), "\n")
cat("Mean R²:", format(mean(all_models$R2), digits = 4), "\n")

cat("\n========================================\n")
cat("STEP 4a Complete\n")
cat("========================================\n")