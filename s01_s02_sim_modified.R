#!/usr/bin/env Rscript

####################################################################################
# Combined Expression and Read Simulation Script
# Part 1: Simulates ground-truth isoform expression for assigned genes
# Part 2: Converts simulated expression to FASTQ files for one sample
####################################################################################

####################################################################################
# change library to local
####################################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3", myPaths)
.libPaths(myPaths)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
gene_start <- as.numeric(args[1])        # Starting gene index for this job
gene_end <- as.numeric(args[2])          # Ending gene index for this job
sample <- as.numeric(args[3])            # Sample index (0 if > 500)
param_row_reads <- as.numeric(args[4])   # Read parameter row
DIR_PROJECT <- as.character(args[5])     # Project directory
DIR_GENOTYPES <- as.character(args[6])   # Genotype directory

####################################################################################
# load dependencies silently
####################################################################################

load_silent <- function(packages) {
  invisible(lapply(packages, function(pkg) {
    suppressPackageStartupMessages(suppressWarnings(library(pkg, character.only = TRUE)))
  }))
}

load_silent(c("MASS", "corpcor", "mvnfast", "data.table", "dplyr",
              "VariantAnnotation", "rtracklayer", "Biostrings", 
              "Matrix", "Rsubread", "stringr"))

####################################################################################
# PART 1: EXPRESSION SIMULATION
####################################################################################

cat("=========================================\n")
cat("PART 1: EXPRESSION SIMULATION\n")
cat("=========================================\n")

####################################################################################
# helper functions
####################################################################################

vcf2Geno <- function(vcf){
    apply(vcf,c(1,2),FUN=function(x){
        as.numeric(substr(x,1,1))+as.numeric(substr(x,3,3))
    })
}

simulate_count <- function(gene_num, anno, anno_gene, param_space, vcf_dir, tmp_dir, expr_dir) {

  set.seed(gene_num)
  
  query <- anno_gene$query[gene_num]
  start <- anno_gene$start[gene_num]
  end <- anno_gene$end[gene_num]
  chr <- anno_gene$chr[gene_num]
  gene <- anno_gene$GeneID[gene_num]
  
  cat("Simulating gene", gene_num, "-", gene, "\n")

  # use query to subset the VCF with tabix
  subset_vcf_path <- paste0(tmp_dir, "/subset_", chr, "_", start, "_", end, ".vcf")

  if (!file.exists(paste0(subset_vcf_path,".gz"))) {
    system(paste("tabix -h", paste0(vcf_dir, "/genos_1kg_eur_500_snps_maf_0.01_chr", chr, ".vcf.gz"), 
                 query, ">", subset_vcf_path))
    system(paste("bgzip", subset_vcf_path))
    system(paste0("tabix -p vcf ",subset_vcf_path,".gz"))
  }

  vcf_file <- paste0(tmp_dir,"/subset_", chr,"_",start,"_",end, ".vcf.gz")
  vcf <- readVcf(vcf_file,"hg38")
  genotypes <- geno(vcf)$GT

  vcf_gr <- rowRanges(vcf)
  alt_alleles <- as.character(unlist(mcols(vcf_gr)$ALT))
  snp_ids <- paste0(seqnames(vcf_gr), "_", start(vcf_gr), "_", 
                    mcols(vcf_gr)$REF, "_", alt_alleles)

  n_snp <- nrow(genotypes)

  # Get parameters for this gene
  gene_params <- param_space[param_space$GeneID==gene,]
  n_causal <- gene_params$n_causal
  prop_shared <- gene_params$prop_shared
  h2_g <- gene_params$h2_g
  cov_E_min <- gene_params$cov_E_min
  cov_E_max <- gene_params$cov_E_max
  cov_B_min <- gene_params$cov_B_min
  cov_B_max <- gene_params$cov_B_max

  M <- sum(anno$GeneID == gene)
  K <- n_causal
  S <- ceiling(prop_shared * K)
  D <- K - S
  W <- S + D * M
  B <- matrix(0, nrow=W, ncol=M)

  causal_ind <- sample(1:n_snp, W)
  X <- genotypes[causal_ind, ]
  X <- vcf2Geno(X)
  X <- t(X)
  X <- scale(X, center=TRUE, scale=TRUE)
  N <- nrow(X)

  cov_B <- diag(1, M)
  cov_B[lower.tri(cov_B)] <- runif(M * (M-1) / 2, cov_B_min, cov_B_max)
  cov_B[upper.tri(cov_B)] <- t(cov_B)[upper.tri(cov_B)]
  cov_B <- cov2cor(cov_B)
  if (!is.positive.definite(cov_B)) cov_B <- make.positive.definite(cov_B)

  if (S > 0) {
      h2_per_snp <- if(length(h2_g)==1) h2_g / K else (h2_g / K)
      shared_beta_mat <- mvrnorm(n = S, mu = rep(0, M), Sigma = cov_B)
      shared_beta_mat <- sweep(shared_beta_mat, 2, sqrt(h2_per_snp), `*`)
      B[1:S, ] <- shared_beta_mat
  }

  if (D > 0) {
      unique_beta_vec <- rnorm(D * M, 0, sqrt(h2_g / K))
      for (m in 1:M) {
          idx_start <- S + (m - 1) * D + 1
          idx_end   <- S + m * D
          B[idx_start:idx_end, m] <- unique_beta_vec[((m - 1) * D + 1):(m * D)]
      }
  }

  genetic_component <- X %*% B
  genetic_var <- diag(cov(genetic_component))
  scaling_factors <- sqrt(h2_g / genetic_var)
  B <- B * scaling_factors
  genetic_var <- diag(cov(X %*% B))

  causal_ids <- snp_ids[causal_ind]
  rownames(B) <- causal_ids

  target_variance <- 1 - h2_g
  cov_E <- diag(1, M)
  cov_E[upper.tri(cov_E)] <- runif(M * (M-1) / 2, cov_E_min, cov_E_max)
  cov_E[lower.tri(cov_E)] <- t(cov_E)[lower.tri(cov_E)]
  cov_E <- cov_E * target_variance / mean(diag(cov_E))
  if (!is.positive.definite(cov_E)) cov_E <- make.positive.definite(cov_E)

  E <- mvrnorm(n=N, mu=rep(0, M), Sigma=cov_E)
  Y <- X%*%B + E

  total_var <- diag(cov(Y))
  heritability <- genetic_var / total_var
  transcript_ids <- anno$Transcript[anno$GeneID == gene]

  Sigma_G <- cov(X %*% B)
  Sigma_Y <- cov(Y)
  one_vec <- rep(1, M)
  gene_h2 <- as.numeric( (t(one_vec) %*% Sigma_G %*% one_vec) /
                         (t(one_vec) %*% Sigma_Y %*% one_vec) )
  
  result <- list(
    Y = Y,
    transcript_ids = transcript_ids,
    B = B,
    isoform_h2 = heritability,
    gene_h2 = gene_h2
  )
  
  save(result, file=paste0(expr_dir, "/expr_", gene, ".RData"))
  return(result)
}

####################################################################################
# Run expression simulation for assigned genes
####################################################################################

window <- 250000
param_space <- fread(paste0(DIR_PROJECT, "/files_for_analysis/gene_level_parameters.txt"))
vcf_dir <- paste0(DIR_GENOTYPES, "/files_for_analysis/1KG_vcf")
anno <- read.table(paste0(DIR_PROJECT, "/files_for_analysis/anno_selected_genes.txt"), header=T)

anno_gene <- anno %>%
  group_by(GeneID) %>%
  summarise(start = unique(start), end = unique(end), 
            seqnames = unique(seqnames), strand = unique(strand))

anno_gene$window_start <- anno_gene$start - window
anno_gene$window_start[anno_gene$window_start<1] <- 1
anno_gene$window_end <- anno_gene$end + window
anno_gene$chr <- as.numeric(substr(anno_gene$seqnames,4,nchar(anno_gene$seqnames)))
anno_gene$query <- paste0(anno_gene$chr,":",anno_gene$window_start,"-",anno_gene$window_end)

tmp_dir <- paste0(DIR_GENOTYPES, "/files_for_analysis/tmp")
expr_dir <- paste0(DIR_PROJECT, "/files_for_analysis/expression")

if(!dir.exists(expr_dir)) dir.create(expr_dir, recursive=T)
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive=T)

# Simulate expression for assigned gene range
for (gene_num in gene_start:gene_end) {
  Y <- simulate_count(gene_num, anno, anno_gene, param_space, vcf_dir, tmp_dir, expr_dir)
}

cat("Expression simulation complete for genes", gene_start, "to", gene_end, "\n")

####################################################################################
# PART 2: READ SIMULATION
####################################################################################

# Only run if sample index is valid (1-500)
if (sample > 0 && sample <= 500) {
  
  cat("\n=========================================\n")
  cat("PART 2: READ SIMULATION\n")
  cat("Sample:", sample, "\n")
  cat("=========================================\n")
  
  param_space_reads <- fread(paste0(DIR_PROJECT, "/files_for_analysis/parameter_space_reads.txt"))
  library_size <- param_space_reads$library_size[param_row_reads]
  read_length <- param_space_reads$read_length[param_row_reads]
  paired_end <- param_space_reads$paired_end[param_row_reads]
  paired_end <- as.logical(paired_end)

  reads_dir <- paste0(DIR_PROJECT, "/files_for_analysis/reads")
  if(!dir.exists(reads_dir)) dir.create(reads_dir)

  fa.file <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.transcripts.fa")
  transcripts <- scanFasta(fa.file)
  tx_id <- str_extract(transcripts$TranscriptID, "^[^|]+")

  if(file.exists(paste0(DIR_PROJECT, "/files_for_analysis/expression/aggregatedY.RData"))){
     load(paste0(DIR_PROJECT, "/files_for_analysis/expression/aggregatedY.RData"))
  } else {
    Y_files <- list.files(paste0(DIR_PROJECT, "/files_for_analysis/expression"), 
                          pattern = "expr_.*\\.RData$", full.names = T)
    
    Y_list <- lapply(Y_files, function(file) {
      load(file)
      return(result)  # Changed from Y to result
    })

    Y_full <- do.call(cbind, lapply(Y_list, `[[`, 1))
    tx_full <- unlist(lapply(Y_list, `[[`, 2))
    save(Y_full, tx_full, file=paste0(DIR_PROJECT, "/files_for_analysis/expression/aggregatedY.RData"))
  }

  Y <- data.frame(expr=Y_full[sample,],
                  tx=str_extract(tx_full, "^[^|]+"))

  Y$expr <- exp(Y$expr)
  id <- rownames(Y_full)[sample]

  TPMs <- data.frame(tx=tx_id)
  TPMs$id <- 1:nrow(TPMs)

  dat <- merge(x=TPMs, y=Y, by="tx", all.x=T, all.y=F)
  dat$expr[is.na(dat$expr)] <- 0
  dat <- dat[order(dat$id),]

  setwd(reads_dir)

  simReads(transcript.file = fa.file,
           expression.levels=dat$expr,
           output.prefix=paste0('sim_',id,"_param_row_reads_",param_row_reads),
           read.length = read_length,         
           library.size = library_size,     
           paired.end = paired_end)
  
  cat("Read simulation complete for sample", sample, "\n")
  
} else {
  cat("\nSkipping read simulation (sample index =", sample, ")\n")
}

cat("\n=========================================\n")
cat("SCRIPT COMPLETE\n")
cat("=========================================\n")