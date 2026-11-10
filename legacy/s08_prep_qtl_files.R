#!/usr/bin/env Rscript

###################################################################
# change library to local
###################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3",myPaths)
.libPaths(myPaths)

####################################################################################
# load dependencies
####################################################################################
library(stringr)
library(VariantAnnotation)
library(rtracklayer)
library(Biostrings)
library(dplyr)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
pass <- as.numeric(args[1])
param_row_reads <- as.numeric(args[2])
annot <- as.character(args[3])

####################################################################################
# begin code
####################################################################################

# bed file should have the following columns:
# Chromosome ID [string]
# Start genomic position of the phenotype (here the TSS of gene1) [integer, 0-based]
# End genomic position of the phenotype (here the TSS of gene1) [integer, 1-based]
# Phenotype ID (here the exon IDs) [string].
# Phenotype group ID (here the gene IDs, multiple exons belong to the same gene) [string]
# Strand orientation [+/-]

se <- readRDS(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/quants/param_row_reads_",param_row_reads,"_salmon_",annot,"_gene.RDS"))

# Get gene-level TPMs
counts_gene

tpm_gene <- assay(se, "abundance")

# Compute proportion of samples > 0.1 TPM
prop_expr_gene <- rowMeans(tpm_gene > 0.1)

# Filter genes
se.filt <- se[prop_expr_gene >= 0.75, ]


gr <- rowRanges(se)
bed <- data.frame(
  Chr    = as.character(seqnames(gr)),
  start  = start(gr) - 1,              # BED is 0-based
  end    = end(gr),                    # inclusive end
  pid    = names(gr),                  # "pid" – rownames like ENSG...
  gid    = mcols(gr)$gene_id,          # gene IDs
  strand = as.character(strand(gr))
)


quant <- data.frame(assay(se))  # will grab the default assay
# quant <- quant %>%
# 	mutate(tx_id=str_extract(rownames(quant), "^[^|]+"))
quant$GeneID <- rownames(quant)

out_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/bed_files/param_row_reads_",param_row_reads)

if(!dir.exists(out_dir)){
	dir.create(out_dir,recursive=T)
}

bed_merged <- merge(x=bed,y=quant,by.x="gid",by.y="GeneID")
# bed_merged <- bed_merged[,c(2,3,4,1,5,6,7:(ncol(bed_merged))),]
# bed_merged$seqnames <- chr
bed_merged <- bed_merged[,c(2,3,4,1,1,6,7:ncol(bed_merged))]
bed_merged <- bed_merged[bed_merged$Chr %in% paste0("chr",1:22),]
bed_merged$Chr <- factor(bed_merged$Chr, levels=paste0("chr",1:22))
bed_merged <- bed_merged[order(bed_merged$start),]
bed_merged <- bed_merged[order(bed_merged$Chr),]

colnames(bed_merged)[1:6] <- c("#Chr","start","end","pid","gid","strand")

# bed_merged_scaled <- bed_merged
# # extract the scale factor for all rows in one go 
# scale_factors <- sapply(bed_merged$pid, function(tx) tx_scale_factor[grep(tx, tx_scale_names)])

# # divide only the columns from 7 onwards by the corresponding scale factor
# bed_merged_scaled[,-1:-6] <- bed_merged[,-1:-6] / scale_factors

write.table(bed_merged,file=paste0(out_dir,"/salmon_",annot,"_gene_counts.bed"),sep="\t",quote=F,col.names=T,row.names=F)
