#!/usr/bin/env Rscript
# Pipeline steps, called by bin/isoqtl-sim:
#   Rscript R/run_step.R <step> <config.yaml> [index]

args_all <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args_all[grep("--file=", args_all)]))
source(file.path(script_dir, "sim_core.R"))

args <- commandArgs(trailingOnly = TRUE)
step <- args[1]
cfg <- yaml::read_yaml(args[2])
index <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_

suppressPackageStartupMessages(library(data.table))

out <- function(...) file.path(cfg$output_dir, ...)
dir_for <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
# Write to a temporary name first so a killed job never leaves a partial output.
save_rds <- function(x, path) {
  dir_for(path)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(x, tmp)
  file.rename(tmp, path)
}
strip_version <- function(id) sub("\\..*$", "", id)
genes_table <- function() fread(out("params", "genes.tsv"))
transcripts_table <- function() fread(out("params", "transcripts.tsv"))
samples <- function() readLines(out("params", "samples.txt"))

chunk_rows <- function(genes) {
  n <- cfg$genes_per_job
  seq(from = (index - 1) * n + 1, to = min(index * n, nrow(genes)))
}

gtf_attr <- function(x, key) {
  m <- regmatches(x, regexpr(paste0(key, ' "[^"]+"'), x))
  out <- rep(NA_character_, length(x))
  out[grepl(paste0(key, ' "'), x)] <- sub(paste0(key, ' "([^"]+)"'), "\\1", m)
  out
}

read_gtf_genes <- function(path) {
  gtf <- fread(cmd = sprintf("gzip -cdf '%s' | grep -v '^#' | awk -F'\\t' '$3==\"gene\"'", path),
               sep = "\t", header = FALSE, quote = "")
  data.table(gene_id = gtf_attr(gtf$V9, "gene_id"), gene_type = gtf_attr(gtf$V9, "gene_type"),
             seqnames = gtf$V1, start = gtf$V4, end = gtf$V5, strand = gtf$V7)
}

# Map transcript to gene using GENCODE FASTA headers:
# >ENST...|ENSG...|...  Falls back to a two-column tx2gene for other annotations.
read_fasta_transcripts <- function(path) {
  headers <- fread(cmd = sprintf("gzip -cdf '%s' | grep '^>' | cut -c2-", path),
                   sep = "\n", header = FALSE, quote = "")$V1
  parts <- tstrsplit(headers, "|", fixed = TRUE, keep = 1:2)
  data.table(header = headers, tx_id = parts[[1]], gene_id = parts[[2]])
}

vcf_path <- function(chr) sub("{chr}", chr, cfg$genotypes$vcf_template, fixed = TRUE)
region_name <- function(chr, start, end) {
  sprintf("%s%s:%d-%d", cfg$genotypes$chr_prefix %||% "", chr, as.integer(start), as.integer(end))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Read a genotype window. Returns samples x SNPs dosages with biallelic,
# polymorphic SNPs only, in the order of params/samples.txt.
read_window <- function(chr, start, end, ids) {
  vcf <- vcf_path(chr)
  header <- system2("tabix", c("-H", shQuote(vcf)), stdout = TRUE)
  cols <- strsplit(header[length(header)], "\t")[[1]]
  lines <- system2("tabix", c(shQuote(vcf), region_name(chr, start, end)), stdout = TRUE)
  if (!length(lines)) return(NULL)
  dt <- fread(text = lines, sep = "\t", header = FALSE, quote = "", colClasses = "character")
  setnames(dt, cols)
  dt <- dt[!grepl(",", ALT, fixed = TRUE) & nchar(REF) == 1 & nchar(ALT) == 1]
  if (!nrow(dt)) return(NULL)
  missing_ids <- setdiff(ids, cols)
  if (length(missing_ids)) stop("Samples missing from VCF: ", paste(head(missing_ids), collapse = ", "))
  gt <- as.matrix(dt[, ..ids])
  # plink2 and 1000G output may carry extra FORMAT fields after GT.
  gt <- matrix(substr(gt, 1, 3), nrow = nrow(gt), dimnames = dimnames(gt))
  G <- t(gt_to_dosage(gt))
  colnames(G) <- paste(dt[["#CHROM"]], dt$POS, dt$REF, dt$ALT, sep = "_")
  rownames(G) <- ids
  G <- G[, apply(G, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  if (anyNA(G)) G[is.na(G)] <- 0
  G
}

gene_window <- function(g, size) {
  list(start = max(1, g$start - size), end = g$end + size)
}

# ---- prepare ---------------------------------------------------------------
step_prepare <- function() {
  set.seed(cfg$seed)
  ref <- cfg$reference[[cfg$truth_annotation]]
  genes <- read_gtf_genes(ref$gtf)
  autosomes <- paste0("chr", 1:22)
  genes <- genes[gene_type == "protein_coding" & sub("^chr", "", seqnames) %in% 1:22]
  genes[, chr := as.integer(sub("^chr", "", seqnames))]
  tx <- read_fasta_transcripts(ref$fasta)
  tx <- unique(tx[gene_id %in% genes$gene_id], by = "tx_id")
  counts <- tx[, .(n_isoforms = .N), by = gene_id]
  genes <- merge(genes, counts, by = "gene_id")
  genes <- genes[n_isoforms >= (cfg$simulation$min_isoforms %||% 1)]
  n <- cfg$simulation$n_genes
  if (!identical(n, "all")) genes <- genes[sample(.N, min(as.integer(n), .N))]
  setorder(genes, chr, start)

  s <- cfg$simulation
  genes[, prop_shared := sample(unlist(s$prop_shared), .N, replace = TRUE)]
  genes[, h2_g := if (is.list(s$h2_g)) rbeta(.N, s$h2_g$beta[[1]], s$h2_g$beta[[2]]) else s$h2_g]
  genes[, n_causal := s$n_causal]
  fwrite(genes, dir_for(out("params", "genes.tsv")), sep = "\t")
  fwrite(tx[gene_id %in% genes$gene_id], out("params", "transcripts.tsv"), sep = "\t")

  header <- system2("tabix", c("-H", shQuote(vcf_path(genes$chr[1]))), stdout = TRUE)
  vcf_ids <- strsplit(header[length(header)], "\t")[[1]][-(1:9)]
  ids <- if (!is.null(cfg$genotypes$sample_file)) {
    wanted <- readLines(cfg$genotypes$sample_file)
    intersect(wanted, vcf_ids)
  } else vcf_ids
  if (!is.null(cfg$genotypes$n_samples)) {
    n_s <- as.integer(cfg$genotypes$n_samples)
    if (n_s > length(ids)) stop("Requested ", n_s, " samples, VCF has ", length(ids))
    ids <- sort(sample(ids, n_s))
  }
  writeLines(ids, out("params", "samples.txt"))
  cat(sprintf("Prepared %d genes, %d transcripts, %d samples\n", nrow(genes), sum(genes$n_isoforms), length(ids)))
}

# ---- simulate: expression plus trait/GWAS per gene -------------------------
step_simulate <- function() {
  genes <- genes_table(); tx <- transcripts_table(); ids <- samples(); s <- cfg$simulation
  status <- character()
  for (i in chunk_rows(genes)) {
    g <- genes[i]
    path <- out("simulated", paste0(g$gene_id, ".rds"))
    if (file.exists(path)) {
      status[g$gene_id] <- readRDS(path)$status
      next
    }
    set.seed(cfg$seed + i)
    w <- gene_window(g, cfg$windows$train)
    G <- read_window(g$chr, w$start, w$end, ids)
    res <- list(gene_id = g$gene_id, status = "no_snps")
    if (!is.null(G)) {
      pos <- as.integer(sapply(strsplit(colnames(G), "_"), `[`, 2))
      in_sim <- which(pos >= g$start - cfg$windows$causal & pos <= g$end + cfg$windows$causal)
      sim <- if (length(in_sim)) simulate_gene(
        G[, in_sim, drop = FALSE], g$n_isoforms, g$n_causal, g$prop_shared, g$h2_g,
        unlist(s$cov_B), unlist(s$cov_E)) else NULL
      if (is.null(sim)) {
        res$status <- "too_few_snps"
      } else {
        sim$causal <- in_sim[sim$causal]
        tx_ids <- tx[gene_id == g$gene_id, tx_id]
        colnames(sim$Y) <- tx_ids
        X <- scale(G)
        gwas <- list()
        for (sc in cfg$trait$scenarios) {
          med_iso <- sample(g$n_isoforms, 1)
          sign <- sample(c(-1, 1), 1)
          beta <- trait_effects(G, sim, sc$h2ge, sc$mediator, isoform = med_iso, sign = sign)
          gwas[[sc$name]] <- list(
            z = simulate_gwas_z(X, beta, sc$n_gwas, sc$h2ge), h2ge = sc$h2ge,
            n_gwas = sc$n_gwas, mediator = sc$mediator, mediator_isoform = tx_ids[med_iso], sign = sign)
        }
        res <- c(sim, list(gene_id = g$gene_id, status = "ok", tx_ids = tx_ids,
                           causal_ids = colnames(G)[sim$causal], snp_ids = colnames(G), gwas = gwas))
      }
    }
    save_rds(res, path)
    status[g$gene_id] <- res$status
  }
  # One small table per chunk, so summarize can count the genes that were simulated.
  fwrite(data.table(gene_id = names(status), status = unname(status)),
         dir_for(out("simulated_status", sprintf("chunk_%05d.tsv", index))), sep = "\t")
}

# Status of every gene: "ok", "no_snps" or "too_few_snps".
simulation_status <- function() {
  files <- list.files(out("simulated_status"), pattern = "\\.tsv$", full.names = TRUE)
  if (!length(files)) stop("No simulation status found in ", out("simulated_status"),
                           "; rerun the simulate step")
  unique(rbindlist(lapply(files, fread)), by = "gene_id")
}

# ---- truth_matrix: combine simulated isoform expression, after all genes finish
step_truth_matrix <- function() {
  genes <- genes_table()
  paths <- out("simulated", paste0(genes$gene_id, ".rds"))
  absent <- paths[!file.exists(paths)]
  if (length(absent)) stop(length(absent), " simulated genes missing; rerun the simulate step")
  ys <- lapply(paths, function(p) readRDS(p)$Y)
  Y <- do.call(cbind, ys[!sapply(ys, is.null)])
  save_rds(Y, out("simulated_all", "isoform_log_expression.rds"))
}

# ---- reads: one sample per array index --------------------------------------
step_reads <- function() {
  suppressPackageStartupMessages(library(Rsubread))
  ids <- samples(); id <- ids[index]; r <- cfg$reads
  fa <- cfg$reference[[cfg$truth_annotation]]$fasta
  prefix <- out("reads", id)
  dir_for(prefix)
  Y <- readRDS(out("simulated_all", "isoform_log_expression.rds"))
  headers <- read_fasta_transcripts(fa)
  # simReads rejects duplicate transcript names, so keep each sequence once via
  # the header-level vector, in FASTA order.
  expr <- numeric(nrow(headers))
  m <- match(headers$tx_id, colnames(Y))
  expr[!is.na(m)] <- exp(Y[id, m[!is.na(m)]])
  expr[duplicated(headers$tx_id)] <- 0
  expr <- expr / sum(expr) * 1e6
  set.seed(cfg$seed + 1000000L + index)
  simReads(transcript.file = fa, expression.levels = expr, output.prefix = prefix,
           read.length = r$read_length, library.size = r$library_size,
           paired.end = isTRUE(r$paired_end))
}

# ---- aggregate: per-sample quantifications to TMM-normalized matrices --------
tmm_log2 <- function(counts) {
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts), method = "TMM")
  log2(edgeR::cpm(dge, normalized.lib.sizes = TRUE) + 1)
}

keep_expressed <- function(tpm) {
  rowMeans(tpm > 0.1) >= cfg$quant$min_prop_expressed
}

step_aggregate_salmon <- function(annot) {
  suppressPackageStartupMessages(library(tximport))
  ids <- samples()
  files <- out("salmon", annot, ids, "quant.sf")
  names(files) <- ids
  absent <- files[!file.exists(files)]
  if (length(absent)) stop(length(absent), " of ", length(files), " quant.sf files missing for ", annot)
  tx <- read_fasta_transcripts(cfg$reference[[annot]]$fasta)
  tx <- unique(tx, by = "header")
  # Salmon names transcripts by the full FASTA header, or by the transcript ID alone
  # when the index was built with --gencode. Map both forms.
  tx2gene <- unique(data.frame(header = c(tx$header, tx$tx_id),
                               gene = rep(strip_version(tx$gene_id), 2)))
  iso <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
  gene <- summarizeToGene(iso, tx2gene = tx2gene, countsFromAbundance = "no")
  rownames(iso$counts) <- rownames(iso$abundance) <- sub("\\|.*$", "", rownames(iso$counts))
  iso_keep <- keep_expressed(iso$abundance)
  gene_keep <- keep_expressed(gene$abundance)
  tx_map <- data.table(tx_id = sub("\\|.*$", "", tx$header), gene = strip_version(tx$gene_id))
  save_rds(list(expr = tmm_log2(iso$counts[iso_keep, , drop = FALSE]), tx_map = tx_map),
           out("quants", paste0("salmon_", annot, "_isoform.rds")))
  save_rds(list(expr = tmm_log2(gene$counts[gene_keep, , drop = FALSE])),
           out("quants", paste0("salmon_", annot, "_gene.rds")))
  write_bed(tmm_log2(gene$counts[gene_keep, , drop = FALSE]), annot, paste0("salmon_", annot))
}

step_featurecounts <- function(annot) {
  suppressPackageStartupMessages(library(Rsubread))
  ids <- samples()
  gtf <- read_gtf_genes(cfg$reference[[annot]]$gtf)
  anno <- data.frame(GeneID = gtf$gene_id, Chr = gtf$seqnames, Start = gtf$start, End = gtf$end, Strand = gtf$strand)
  bams <- out("star", paste0(ids, "_Aligned.sortedByCoord.out.bam"))
  fc <- featureCounts(bams, annot.ext = anno, isPairedEnd = isTRUE(cfg$reads$paired_end),
                      countReadPairs = isTRUE(cfg$reads$paired_end), nthreads = cfg$threads)
  counts <- fc$counts
  colnames(counts) <- ids
  rownames(counts) <- strip_version(rownames(counts))
  counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
  keep <- rowMeans(edgeR::cpm(counts) > 0.1) >= cfg$quant$min_prop_expressed
  expr <- tmm_log2(counts[keep, , drop = FALSE])
  save_rds(list(expr = expr), out("quants", paste0("star_", annot, "_gene.rds")))
  write_bed(expr, annot, paste0("star_", annot))
}

# QTL-ready phenotype BED: normalized expression, TSS coordinates, autosomes only.
write_bed <- function(expr, annot, label) {
  gtf <- read_gtf_genes(cfg$reference[[annot]]$gtf)
  gtf[, gene := strip_version(gene_id)]
  gtf <- unique(gtf[seqnames %in% paste0("chr", 1:22)], by = "gene")
  gtf <- gtf[gene %in% rownames(expr)]
  gtf[, tss := ifelse(strand == "+", start, end)]
  gtf[, chrnum := as.integer(sub("chr", "", seqnames))]
  setorder(gtf, chrnum, tss)
  bed <- data.table(`#Chr` = gtf$seqnames, start = gtf$tss - 1L, end = gtf$tss,
                    pid = gtf$gene, gid = gtf$gene, strand = gtf$strand,
                    expr[gtf$gene, , drop = FALSE])
  fwrite(bed, dir_for(out("bed", paste0(label, "_gene.bed"))), sep = "\t")
}

# ---- train_test: prediction models and TWAS for one expression source --------
cv_folds <- function(n, k) sample(rep_len(seq_len(k), n))

fit_method <- function(method, X, y) {
  if (method == "enet") {
    f <- glmnet::cv.glmnet(X, y, alpha = 0.5, nfolds = 5, standardize = FALSE)
    as.numeric(coef(f, s = "lambda.min"))[-1]
  } else if (method == "lasso") {
    f <- glmnet::cv.glmnet(X, y, alpha = 1, nfolds = 5, standardize = FALSE)
    as.numeric(coef(f, s = "lambda.min"))[-1]
  } else if (method == "blup") {
    as.numeric(rrBLUP::mixed.solve(y, Z = X)$u)
  } else if (method == "susie") {
    as.numeric(coef(susieR::susie(X, y, L = 10))[-1])
  } else stop("Unknown method: ", method)
}

# Choose among methods by out-of-sample R2, then refit on all samples.
train_best <- function(X, y, methods) {
  folds <- cv_folds(nrow(X), 5)
  scores <- sapply(methods, function(m) {
    pred <- rep(NA_real_, length(y))
    for (k in 1:5) {
      tr <- folds != k
      w <- tryCatch(fit_method(m, X[tr, , drop = FALSE], y[tr]), error = function(e) NULL)
      if (is.null(w)) return(NA_real_)
      pred[!tr] <- as.numeric(X[!tr, , drop = FALSE] %*% w)
    }
    if (var(pred) == 0) return(0)
    cor(pred, y)^2
  })
  if (all(is.na(scores))) return(NULL)
  best <- methods[which.max(scores)]
  list(method = best, cv_r2 = max(scores, na.rm = TRUE), weights = fit_method(best, X, y))
}

# isoTWAS models for one gene's isoforms, trained jointly with the isotwas package.
#   Y: samples x isoforms, standardized, with transcript IDs as column names.
# Returns one fit per isoform in the form train_best returns. The simulation has no
# replicate quantifications, so Y stands in for them and the isoform covariance is
# estimated from Y itself (omega_est "mean").
train_isotwas <- function(X, Y, opts, seed) {
  res <- isotwas::compute_isotwas(
    X = X, Y = Y, Y.rep = Y, R = 1, id = rownames(Y),
    omega_est = opts$omega_est %||% "mean", method = unlist(opts$methods),
    run_all = FALSE, nfolds = 5, verbose = FALSE, seed = seed)
  if (nrow(res$failed_methods)) cat(sprintf(
    "  isotwas method(s) failed: %s\n",
    paste(unique(paste0(res$failed_methods$method, " (", res$failed_methods$message, ")")),
          collapse = "; ")))
  sapply(colnames(Y), function(t) {
    m <- res$best_models$transcripts[[t]]
    w <- numeric(ncol(X))
    at <- match(m$weights$SNP, colnames(X))
    if (anyNA(at)) stop("isotwas returned SNPs that are not in the genotype matrix")
    w[at] <- m$weights$Weight
    list(method = m$method %||% NA_character_, cv_r2 = m$r2, weights = w)
  }, simplify = FALSE)
}

step_train_test <- function(source_name) {
  genes <- genes_table(); ids <- samples(); methods <- unlist(cfg$models$methods)
  isotwas_opts <- cfg$models$isotwas
  use_isotwas <- length(unlist(isotwas_opts$methods)) > 0
  if (use_isotwas) suppressPackageStartupMessages(library(isotwas))
  quant <- function(level) {
    p <- out("quants", paste0(source_name, "_", level, ".rds"))
    if (file.exists(p)) readRDS(p) else NULL
  }
  is_truth <- source_name == "truth"
  q_gene <- if (is_truth) NULL else quant("gene")
  q_iso <- if (is_truth) NULL else quant("isoform")

  results <- list()
  for (i in chunk_rows(genes)) {
    g <- genes[i]
    sim <- readRDS(out("simulated", paste0(g$gene_id, ".rds")))
    if (!identical(sim$status, "ok")) next
    set.seed(cfg$seed + 2000000L + i)
    w <- gene_window(g, cfg$windows$train)
    G <- read_window(g$chr, w$start, w$end, ids)
    X <- scale(G)
    gene_key <- strip_version(g$gene_id)

    features <- list()
    if (is_truth) {
      features[[g$gene_id]] <- list(level = "gene", y = gene_expression(sim$Y))
      # A single-isoform gene's isoform is its gene expression, so the isoform-level
      # test reuses the gene-level model. Both levels then cover the same genes.
      for (t in colnames(sim$Y))
        features[[t]] <- list(level = "isoform", y = sim$Y[, t],
                              same_as = if (g$n_isoforms == 1) g$gene_id)
    } else {
      if (!is.null(q_gene) && gene_key %in% rownames(q_gene$expr))
        features[[g$gene_id]] <- list(level = "gene", y = q_gene$expr[gene_key, ids])
      if (!is.null(q_iso)) {
        for (t in intersect(q_iso$tx_map[gene == gene_key, tx_id], rownames(q_iso$expr)))
          features[[t]] <- list(level = "isoform", y = q_iso$expr[t, ids])
      }
    }

    rows <- list()
    # One row per scenario: the TWAS test of this model against that scenario's GWAS.
    add_rows <- function(fit, level, name) {
      for (sc in names(sim$gwas)) {
        tw <- twas_test(fit$weights, sim$gwas[[sc]]$z, X)
        rows[[length(rows) + 1]] <<- data.table(
          gene_id = g$gene_id, n_isoforms = g$n_isoforms, h2_g = g$h2_g, gene_h2 = sim$gene_h2,
          source = source_name, level = level, feature = name, method = fit$method,
          cv_r2 = fit$cv_r2, scenario = sc, h2ge = sim$gwas[[sc]]$h2ge, n_gwas = sim$gwas[[sc]]$n_gwas,
          mediator = sim$gwas[[sc]]$mediator,
          is_mediator = level == "gene" || name == sim$gwas[[sc]]$mediator_isoform,
          twas_z = tw[["z"]], twas_p = tw[["p"]])
      }
    }
    fits <- list()
    for (name in names(features)) {
      f <- features[[name]]
      if (!is.null(f$same_as)) {
        fit <- fits[[f$same_as]]
      } else {
        y <- as.numeric(scale(f$y))
        fit <- if (all(is.finite(y))) train_best(X, y, methods)
      }
      if (is.null(fit)) next
      fits[[name]] <- fit
      add_rows(fit, f$level, name)
    }

    if (use_isotwas) {
      iso <- names(features)[vapply(features, function(f) f$level == "isoform", logical(1))]
      Y <- scale(do.call(cbind, lapply(features[iso], function(f) as.numeric(f$y))))
      dimnames(Y) <- list(ids, iso)
      Y <- Y[, colSums(!is.finite(Y)) == 0, drop = FALSE]
      iso_fits <- if (ncol(Y) >= 2) {
        tryCatch(train_isotwas(X, Y, isotwas_opts, cfg$seed + 3000000L + i), error = function(e) {
          cat(sprintf("  isotwas failed for %s: %s\n", g$gene_id, conditionMessage(e)))
          NULL
        })
      } else {
        # With one isoform there is nothing to model jointly: use its univariate model.
        fits[intersect(colnames(Y), names(fits))]
      }
      for (name in names(iso_fits)) add_rows(iso_fits[[name]], "isoform_isotwas", name)
    }
    if (length(rows)) results[[length(results) + 1]] <- rbindlist(rows)
  }
  path <- out("twas", source_name, sprintf("chunk_%05d.tsv.gz", index))
  dir_for(path)
  fwrite(if (length(results)) rbindlist(results) else data.table(), path, sep = "\t")
}

# ---- summarize: power, calibration and prediction accuracy -------------------
step_summarize <- function() {
  files <- list.files(out("twas"), pattern = "\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
  files <- files[file.size(files) > 40]
  if (!length(files)) stop("No TWAS results found in ", out("twas"))
  d <- rbindlist(lapply(files, fread), fill = TRUE)
  # A model with all-zero weights predicts nothing, so it is not a model. Dropping
  # these also keeps their p of 1 out of the isoform combination.
  null_models <- is_null_model(d$twas_z, d$twas_p)
  cat(sprintf("%d of %d trained models have all-zero weights and are not tested\n",
              uniqueN(d[null_models, .(source, level, feature)]),
              uniqueN(d[, .(source, level, feature)])))
  d <- d[!null_models]
  # The denominator is the genes that were simulated. Genes skipped for lack of
  # SNPs carry no signal, so they are left out and reported separately.
  status <- simulation_status()
  n_genes <- sum(status$status == "ok")
  n_skipped <- nrow(status) - n_genes
  if (n_skipped) cat(sprintf("%d of %d genes were not simulated (%s) and are excluded\n",
                             n_skipped, nrow(status),
                             paste(setdiff(unique(status$status), "ok"), collapse = ", ")))
  alpha <- cfg$summary$alpha

  gene_rows <- d[level == "gene", .(p = twas_p[1], cv_r2 = cv_r2[1]),
                 by = .(source, scenario, h2ge, n_gwas, mediator, gene_id, n_isoforms)]
  gene_rows[, level := "gene"]
  # Isoform-level gene tests: ACAT across the isoforms modeled for each gene, once
  # for the univariate isoform models and once for the isoTWAS models if trained.
  iso_rows <- d[level %in% c("isoform", "isoform_isotwas"),
                .(p = acat(twas_p), cv_r2 = mean(cv_r2)),
                by = .(source, scenario, h2ge, n_gwas, mediator, gene_id, n_isoforms, level)]
  iso_rows[, level := ifelse(level == "isoform", "isoform_acat", "isotwas_acat")]
  tests <- rbindlist(list(gene_rows, iso_rows), use.names = TRUE)
  summary <- tests[, .(
    n_genes_simulated = n_genes,
    n_genes_not_simulated = n_skipped,
    n_genes_modeled = uniqueN(gene_id),
    # Simulated genes with no model cannot be detected, so they count against power.
    rejection_rate = sum(p < alpha) / n_genes,
    rejection_rate_bonferroni = sum(p < alpha / n_genes) / n_genes,
    mean_cv_r2 = mean(cv_r2)
  ), by = .(source, level, scenario, h2ge, n_gwas, mediator)]
  setorder(summary, source, level, h2ge, n_gwas)
  fwrite(summary, dir_for(out("summary", "power_summary.tsv")), sep = "\t")
  fwrite(tests, out("summary", "gene_tests.tsv.gz"), sep = "\t")
  print(summary)
}

switch(step,
  prepare = step_prepare(),
  simulate = step_simulate(),
  truth_matrix = step_truth_matrix(),
  reads = step_reads(),
  aggregate_salmon = step_aggregate_salmon(args[4]),
  featurecounts = step_featurecounts(args[4]),
  train_test = step_train_test(args[4]),
  summarize = step_summarize(),
  stop("Unknown step: ", step))
