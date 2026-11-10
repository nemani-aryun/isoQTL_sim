# Statistical checks for R/sim_core.R. Run: Rscript tests/test_sim_core.R
# Each check compares simulation output with a value known from theory.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("--file=", "", args[grep("--file=", args)]))
source(file.path(here, "..", "R", "sim_core.R"))

failures <- 0L
check <- function(label, ok, detail = "") {
  cat(sprintf("[%s] %s %s\n", if (ok) "PASS" else "FAIL", label, detail))
  if (!ok) failures <<- failures + 1L
}

set.seed(20261104)

# Genotypes with LD: haplotypes from a latent AR(1) Gaussian, thresholded.
make_genotypes <- function(n = 500, p = 300, rho = 0.9) {
  hap <- function() {
    z <- matrix(rnorm(n * p), n, p)
    for (j in 2:p) z[, j] <- rho * z[, j - 1] + sqrt(1 - rho^2) * z[, j]
    (z > qnorm(0.7)) * 1
  }
  G <- hap() + hap()
  G <- G[, apply(G, 2, var) > 0]
  rownames(G) <- paste0("S", seq_len(n))
  G
}
G <- make_genotypes()
X <- scale(G)

# --- gt_to_dosage
gt <- matrix(c("0|0", "0|1", "1/1", "1|0"), 2)
check("gt_to_dosage", all(gt_to_dosage(gt) == c(0, 1, 2, 1)))

# --- isoform genetic variance hits the target exactly
errs <- unlist(lapply(c(1, 2, 3, 5, 8), function(M) {
  replicate(50, {
    h2 <- rbeta(1, 2, 20)
    sim <- simulate_gene(G, M, n_causal = 4, prop_shared = sample(c(.5, .75, 1), 1), h2_g = h2)
    Xc <- scale(G[, sim$causal, drop = FALSE])
    max(abs(apply(Xc %*% sim$B, 2, var) - h2) / h2)
  })
}))
check("isoform genetic variance equals target h2", max(errs) < 1e-8,
      sprintf("max relative error %.1e", max(errs)))

# --- total expression variance is about 1
vars <- replicate(100, {
  sim <- simulate_gene(G, 3, h2_g = 0.1)
  mean(apply(sim$Y, 2, var))
})
check("total isoform variance near 1", abs(mean(vars) - 1) < 0.02,
      sprintf("mean %.3f", mean(vars)))

# --- edge cases
one <- simulate_gene(G, 1, n_causal = 1, prop_shared = 1, h2_g = 0.1)
check("single isoform and single causal SNP", ncol(one$Y) == 1 && length(one$causal) == 1)
check("h2 = 0 gives zero effects", all(simulate_gene(G, 2, h2_g = 0)$B == 0))
check("too few SNPs returns NULL", is.null(simulate_gene(G[, 1:3], 5, n_causal = 4, prop_shared = 0.5)))

# --- trait effects explain h2ge
sim <- simulate_gene(G, 3, h2_g = 0.2)
for (med in c("gene", "isoform")) {
  b <- trait_effects(G, sim, 0.01, med, isoform = 2)
  check(sprintf("trait variance equals h2ge (%s)", med),
        abs(var(as.numeric(X %*% b)) - 0.01) < 1e-10)
}
check("h2ge = 0 gives zero trait effects", all(trait_effects(G, sim, 0) == 0))

# --- GWAS Z covariance equals LD
beta0 <- numeric(ncol(X))
zs <- replicate(4000, simulate_gwas_z(X, beta0, 1e5, 0))
R <- crossprod(X) / (nrow(X) - 1)
idx <- seq(1, ncol(X), length.out = 40)
check("null GWAS Z variance is 1", abs(mean(apply(zs, 1, var)) - 1) < 0.02)
check("null GWAS Z covariance matches LD",
      mean(abs(cov(t(zs))[idx, idx] - R[idx, idx])) < 0.04)

# --- TWAS with true weights: null calibration and analytic non-centrality
twas_rep <- function(h2ge, n_gwas, n = 2000) {
  replicate(n, {
    s <- simulate_gene(G, 1, n_causal = 4, prop_shared = 1, h2_g = 0.1)
    w <- numeric(ncol(X)); w[s$causal] <- s$B[, 1]
    b <- trait_effects(G, s, h2ge, "gene")
    twas_test(w, simulate_gwas_z(X, b, n_gwas, h2ge), X)["z"]
  })
}
null_z <- twas_rep(0, 1e5)
check("null TWAS Z mean 0", abs(mean(null_z)) < 0.08, sprintf("mean %.3f", mean(null_z)))
check("null TWAS Z variance 1", abs(var(null_z) - 1) < 0.1, sprintf("var %.3f", var(null_z)))
alt_z <- twas_rep(0.001, 1e5)
expected <- sqrt(1e5 * 0.001 / (1 - 0.001))
check("alternative TWAS Z matches sqrt(N * h2ge / (1 - h2ge))",
      abs(mean(alt_z) - expected) < 0.1,
      sprintf("mean %.3f expected %.3f", mean(alt_z), expected))

# --- zero weights
check("zero weights give p = 1", twas_test(numeric(ncol(X)), rnorm(ncol(X)), X)["p"] == 1)

# --- models with no predicted variance are flagged and kept out of ACAT
null_test <- twas_test(numeric(ncol(X)), rnorm(ncol(X)), X)
check("zero-weight model is flagged as a null model", is_null_model(null_test["z"], null_test["p"]))
check("real tests are not flagged", !any(is_null_model(c(0.5, -3, 0), c(0.62, 0.003, 0.999))))
p_mixed <- c(1e-6, 0.5, 1); z_mixed <- c(4.9, 0.67, 0)
check("ACAT ignores null models", abs(acat(p_mixed[!is_null_model(z_mixed, p_mixed)]) - acat(c(1e-6, 0.5))) < 1e-12 &&
        acat(c(1e-6, 0.5)) < 1e-5)

# --- ACAT
check("ACAT of identical p returns p", abs(acat(c(0.05, 0.05, 0.05)) - 0.05) < 1e-8)
p_null <- replicate(5000, acat(runif(5)))
check("ACAT calibrated under independence", abs(mean(p_null < 0.05) - 0.05) < 0.01,
      sprintf("rate %.3f", mean(p_null < 0.05)))

if (failures) {
  cat(failures, "check(s) failed\n")
  quit(status = 1)
}
cat("All checks passed\n")
