# Core simulation model for isoQTL_sim. Pure functions: no paths, no HPC state.
#
# Differences from the original s01_s02_sim_modified.R:
#  * Effect scaling uses sweep(B, 2, ...). The original `B * scaling_factors`
#    recycled along rows, so realized isoform h2 missed its target.
#  * Genes whose window has too few SNPs return NULL instead of aborting a batch.
#  * Trait and GWAS simulation are new (simulate_gwas_z).

# Convert a VCF GT character matrix ("0|1", "1/1") to 0/1/2 dosages.
gt_to_dosage <- function(gt) {
  dosage <- matrix(as.numeric(substr(gt, 1, 1)) + as.numeric(substr(gt, 3, 3)),
                   nrow = nrow(gt), dimnames = dimnames(gt))
  dosage
}

random_correlation <- function(M, low, high) {
  S <- diag(1, M)
  if (M > 1) {
    S[upper.tri(S)] <- runif(M * (M - 1) / 2, low, high)
    S[lower.tri(S)] <- t(S)[lower.tri(S)]
  }
  S
}

ensure_pd <- function(S) {
  if (!corpcor::is.positive.definite(S)) S <- corpcor::make.positive.definite(S)
  S
}

# Simulate isoform expression for one gene.
#   G: samples x SNPs dosage matrix for the simulation window.
#   M: number of isoforms.
# Each isoform has `n_causal` causal SNPs: ceiling(prop_shared * n_causal) shared
# across isoforms with correlated effects, and the rest specific to that isoform.
simulate_gene <- function(G, M, n_causal = 4, prop_shared = 0.75, h2_g = 0.1,
                          cov_B = c(0.7, 0.9), cov_E = c(0, 0.5)) {
  stopifnot(M >= 1, h2_g >= 0, h2_g < 1)
  S <- ceiling(prop_shared * n_causal)
  D <- n_causal - S
  W <- S + D * M
  if (W > ncol(G)) return(NULL)

  causal <- sample(ncol(G), W)
  X <- scale(G[, causal, drop = FALSE])
  N <- nrow(X)

  B <- matrix(0, nrow = W, ncol = M)
  if (S > 0) {
    Sigma_B <- ensure_pd(cov2cor(random_correlation(M, cov_B[1], cov_B[2])))
    B[seq_len(S), ] <- matrix(MASS::mvrnorm(S, rep(0, M), Sigma_B), nrow = S, ncol = M)
  }
  if (D > 0) {
    for (m in seq_len(M)) {
      B[S + (m - 1) * D + seq_len(D), m] <- rnorm(D)
    }
  }

  # Scale each isoform (column) so its in-sample genetic variance equals h2_g.
  genetic_var <- apply(X %*% B, 2, var)
  B <- sweep(B, 2, sqrt(h2_g / genetic_var), `*`)
  B[!is.finite(B)] <- 0
  colnames(B) <- NULL

  Sigma_E <- random_correlation(M, cov_E[1], cov_E[2]) * (1 - h2_g)
  Sigma_E <- ensure_pd(Sigma_E)
  E <- matrix(MASS::mvrnorm(N, rep(0, M), Sigma_E), nrow = N, ncol = M)

  genetic <- X %*% B
  Y <- genetic + E
  rownames(Y) <- rownames(G)

  one <- rep(1, M)
  list(
    Y = Y,
    B = B,
    causal = causal,
    isoform_h2 = unname(apply(genetic, 2, var) / apply(Y, 2, var)),
    gene_h2 = as.numeric(crossprod(one, cov(genetic) %*% one) /
                           crossprod(one, cov(Y) %*% one))
  )
}

# Gene-level expression on the same log scale as the isoforms.
gene_expression <- function(Y) log(rowSums(exp(Y)))

# Per-SNP (standardized genotype) trait effects from a simulated gene.
#   mediator = "gene":    the trait acts through the sum of isoform genetic values.
#   mediator = "isoform": the trait acts through one isoform (`isoform`).
# Effects are scaled so the mediated component explains h2ge of trait variance.
trait_effects <- function(G, sim, h2ge, mediator = c("gene", "isoform"),
                          isoform = 1L, sign = 1) {
  mediator <- match.arg(mediator)
  beta <- numeric(ncol(G))
  if (h2ge <= 0) return(beta)
  weights <- if (mediator == "gene") rep(1, ncol(sim$B)) else
    as.numeric(seq_len(ncol(sim$B)) == isoform)
  b <- as.numeric(sim$B %*% weights)
  X <- scale(G[, sim$causal, drop = FALSE])
  s2 <- var(as.numeric(X %*% b))
  if (s2 <= 0) return(beta)
  beta[sim$causal] <- sign * b * sqrt(h2ge / s2)
  beta
}

# Simulate GWAS marginal Z-scores for a trait with effects `beta`.
# Z ~ N(sqrt(Ngwas / (1 - h2ge)) * R %*% beta, R), where R is the LD of the real
# panel in X. Noise is drawn as X'e / sqrt(n - 1), whose covariance is exactly R,
# so no ridge term is needed even though p > n.
simulate_gwas_z <- function(X, beta, n_gwas, h2ge) {
  n <- nrow(X)
  mean_z <- sqrt(n_gwas / (1 - h2ge)) *
    as.numeric(crossprod(X, X %*% beta)) / (n - 1)
  noise <- as.numeric(crossprod(X, rnorm(n))) / sqrt(n - 1)
  mean_z + noise
}

# TWAS statistic from weights, GWAS Z and panel genotypes X (standardized).
# Returns Z = 0 and p = 1 when the weights carry no predicted variance.
twas_test <- function(weights, z, X) {
  pred <- as.numeric(X %*% weights)
  v <- sum(pred^2) / (nrow(X) - 1)
  if (!is.finite(v) || v <= 0) return(c(z = 0, p = 1))
  tz <- sum(weights * z) / sqrt(v)
  c(z = tz, p = 2 * pnorm(-abs(tz)))
}

# TRUE for TWAS results from a model with no predicted variance, which twas_test
# reports as Z = 0 and p = 1. These are not tests and must be dropped before
# combining: a p of 1 outweighs every other p-value in a Cauchy combination.
is_null_model <- function(z, p) z == 0 & p == 1

# Cauchy combination of p-values, used to combine isoform-level TWAS tests.
acat <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(1)
  p <- pmin(pmax(p, 1e-300), 1 - 1e-15)
  stat <- mean(tan((0.5 - p) * pi))
  pcauchy(stat, lower.tail = FALSE)
}
