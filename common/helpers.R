# =============================================================================
#  Shared helpers for every chapter script (base R only).
#  Each chapter script calls  source("common/helpers.R")  and sets FIG_DIR.
# =============================================================================

PALETTE <- c(navy = "#1F4E79", red = "#C0392B", green = "#2E8B57", gold = "#D68910",
             purple = "#7D3C98", teal = "#17A589", grey = "#5D6D7E")

# Open a PDF device with a clean, consistent theme, draw, and close it.
save_fig <- function(name, draw, width = 6.5, height = 4) {
  pdf(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)
  on.exit(invisible(dev.off()))
  par(mar = c(4.2, 4.4, 1.2, 1), las = 1, bty = "l", cex.lab = 1.1,
      col.axis = "#34495E", fg = "#5D6D7E", tcl = -0.3, mgp = c(2.6, 0.6, 0))
  draw()
}

section <- function(title) cat(sprintf("\n-- %s %s\n", title, strrep("-", 70 - nchar(title))))
say     <- function(fmt, ...) cat(sprintf(paste0(fmt, "\n"), ...))

# Equal-tailed 95% interval of a distribution given its quantile function.
ci95 <- function(qfun, ...) qfun(c(0.025, 0.975), ...)
fmt_ci <- function(ci, digits = 3) sprintf(paste0("(%.", digits, "f, %.", digits, "f)"), ci[1], ci[2])


# ---- MCMC diagnostics (no packages needed) -----------------------------------

# Effective sample size via Geyer's initial positive sequence of autocorrelations.
ess <- function(x) {
  x <- as.numeric(x)
  if (var(x) == 0) return(NA_real_)
  rho   <- acf(x, lag.max = min(length(x) - 1, 5000), plot = FALSE)$acf[, 1, 1]
  pairs <- rho[seq(1, length(rho) - 1, by = 2)] + rho[seq(2, length(rho), by = 2)]
  keep  <- if (any(pairs <= 0)) seq_len(which(pairs <= 0)[1] - 1) else seq_along(pairs)
  length(x) / max(1, -1 + 2 * sum(pairs[keep]))
}

# Gelman-Rubin potential scale reduction factor; `chains` is iterations x chains.
rhat <- function(chains) {
  n <- nrow(chains)
  W <- mean(apply(chains, 2, var))
  B <- n * var(colMeans(chains))
  sqrt(((n - 1) / n * W + B / n) / W)
}

# Overlaid trace plot of several chains (iterations x chains matrix).
trace_plot <- function(chains, ylab, main = "") {
  matplot(chains, type = "l", lty = 1, lwd = 0.6, col = adjustcolor(PALETTE, 0.75),
          xlab = "iteration", ylab = ylab, main = main)
}

# Posterior summary row: mean, sd and 95% interval of a vector of draws.
post_summary <- function(x, digits = 4) {
  round(c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))), digits)
}


# ---- Reusable samplers ----------------------------------------------------------

# Gibbs sampler for Bayesian linear regression
#   y ~ N(X beta, sigma^2),  beta_j ~ N(0, prior_sd^2),  sigma^2 ~ InvGamma(a, b).
# Returns a matrix of draws (columns: coefficients, then "sigma").
gibbs_lm <- function(y, X, prior_sd = 100, a = 0.01, b = 0.01, iters = 10000, burn = 1000) {
  n <- length(y); p <- ncol(X)
  XtX <- crossprod(X); Xty <- crossprod(X, y)
  beta <- qr.solve(X, y); sigma2 <- mean((y - X %*% beta)^2)
  coef_names <- if (is.null(colnames(X))) paste0("beta", seq_len(p) - 1) else colnames(X)
  draws <- matrix(NA, iters, p + 1, dimnames = list(NULL, c(coef_names, "sigma")))
  for (it in seq_len(iters)) {
    prec  <- XtX / sigma2 + diag(1 / prior_sd^2, p)
    R     <- chol(prec)                                   # prec = R'R
    mean  <- backsolve(R, forwardsolve(t(R), Xty / sigma2))
    beta  <- mean + backsolve(R, rnorm(p))                # draw from N(mean, prec^-1)
    sse   <- sum((y - X %*% beta)^2)
    sigma2 <- 1 / rgamma(1, a + n / 2, b + sse / 2)
    draws[it, ] <- c(beta, sqrt(sigma2))
  }
  draws[-seq_len(burn), , drop = FALSE]
}

# Random-walk Metropolis with a multivariate normal proposal (block update).
#   log_post: function(theta) returning the log posterior up to a constant
#   prop_cov: proposal covariance, e.g. (2.4^2 / p) * vcov(glm_fit)
metropolis_block <- function(log_post, init, prop_cov, iters = 20000, burn = 2000) {
  p <- length(init); L <- t(chol(prop_cov))
  theta <- init; lp <- log_post(theta); accepted <- 0
  draws <- matrix(NA, iters, p, dimnames = list(NULL, names(init)))
  for (it in seq_len(iters)) {
    cand <- theta + drop(L %*% rnorm(p))
    lp_cand <- log_post(cand)
    if (log(runif(1)) < lp_cand - lp) { theta <- cand; lp <- lp_cand; accepted <- accepted + 1 }
    draws[it, ] <- theta
  }
  structure(draws[-seq_len(burn), , drop = FALSE], acceptance = accepted / iters)
}

# Inverse-Gaussian random numbers (Michael, Schucany & Haas, 1976).
rinvgauss <- function(n, mean, shape) {
  nu <- rnorm(n)^2
  x  <- mean + mean^2 * nu / (2 * shape) - mean / (2 * shape) * sqrt(4 * mean * shape * nu + mean^2 * nu^2)
  ifelse(runif(n) <= mean / (mean + x), x, mean^2 / x)
}

# Summary table (mean, sd, 95% interval) for every column of a draws matrix.
summary_table <- function(draws, digits = 3) {
  t(apply(draws, 2, function(x) round(c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))), digits)))
}
