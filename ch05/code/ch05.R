# =============================================================================
#  Chapter 5 -- Hypothesis testing
#  Reproduces every number and figure in ch05/ch05.tex.
#
#  Run from the repository root:   Rscript ch05/code/ch05.R
#  Needs base R + MASS (Boston data).
# =============================================================================

source(file.path("common", "helpers.R"))
set.seed(2026)
FIG_DIR <- file.path("ch05", "figures")
logsumexp <- function(x) max(x) + log(sum(exp(x - max(x))))


# ---- Exercise 1: one-sided test for a Poisson rate --------------------------------

section("Ex 1  Poisson rate, prior median 1")

# If lambda ~ Gamma(0.1, b) then b * lambda ~ Gamma(0.1, 1), so the prior median is m1 / b.
b1 <- qgamma(0.5, 0.1, 1)                    # makes the prior median equal to 1
say("b = %.6g   check P(lambda <= 1) = %.4f", b1, pgamma(1, 0.1, b1))
cases1 <- data.frame(N = c(10, 20, 50, 100), Y = c(12, 24, 60, 120))
cases1$post_alt <- pgamma(1, 0.1 + cases1$Y, b1 + cases1$N, lower.tail = FALSE)
cases1$BF <- cases1$post_alt / (1 - cases1$post_alt)          # prior odds are 1
print(transform(cases1, post_alt = round(post_alt, 4), BF = signif(BF, 4)), row.names = FALSE)


# ---- Exercise 2: Moderna trial, which arm has the lower rate? -------------------

section("Ex 2  Moderna: Bayes factor for the direction of the effect")

yP <- 185; nP <- 14073; yV <- 11; nV <- 14134

# log P(theta_V > theta_P | data) = log of integral of f_P(x) * P(theta_V > x) dx, on a fine grid.
log_prob_V_above_P <- function(a0, b0) {
  aP <- a0 + yP; bP <- b0 + nP - yP; aV <- a0 + yV; bV <- b0 + nV - yV
  x  <- seq(1e-6, 0.06, length.out = 200001); dx <- diff(x)[1]
  lf <- dbeta(x, aP, bP, log = TRUE) +
    suppressWarnings(pbeta(x, aV, bV, lower.tail = FALSE, log.p = TRUE))   # extreme tail: ~3 digits
  logsumexp(lf) + log(dx)
}
priors2 <- list("Beta(1,1) uniform" = c(1, 1), "Beta(0.5,0.5) Jeffreys" = c(0.5, 0.5),
                "Beta(1,99) mean 1%" = c(1, 99), "Beta(10,990) mean 1%, strong" = c(10, 990))
for (nm in names(priors2)) {
  lp <- log_prob_V_above_P(priors2[[nm]][1], priors2[[nm]][2])
  say("%-30s  P(H1: placebo lower | data) = %.3e   log10 BF(H2:H1) = %.1f", nm, exp(lp),
      (log1p(-exp(lp)) - lp) / log(10))
}
# Independent check (uniform prior): for small rates Beta(a, n) ~ Gamma(a, 1) / n, so
# P(theta_V > theta_P) ~ P(G_V / (G_V + G_P) > nV / (nV + nP)) with G_V/(G_V+G_P) ~ Beta(aV, aP).
say("gamma-approximation check (uniform prior): %.3e",
    pbeta(nV / (nV + nP), 1 + yV, 1 + yP, lower.tail = FALSE))


# ---- Exercise 3: Mr. October, one rate or two? -----------------------------------

section("Ex 3  separate vs common home-run rate")

Y3 <- c(563, 10); N3 <- c(2820, 27)
# log of integral_0^c Poisson(Y; N lambda) d lambda = -log N + log pgamma(c, Y + 1, N)
log_m_sep <- function(c) sum(-log(c) - log(N3) + pgamma(c, Y3 + 1, N3, log.p = TRUE))
log_m_com <- function(c) {
  S <- sum(Y3); M <- sum(N3)
  -log(c) + sum(Y3 * log(N3) - lgamma(Y3 + 1)) + lgamma(S + 1) - (S + 1) * log(M) +
    pgamma(c, S + 1, M, log.p = TRUE)
}
for (c3 in c(1, 10)) {
  lbf <- log_m_sep(c3) - log_m_com(c3)
  say("c = %2d   BF(H1 separate : H2 common) = %.3f   P(H1 | data) = %.3f", c3, exp(lbf), plogis(lbf))
}


# ---- Exercise 4: deforestation strategies ----------------------------------------

section("Ex 4  deforestation strategies")

s1 <- c(-11.5, -8.4, -6.8, -5.6, -0.4, 1.9, 2.8)
s2 <- c(-7.0, -4.3, -3.0, -3.0, -0.1, 0.0, 2.2, 4.3, 10.8)
say("means %.3f %.3f   sds %.3f %.3f", mean(s1), mean(s2), sd(s1), sd(s2))

# (a) separate Jeffreys priors: mu_j | Y ~ t_{n_j - 1}(ybar_j, s_j^2 / n_j)
S4 <- 1e6
mu1 <- mean(s1) + sd(s1) / sqrt(7) * rt(S4, 6)
mu2 <- mean(s2) + sd(s2) / sqrt(9) * rt(S4, 8)
say("(a) P(mu1 < mu2 | Y) = %.4f   diff mu2 - mu1: mean %.3f  95%% CI %s",
    mean(mu1 < mu2), mean(mu2 - mu1), fmt_ci(quantile(mu2 - mu1, c(0.025, 0.975))))

# (b) regression with Jeffreys prior: beta2 | Y ~ t_{14}(bhat2, se^2)
y4 <- c(s1, s2); x4 <- rep(0:1, c(7, 9))
fit4 <- summary(lm(y4 ~ x4))
b2 <- fit4$coef[2, 1]; se2 <- fit4$coef[2, 2]
say("(b) beta2 = %.3f  se %.3f   95%% CI %s   P(beta2 > 0) = %.4f",
    b2, se2, fmt_ci(b2 + qt(c(0.025, 0.975), 14) * se2), pt(b2 / se2, 14))

# (c) beta2 = gamma * delta, delta ~ N(0, W sigma^2): closed-form Bayes factor
#     (flat prior on beta1, p(sigma^2) ∝ 1/sigma^2).
bf_point_null <- function(y, x, W) {
  yc <- y - mean(y); xc <- x - mean(x); n <- length(y)
  S1 <- sum(yc^2)
  S2 <- S1 - sum(xc * yc)^2 / (sum(xc^2) + 1 / W)
  (1 + W * sum(xc^2))^(-1 / 2) * (S2 / S1)^(-(n - 1) / 2)
}
for (W in c(0.1, 1, 10, 100)) {
  bf <- bf_point_null(y4, x4, W)
  say("(c) W = %5.1f   BF(H2 : H1) = %.4f   P(H2: mu1 != mu2 | Y) = %.4f", W, bf, bf / (1 + bf))
}


# ---- Exercise 5: Bayes factor for a normal mean (equation 5.10) ------------------

section("Ex 5  Bayes factor for a normal mean")

bf_normal <- function(y, tau) (1 + tau^2)^(-1 / 2) * exp(y^2 / 2 * tau^2 / (1 + tau^2))
threshold <- function(tau, k) sqrt((1 + tau^2) / tau^2 * 2 * (log(k) + 0.5 * log(1 + tau^2)))
for (tau in c(1, 10)) say("tau = %2d: BF > 1 when |y| > %.3f;  BF > 10 when |y| > %.3f;  BF > 100 when |y| > %.3f;  BF at y = 0: %.3f",
                          tau, threshold(tau, 1), threshold(tau, 10), threshold(tau, 100), bf_normal(0, tau))

save_fig("ex05_bf", width = 6.5, height = 4, function() {
  y <- seq(-5, 5, length.out = 1001)
  plot(y, log10(bf_normal(y, 1)), type = "l", lwd = 2.5, col = PALETTE["navy"], ylim = c(-1.5, 4),
       xlab = "observed y", ylab = expression(log[10] ~ BF(y)))
  lines(y, log10(bf_normal(y, 10)), lwd = 2.5, col = PALETTE["red"])
  abline(h = c(0, 1), lty = c(2, 3), col = PALETTE["grey"])
  legend("top", c(expression(tau == 1), expression(tau == 10)), lwd = 2.5,
         col = PALETTE[c("navy", "red")], bty = "n", horiz = TRUE)
})


# ---- Exercise 6: A/B/C test for app designs --------------------------------------

section("Ex 6  A/B/C test")

counts6 <- rbind(A = c(37, 40, 15, 7, 1, 0, 0), B = c(12, 29, 28, 16, 8, 6, 1), C = c(27, 49, 18, 5, 0, 1, 0))
totals6 <- drop(counts6 %*% 0:6)
say("total opens: A %d  B %d  C %d  (each n = 100)", totals6["A"], totals6["B"], totals6["C"])
shape6 <- totals6 + 0.5; rate6 <- 100                         # Jeffreys prior
draws6 <- sapply(names(totals6), function(g) rgamma(1e6, shape6[g], rate6))
print(round(t(apply(draws6, 2, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975))))), 3))
best6 <- table(factor(colnames(draws6)[max.col(draws6)], levels = names(totals6))) / nrow(draws6)
print(round(best6, 5))
say("dispersion check (variance/mean): %s",
    paste(sprintf("%s %.2f", names(totals6), apply(counts6, 1, function(cn) {
      v <- rep(0:6, cn); var(v) / mean(v) })), collapse = "  "))

save_fig("ex06_rates", width = 6.5, height = 3.8, function() {
  lam <- seq(0.6, 2.6, length.out = 800)
  plot(NA, xlim = range(lam), ylim = c(0, 4.5), xlab = expression("mean opens per month " * lambda),
       ylab = "Posterior density")
  for (k in 1:3) lines(lam, dgamma(lam, shape6[k], rate6), lwd = 2.5, col = PALETTE[k])
  legend("topright", paste("design", names(totals6)), lwd = 2.5, col = PALETTE[1:3], bty = "n")
})


# ---- Exercise 7: thyroid cancer rates and the Bayesian FDR -----------------------

section("Ex 7  thyroid cancer, Bayesian FDR")

N7 <- c(179611, 125453, 41646, 98826, 68642, 71844, 67702, 45954, 96325, 200201,
        30557, 81858, 103475, 166203, 52032, 85767, 145591, 139176, 87842, 97930)
Y7 <- c(20, 18, 7, 22, 10, 6, 9, 9, 18, 27, 5, 15, 12, 14, 12, 15, 20, 15, 24, 11)
lambda0 <- 13 / 1e5

# (a) Gamma(a, b): the ratio q95 / q50 depends only on a; solve it = 10, then set b.
a7 <- uniroot(function(a) qgamma(0.95, a, 1) / qgamma(0.5, a, 1) - 10, c(0.05, 5))$root
b7 <- qgamma(0.5, a7, 1) / lambda0
say("a = %.4f  b = %.1f   check: P(lambda > lambda0) = %.3f   q95 / lambda0 = %.2f",
    a7, b7, pgamma(lambda0, a7, b7, lower.tail = FALSE), qgamma(0.95, a7, b7) / lambda0)

q7 <- pgamma(lambda0, a7 + Y7, b7 + N7, lower.tail = FALSE)
tab7 <- data.frame(region = 1:20, N = N7, Y = Y7, rate_per_100k = round(1e5 * Y7 / N7, 1),
                   post_mean_per_100k = round(1e5 * (a7 + Y7) / (b7 + N7), 1), q = round(q7, 4))
print(tab7[order(-q7), ], row.names = FALSE)

# (b) Bayesian FDR of the rule q_i > t is the average of (1 - q_i) over the regions flagged.
fdr_at <- function(t) { flag <- q7 > t; if (!any(flag)) 0 else mean(1 - q7[flag]) }
cands <- sort(unique(c(0, q7)))
fdrs  <- sapply(cands, fdr_at)
t7 <- min(cands[fdrs <= 0.10])
say("chosen t = %.4f (smallest threshold with BFDR <= 0.10)   BFDR = %.4f   regions flagged: %s",
    t7, fdr_at(t7), paste(which(q7 > t7), collapse = ", "))
print(data.frame(t = round(cands, 4), n_flagged = sapply(cands, function(t) sum(q7 > t)),
                 BFDR = round(fdrs, 4))[cands > 0.5, ], row.names = FALSE)


# ---- Exercise 8: frequentist properties of a Bayesian test ------------------------

section("Ex 8  type I error and power")

# Flat prior: lambda | Y ~ Gamma(S + 1, n), S = sum(Y) ~ Poisson(n * lambda_true).
# Reject when P(lambda > 1 | Y) > 0.95, i.e. when S >= s_crit(n).
s_crit <- function(n) min(which(pgamma(1, 0:(20 * n) + 1, n, lower.tail = FALSE) > 0.95)) - 1
reject_prob <- function(lambda, n) ppois(s_crit(n) - 1, n * lambda, lower.tail = FALSE)

lam8 <- seq(1, 2, by = 0.05)
pow10 <- sapply(lam8, reject_prob, n = 10)
say("n = 10: critical S = %d", s_crit(10))
print(round(data.frame(lambda = lam8, P_reject = pow10)[seq(1, 21, by = 2), ], 4), row.names = FALSE)
n8 <- 10:100
type1 <- sapply(n8, function(n) reject_prob(1, n))
power15 <- sapply(n8, function(n) reject_prob(1.5, n))
say("type I error over n = 10..100: min %.4f  max %.4f  mean %.4f", min(type1), max(type1), mean(type1))
say("power at lambda = 1.5: n = 10 %.3f, n = 20 %.3f, n = 50 %.3f, n = 100 %.4f",
    power15[n8 == 10], power15[n8 == 20], power15[n8 == 50], power15[n8 == 100])

save_fig("ex08_power", width = 7.5, height = 3.2, function() {
  par(mfrow = c(1, 3), mar = c(4, 4.2, 2, 0.5))
  plot(lam8, pow10, type = "b", pch = 19, cex = 0.6, col = PALETTE["navy"], ylim = c(0, 1),
       xlab = expression("true " * lambda), ylab = "P(reject)", main = "(a) n = 10")
  abline(h = 0.05, lty = 2, col = PALETTE["red"])
  plot(n8, type1, type = "l", lwd = 2, col = PALETTE["navy"], ylim = c(0, 0.1),
       xlab = "n", ylab = "type I error", main = expression("(b) " * lambda == 1))
  abline(h = 0.05, lty = 2, col = PALETTE["red"])
  plot(n8, power15, type = "l", lwd = 2, col = PALETTE["navy"], ylim = c(0, 1),
       xlab = "n", ylab = "power", main = expression("(c) " * lambda == 1.5))
})


# ---- Exercise 9: SSVS for the Boston housing data --------------------------------

section("Ex 9  SSVS for Boston")

boston <- MASS::Boston
Xs <- scale(as.matrix(boston[, setdiff(names(boston), "medv")]))
yc <- boston$medv - mean(boston$medv)
n9 <- nrow(Xs); p9 <- ncol(Xs); W9 <- 1
XtX9 <- crossprod(Xs); Xty9 <- crossprod(Xs, yc); yty9 <- sum(yc^2)

# log marginal likelihood of model gamma, integrating out beta_gamma ~ N(0, W sigma^2 I),
# the intercept (flat) and sigma^2 (prior ∝ 1/sigma^2); constants common to all models dropped.
log_marg <- function(gam) {
  k <- sum(gam)
  if (k == 0) return(-(n9 - 1) / 2 * log(yty9))
  A <- XtX9[gam, gam, drop = FALSE] + diag(1 / W9, k)
  R <- chol(A); v <- forwardsolve(t(R), Xty9[gam])
  -0.5 * k * log(W9) - sum(log(diag(R))) - (n9 - 1) / 2 * log(yty9 - sum(v^2))
}

# Exact posterior over all 2^13 models (prior: gamma_j iid Bernoulli(0.5)).
all_models <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), p9)))
lm_all <- apply(all_models, 1, log_marg)
post_models <- exp(lm_all - logsumexp(lm_all))
incl_exact <- colSums(all_models * post_models)
top <- order(-post_models)[1:5]

# SSVS Gibbs sampler: update each gamma_j from its full conditional (beta integrated out).
ssvs <- function(iters = 6000, burn = 1000) {
  gam <- rep(TRUE, p9); out <- matrix(NA, iters, p9)
  for (it in seq_len(iters)) {
    for (j in sample(p9)) {
      g1 <- gam; g1[j] <- TRUE; g0 <- gam; g0[j] <- FALSE
      gam[j] <- runif(1) < plogis(log_marg(g1) - log_marg(g0))
    }
    out[it, ] <- gam
  }
  out[-seq_len(burn), ]
}
draws9 <- ssvs()
incl_ssvs <- colMeans(draws9)
keys <- apply(draws9, 1, function(g) paste(colnames(Xs)[g], collapse = " + "))
freq <- sort(table(keys), decreasing = TRUE)
print(data.frame(variable = colnames(Xs), SSVS = round(incl_ssvs, 3), exact = round(incl_exact, 3)), row.names = FALSE)
cat("Most frequent models in the SSVS chain:\n"); print(round(head(freq, 3) / nrow(draws9), 3))
cat("Top models by exact enumeration:\n")
for (i in top) say("  %.3f   %s", post_models[i], paste(colnames(Xs)[all_models[i, ]], collapse = " + "))

save_fig("ex09_inclusion", width = 6.5, height = 4, function() {
  o <- order(incl_exact)
  par(mar = c(4.2, 6, 1, 1))
  barplot(incl_exact[o], names.arg = colnames(Xs)[o], horiz = TRUE, las = 1, col = PALETTE["navy"],
          border = NA, xlim = c(0, 1), xlab = "posterior inclusion probability")
  points(incl_ssvs[o], seq(0.7, by = 1.2, length.out = p9), pch = 19, col = PALETTE["red"])
  abline(v = 0.5, lty = 2, col = PALETTE["grey"])
  legend("bottomright", c("exact (all 8,192 models)", "SSVS sampler"), fill = c(PALETTE["navy"], NA),
         border = NA, pch = c(NA, 19), col = c(NA, PALETTE["red"]), bty = "n")
})
