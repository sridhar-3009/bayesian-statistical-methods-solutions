# =============================================================================
#  Chapter 9 -- Statistical properties of Bayesian methods
#  Reproduces every number and figure in ch09/ch09.tex.
#
#  Run from the repository root:   Rscript ch09/code/ch09.R
#  Needs base R only.
# =============================================================================

source(file.path("common", "helpers.R"))
set.seed(2026)
FIG_DIR <- file.path("ch09", "figures")


# ---- Exercise 1: decision-theoretic threshold for a one-sided test ---------------------

section("Ex 1  choosing the threshold c")

# Y | mu ~ N(mu, 2), mu ~ N(0, 2)  =>  mu | Y ~ N(Y / 2, 1)  and  P(H1 | Y) = pnorm(Y / 2).
type1 <- function(c) 1 - pnorm(sqrt(2) * qnorm(c))        # worst case mu = 0: Y ~ N(0, 2)
say("c = 0.5   -> type I error %.4f", type1(0.5))
say("c = 10/11 -> type I error %.4f", type1(10 / 11))
c_05 <- pnorm(qnorm(0.95) / sqrt(2))
say("c giving type I error 0.05: %.4f (check %.4f)", c_05, type1(c_05))
set.seed(1); y_null <- rnorm(1e6, 0, sqrt(2)); set.seed(2026)
say("simulation check at c = %.4f: %.4f", c_05, mean(pnorm(y_null / 2) > c_05))

save_fig("ex01_type1", width = 6, height = 3.8, function() {
  cc <- seq(0.5, 0.999, length.out = 500)
  plot(cc, type1(cc), type = "l", lwd = 2.5, col = PALETTE["navy"], xlab = "threshold c",
       ylab = "type I error (at mu = 0)")
  abline(h = 0.05, v = c_05, lty = 2, col = PALETTE["red"])
  points(c(0.5, 10 / 11), type1(c(0.5, 10 / 11)), pch = 19, col = PALETTE["gold"])
})


# ---- Exercise 2: should we run the confirmatory trial? --------------------------------

section("Ex 2  trial decision")

# Reconstructed set-up (see ch09.tex): power 0.8 if the drug is better, false confirmation
# probability alpha otherwise, payoff c*X if the trial confirms, cost X.
run_trial <- function(theta, c, alpha = 0.05) c * (0.8 * theta + alpha * (1 - theta)) > 1
theta_min <- function(c, alpha = 0.05) pmin(1, pmax(0, (1 / c - alpha) / (0.8 - alpha)))
for (cc in c(1.25, 1.5, 2, 5, 20)) say("c = %5.2f: run the trial if theta > %.3f", cc, theta_min(cc))
save_fig("ex02_decision", width = 6, height = 3.8, function() {
  cc <- seq(1, 20, length.out = 500)
  plot(cc, theta_min(cc), type = "l", lwd = 2.5, col = PALETTE["navy"], log = "x", ylim = c(0, 1),
       xlab = "payoff multiple c (log scale)", ylab = expression("run the trial if " * theta > "curve"))
  lines(cc, theta_min(cc, alpha = 0), lwd = 2, lty = 2, col = PALETTE["red"])
  legend("topright", c("false confirmation 5%", "no false confirmation"), lty = 1:2, lwd = 2,
         col = PALETTE[c("navy", "red")], bty = "n")
})


# ---- Exercise 3: bias-variance tradeoff for a proportion -------------------------------

section("Ex 3  sample proportion vs Jeffreys posterior mean")

mse1 <- function(th, n) th * (1 - th) / n
mse2 <- function(th, n) (n * th * (1 - th) + (0.5 - th)^2) / (n + 1)^2
for (n in c(5, 20, 100)) {
  lim <- uniroot(function(t) mse2(t, n) - mse1(t, n), c(1e-9, 0.5))$root
  say("n = %3d: Jeffreys mean has smaller MSE for theta in (%.3f, %.3f)", n, lim, 1 - lim)
}
say("limit as n -> infinity: theta(1-theta) > 1/12, i.e. theta in (%.3f, %.3f)", (1 - sqrt(2 / 3)) / 2, (1 + sqrt(2 / 3)) / 2)
save_fig("ex03_mse", width = 7.5, height = 3.2, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.4, 2, 0.6))
  th <- seq(0, 1, length.out = 401)
  for (n in c(5, 20)) {
    plot(th, mse1(th, n), type = "l", lwd = 2.5, col = PALETTE["navy"], ylim = c(0, max(mse1(th, n), mse2(th, n))),
         xlab = expression(theta[0]), ylab = "MSE", main = paste("n =", n))
    lines(th, mse2(th, n), lwd = 2.5, lty = 2, col = PALETTE["red"])
    legend("top", c("Y/n", "(Y+1/2)/(n+1)"), lwd = 2.5, lty = 1:2, col = PALETTE[c("navy", "red")], bty = "n", cex = 0.8)
  }
})


# ---- Exercise 4: posterior mean and mode of a normal variance --------------------------

section("Ex 4  estimators of sigma^2")

a4 <- 0.001; b4 <- 0.001; n4 <- 10; ss4 <- 200
A <- a4 + n4 / 2; B <- b4 + ss4 / 2
say("posterior InvGamma(%.3f, %.3f): mean %.3f   mode %.3f   (MLE ss/n = %.1f)", A, B, B / (A - 1), B / (A + 1), ss4 / n4)
dinvgamma <- function(x, a, b) exp(a * log(b) - lgamma(a) - (a + 1) * log(x) - b / x)
save_fig("ex04_posterior", width = 6, height = 3.8, function() {
  x <- seq(2, 80, length.out = 800)
  plot(x, dinvgamma(x, A, B), type = "l", lwd = 2.5, col = PALETTE["navy"], xlab = expression(sigma^2),
       ylab = "Posterior density")
  abline(v = c(B / (A - 1), B / (A + 1)), lty = c(1, 2), col = PALETTE[c("red", "gold")], lwd = 2)
  legend("topright", c("posterior mean", "posterior mode"), lty = 1:2, lwd = 2, col = PALETTE[c("red", "gold")], bty = "n")
})
# MSE as n grows (true sigma^2 = 20), exact formulas via SS ~ sigma^2 chi^2_n
mse_est <- function(n, s2 = 20, which = "mean") {
  m <- if (which == "mean") 1 / (a4 + n / 2 - 1) else 1 / (a4 + n / 2 + 1)
  est_mean <- m * (b4 + n * s2 / 2); est_var <- m^2 * (s2^2 * 2 * n) / 4
  est_var + (est_mean - s2)^2
}
for (n in c(10, 100, 1000)) say("n = %4d: MSE(mean) %.2f   MSE(mode) %.2f", n, mse_est(n), mse_est(n, which = "mode"))


# ---- Exercise 5: frequentist coverage of Jeffreys credible intervals -------------------

section("Ex 5  coverage of equal-tailed Jeffreys intervals")

coverage5 <- function(n, th) {
  y <- 0:n
  lo <- qbeta(0.025, y + 0.5, n - y + 0.5); hi <- qbeta(0.975, y + 0.5, n - y + 0.5)
  sum(dbinom(y, n, th) * (th >= lo & th <= hi))
}
th5 <- seq(0.05, 0.5, by = 0.05); n5 <- c(1, 5, 10, 25)
cov5 <- sapply(n5, function(n) sapply(th5, coverage5, n = n))
dimnames(cov5) <- list(theta = th5, n = n5)
print(round(cov5, 3))
sim5 <- sapply(n5, function(n) {                       # Monte Carlo check (10,000 datasets each)
  sapply(th5, function(th) {
    y <- rbinom(1e4, n, th)
    mean(th >= qbeta(0.025, y + 0.5, n - y + 0.5) & th <= qbeta(0.975, y + 0.5, n - y + 0.5))
  })
})
say("max |simulated - exact| coverage: %.4f", max(abs(sim5 - cov5)))
say("average coverage by n: %s", paste(sprintf("n=%d %.3f", n5, colMeans(cov5)), collapse = "  "))
save_fig("ex05_coverage", width = 6.5, height = 4, function() {
  matplot(th5, cov5, type = "b", pch = 19, lty = 1, lwd = 2, col = PALETTE[1:4], ylim = c(0.7, 1),
          xlab = expression("true " * theta[0]), ylab = "coverage of 95% interval")
  abline(h = 0.95, lty = 2, col = PALETTE["grey"])
  legend("bottomright", paste("n =", n5), col = PALETTE[1:4], lwd = 2, pch = 19, bty = "n")
})


# ---- Exercise 6: Bayesian re-analysis of a published logistic regression ---------------

section("Ex 6  re-analysing a published estimate")

b6 <- 1.8; se6 <- (3.17 - 0.42) / (2 * qnorm(0.975))
say("implied SE %.4f   z = %.3f   implied two-sided p = %.4f (reported 0.01)", se6, b6 / se6, 2 * pnorm(-b6 / se6))
for (tau in c(Inf, 2, 1, 0.5)) {
  if (is.infinite(tau)) { m <- b6; s <- se6 } else { v <- 1 / (1 / se6^2 + 1 / tau^2); m <- v * b6 / se6^2; s <- sqrt(v) }
  bf <- if (is.infinite(tau)) NA else dnorm(b6, 0, sqrt(se6^2 + tau^2)) / dnorm(b6, 0, se6)
  say("prior N(0, %s^2): posterior mean %.3f  95%% CI %s   P(beta > 0) = %.4f   OR mean %.2f   BF(H1:H0) %s",
      format(tau), m, fmt_ci(m + c(-1, 1) * qnorm(0.975) * s, 2), pnorm(m / s), exp(m + s^2 / 2),
      if (is.na(bf)) "-" else sprintf("%.2f", bf))
}
save_fig("ex06_posterior", width = 6.5, height = 3.8, function() {
  x <- seq(-1, 4.5, length.out = 600)
  plot(x, dnorm(x, b6, se6), type = "l", lwd = 2.5, col = PALETTE["navy"], xlab = expression(beta[1]),
       ylab = "Posterior density", ylim = c(0, 0.85))
  for (k in 1:2) {
    tau <- c(1, 0.5)[k]; v <- 1 / (1 / se6^2 + 1 / tau^2)
    lines(x, dnorm(x, v * b6 / se6^2, sqrt(v)), lwd = 2, lty = k + 1, col = PALETTE[k + 1])
  }
  abline(v = 0, lty = 3, col = PALETTE["grey"])
  legend("topright", c("flat prior", "N(0, 1) prior", "N(0, 0.5^2) prior"), lwd = 2, lty = 1:3,
         col = PALETTE[1:3], bty = "n")
})


# ---- Exercise 7: sample size for a target posterior SD ---------------------------------

section("Ex 7  sample size")

post_sd_expected <- function(n, th = 0.4) {           # average over the binomial distribution of Y
  y <- 0:n; a <- y + 1; b <- n - y + 1
  sum(dbinom(y, n, th) * sqrt(a * b / ((a + b)^2 * (a + b + 1))))
}
n_approx <- ceiling(0.4 * 0.6 / 0.01^2 - 3)
say("approximation theta(1-theta)/(n+3) = 0.01^2 gives n = %d", n_approx)
say("expected posterior sd at n = %d: %.5f;  smallest n with expected sd <= 0.01: %d",
    n_approx, post_sd_expected(n_approx), min(which(sapply(2300:2450, post_sd_expected) <= 0.01)) + 2299)
for (th in c(0.3, 0.5)) say("if theta were %.1f instead: posterior sd at n = %d is %.5f", th, n_approx, post_sd_expected(n_approx, th))
