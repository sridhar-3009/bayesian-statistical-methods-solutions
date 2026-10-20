# =============================================================================
#  Chapter 2 -- From prior information to posterior inference
#  Reproduces every number and figure in ch02/ch02.tex.
#
#  Run from the repository root:   Rscript ch02/code/ch02.R
#  Output:  figures -> ch02/figures/*.pdf     numbers -> console
#  Needs base R only (no extra packages).
# =============================================================================

source(file.path("common", "helpers.R"))
set.seed(2026)
FIG_DIR <- file.path("ch02", "figures")


# ---- Exercise 1: frequentist coverage of the normal-normal interval ---------

section("Ex 1  coverage of the posterior interval")

# Coverage of  n*Ybar/(n+m) +/- 1.96 sigma/sqrt(n+m)  when the true mean is mu.
coverage <- function(mu, m, n = 10, sigma = 1) {
  half  <- 1.96 * sqrt((n + m) / n)
  shift <- m * mu / (sqrt(n) * sigma)
  pnorm(half + shift) - pnorm(-half + shift)
}

mu_grid <- seq(-3, 3, length.out = 601)
m_vals  <- c(0, 1, 5, 20)

save_fig("ex01_coverage", function() {
  plot(NA, xlim = range(mu_grid), ylim = c(0, 1), xlab = expression("true " * mu),
       ylab = "Frequentist coverage")
  for (k in seq_along(m_vals))
    lines(mu_grid, coverage(mu_grid, m_vals[k]), lwd = 2.5, col = PALETTE[k])
  abline(h = 0.95, lty = 2, col = PALETTE["grey"])
  legend("bottom", legend = paste("m =", m_vals), col = PALETTE[seq_along(m_vals)],
         lwd = 2.5, bty = "n", horiz = TRUE)
})

for (m in m_vals) say("m = %2d   coverage at mu = 0: %.4f   at mu = 2: %.4f",
                      m, coverage(0, m), coverage(2, m))

# Monte Carlo check for m = 0
reps <- 1e5; n <- 10; mu_true <- 1.7
ybar <- rnorm(reps, mu_true, 1 / sqrt(n))
say("simulated coverage, m = 0, mu = %.1f: %.4f", mu_true,
    mean(abs(ybar - mu_true) <= 1.96 / sqrt(n)))


# ---- Exercise 2: Mr. October -------------------------------------------------

section("Ex 2  Reggie Jackson home-run rates")

prior2 <- c(a = 0.01, b = 0.01)
seasons <- data.frame(row.names = c("regular", "world_series"),
                      hr = c(563, 10), games = c(2820, 27))
seasons$shape <- prior2["a"] + seasons$hr
seasons$rate  <- prior2["b"] + seasons$games

for (s in rownames(seasons)) with(seasons[s, ],
  say("%-13s Gamma(%.2f, %.2f)  mean = %.4f  sd = %.4f  95%% CI = %s",
      s, shape, rate, shape / rate, sqrt(shape) / rate, fmt_ci(ci95(qgamma, shape, rate), 4)))

lam_reg <- rgamma(1e6, seasons["regular", "shape"], seasons["regular", "rate"])
lam_ws  <- rgamma(1e6, seasons["world_series", "shape"], seasons["world_series", "rate"])
ratio   <- lam_ws / lam_reg
say("P(lambda_WS > lambda_reg) = %.4f   ratio median = %.3f  95%% CI = %s",
    mean(lam_ws > lam_reg), median(ratio), fmt_ci(quantile(ratio, c(0.025, 0.975)), 3))

save_fig("ex02_reggie", function() {
  lam <- seq(0, 0.8, length.out = 1601)
  d_reg <- dgamma(lam, seasons["regular", "shape"], seasons["regular", "rate"])
  d_ws  <- dgamma(lam, seasons["world_series", "shape"], seasons["world_series", "rate"])
  plot(lam, d_reg / max(d_reg), type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression("home runs per game " * lambda), ylab = "Posterior (scaled to max 1)")
  lines(lam, d_ws / max(d_ws), lwd = 2.5, col = PALETTE["red"])
  legend("topright", c("regular season", "World Series"), lwd = 2.5,
         col = PALETTE[c("navy", "red")], bty = "n")
})


# ---- Exercise 4: factory injuries --------------------------------------------

section("Ex 4  injury rate")

injuries <- c(4, 5, 2, 7, 6, 6)
injury_priors <- list("Gamma(0.1, 0.1)" = c(0.1, 0.1), "Gamma(0.01, 0.01)" = c(0.01, 0.01),
                      "Gamma(1, 1)" = c(1, 1))
for (nm in names(injury_priors)) {
  a <- injury_priors[[nm]][1] + sum(injuries)
  b <- injury_priors[[nm]][2] + length(injuries)
  say("%-18s -> Gamma(%.2f, %.2f)  mean = %.3f  95%% CI = %s  P(lambda > 4) = %.4f",
      nm, a, b, a / b, fmt_ci(ci95(qgamma, a, b), 2), pgamma(4, a, b, lower.tail = FALSE))
}

save_fig("ex04_injuries", function() {
  lam <- seq(1.5, 9, length.out = 801)
  a <- 0.1 + sum(injuries); b <- 0.1 + length(injuries)
  dens <- dgamma(lam, a, b)
  plot(lam, dens, type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression("injuries per month " * lambda), ylab = "Posterior density")
  over <- lam > 4
  polygon(c(4, lam[over], max(lam)), c(0, dens[over], 0),
          col = adjustcolor(PALETTE["red"], 0.3), border = NA)
  abline(v = 4, lty = 2, col = PALETTE["red"])
})


# ---- Exercise 5: negative binomial / beta ------------------------------------

section("Ex 5  negative binomial likelihood, beta prior")

a5 <- 1 + 5; b5 <- 1 + 10                    # Beta(a + m, b + Y)
say("posterior Beta(%d, %d)  mean = %.4f  95%% CI = %s", a5, b5, a5 / (a5 + b5),
    fmt_ci(ci95(qbeta, a5, b5), 4))

save_fig("ex05_negbin", function() {
  theta <- seq(0, 1, length.out = 801)
  ci <- ci95(qbeta, a5, b5)
  dens <- dbeta(theta, a5, b5)
  plot(theta, dens, type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(theta), ylab = "Posterior density")
  inside <- theta >= ci[1] & theta <= ci[2]
  polygon(c(ci[1], theta[inside], ci[2]), c(0, dens[inside], 0),
          col = adjustcolor(PALETTE["navy"], 0.2), border = NA)
})


# ---- Exercise 6: regression to the mean --------------------------------------

section("Ex 6  heart rate and regression to the mean")

prior_mean <- 64; prior_sd <- 10; error_sd <- 5
shrink <- prior_sd^2 / (prior_sd^2 + error_sd^2)          # weight on the data
post_sd <- sqrt(1 / (1 / prior_sd^2 + 1 / error_sd^2))
say("posterior mean = %.1f + %.1f y   posterior sd = %.3f",
    (1 - shrink) * prior_mean, shrink, post_sd)

save_fig("ex06_shrinkage", width = 5.5, height = 4.5, function() {
  y <- seq(30, 110, length.out = 200)
  plot(y, (1 - shrink) * prior_mean + shrink * y, type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = "measured heart rate y (bpm)", ylab = "E(true heart rate | y)", asp = 1)
  abline(0, 1, lty = 2, col = PALETTE["grey"])
  abline(h = prior_mean, v = prior_mean, lty = 3, col = PALETTE["grey"])
  legend("topleft", c("posterior mean", "no shrinkage (y = x)"), lty = 1:2, lwd = c(2.5, 1),
         col = PALETTE[c("navy", "grey")], bty = "n")
})

# (b) exact value: E[Y1 - Y2 | Y1 > 90] = (1 - shrink) * (E[Y1 | Y1 > 90] - 64)
marg_sd <- sqrt(prior_sd^2 + error_sd^2)
z90 <- (90 - prior_mean) / marg_sd
e_y1 <- prior_mean + marg_sd * dnorm(z90) / pnorm(z90, lower.tail = FALSE)
say("E(Y1 | Y1 > 90) = %.3f   expected spurious 'effect' = %.3f bpm", e_y1, (1 - shrink) * (e_y1 - prior_mean))

people <- 2e6
true_hr <- rnorm(people, prior_mean, prior_sd)
before  <- true_hr + rnorm(people, 0, error_sd)
after   <- true_hr + rnorm(people, 0, error_sd)        # treatment does nothing
enrolled <- before > 90
say("simulation: %d enrolled, mean before - after = %.3f bpm",
    sum(enrolled), mean(before[enrolled] - after[enrolled]))


# ---- Exercise 7: eliciting a wildfire prior ----------------------------------

section("Ex 7  wildfire prior")

# Variance a / b^2 = 100  =>  b = sqrt(a / 100); then solve median = 75 for a.
median_gap <- function(a) qgamma(0.5, a, sqrt(a / 100)) - 75
a7 <- uniroot(median_gap, c(10, 200))$root
b7 <- sqrt(a7 / 100)
say("a = %.3f  b = %.4f   mean = %.3f  sd = %.3f  P(lambda > 75) = %.4f",
    a7, b7, a7 / b7, sqrt(a7) / b7, pgamma(75, a7, b7, lower.tail = FALSE))
say("moment-only choice a = 56.25, b = 0.75:  P(lambda > 75) = %.4f",
    pgamma(75, 56.25, 0.75, lower.tail = FALSE))


# ---- Exercise 8: inverse-gamma variance ---------------------------------------

section("Ex 8  measurement-error variance")

n8 <- 20; ss8 <- 15
for (c8 in c(1, 2)) {
  p <- sapply(c(0.1, 1), function(ab)                 # P(sigma > c) = P(1/sigma^2 < 1/c^2)
    pgamma(1 / c8^2, ab + n8 / 2, ab + ss8 / 2))
  say("c = %d   a=b=0.1: %.3e   a=b=1: %.3e   ratio = %.3f", c8, p[1], p[2], p[1] / p[2])
}


# ---- Exercise 9: induced prior for logistic regression -----------------------

section("Ex 9  induced prior on p(X)")

induced_p <- function(c, S) plogis(rnorm(S, 0, c) + rnorm(S) * rnorm(S, 0, c))
ks_to_uniform <- function(p) {
  p <- sort(p); S <- length(p)
  max(abs(p - seq_len(S) / S))
}

c_grid <- seq(0.5, 3, by = 0.05)
ks <- sapply(c_grid, function(c) ks_to_uniform(induced_p(c, 2e5)))
c_best <- c_grid[which.min(ks)]
say("best c on grid = %.2f  (KS distance %.4f)", c_best, min(ks))

p_final <- induced_p(c_best, 1e6)
say("S = 1e6: KS distance %.4f", ks_to_uniform(p_final))

save_fig("ex09_induced", width = 7.5, height = 3.2, function() {
  par(mfrow = c(1, 3), mar = c(4, 3, 2.2, 0.5))
  for (c in c(0.5, c_best, 10)) {
    hist(if (c == c_best) p_final else induced_p(c, 1e6), breaks = 40, freq = FALSE,
         col = PALETTE["navy"], border = "white", lwd = 0.3, xlab = "p(X)", ylab = "",
         main = sprintf("c = %.2f", c))
    abline(h = 1, lty = 2, col = PALETTE["red"])
  }
})


# ---- Exercise 10: adverse side effects ----------------------------------------

section("Ex 10  adverse side effects")

ase <- rep(c(0, 1, 2, 10), times = c(30, 12, 6, 2))
ase_priors <- list("Gamma(0.01, 0.01)" = c(0.01, 0.01), "Gamma(0.1, 0.1)" = c(0.1, 0.1),
                   "Gamma(1, 1)" = c(1, 1))
for (nm in names(ase_priors)) {
  a <- ase_priors[[nm]][1] + sum(ase); b <- ase_priors[[nm]][2] + length(ase)
  say("%-18s -> Gamma(%.2f, %.2f)  mean = %.4f  95%% CI = %s  P(lambda > 1) = %.4f",
      nm, a, b, a / b, fmt_ci(ci95(qgamma, a, b), 3), pgamma(1, a, b, lower.tail = FALSE))
}

a10 <- 0.01 + sum(ase); b10 <- 0.01 + length(ase)
lambda_hat <- a10 / b10

save_fig("ex10_posterior", width = 6, height = 3.8, function() {
  lam <- seq(0.5, 1.35, length.out = 600)
  plot(lam, dgamma(lam, a10, b10), type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(lambda), ylab = "Posterior density")
  abline(v = 1, lty = 2, col = PALETTE["red"])
})

save_fig("ex10_fit", width = 6.5, height = 3.8, function() {
  k <- 0:10
  observed <- table(factor(ase, levels = k)) / length(ase)
  centers <- barplot(rbind(observed, dpois(k, lambda_hat)), beside = TRUE, names.arg = k,
                     col = PALETTE[c("navy", "gold")], border = NA,
                     xlab = "number of side effects", ylab = "Proportion")
  legend("topright", c("observed", sprintf("Poisson(%.2f)", lambda_hat)),
         fill = PALETTE[c("navy", "gold")], border = NA, bty = "n")
})
say("expected zeros under Poisson(%.3f): %.1f of 50 (observed 30);  P(Y >= 10) = %.1e",
    lambda_hat, 50 * dpois(0, lambda_hat), ppois(9, lambda_hat, lower.tail = FALSE))
say("sample mean = %.2f   sample variance = %.2f", mean(ase), var(ase))


# ---- Exercise 11: smoking cessation (geometric) -------------------------------

section("Ex 11  smoking cessation")

say("prior P(Y >= 7) = 1/8 = %.4f", 1 / 8)
say("posterior Beta(2, 6)   mean = %.4f   95%% CI = %s", 2 / 8, fmt_ci(ci95(qbeta, 2, 6), 3))
say("posterior predictive P(Y* >= 7 | Y = 5) = B(2, 13) / B(2, 6) = %.4f", beta(2, 13) / beta(2, 6))

theta_draws <- rbeta(1e6, 2, 6)
y_star <- rgeom(1e6, theta_draws)
say("Monte Carlo check: %.4f", mean(y_star >= 7))

save_fig("ex11_posterior", width = 6, height = 3.8, function() {
  theta <- seq(0, 1, length.out = 801)
  plot(theta, dbeta(theta, 2, 6), type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression("daily relapse probability " * theta), ylab = "Density")
  lines(theta, dunif(theta), lty = 2, col = PALETTE["grey"])
  legend("topright", c("posterior Beta(2, 6)", "prior Uniform(0, 1)"), lty = 1:2,
         lwd = c(2.5, 1), col = PALETTE[c("navy", "grey")], bty = "n")
})


# ---- Exercise 12: censored exponential survival -------------------------------

section("Ex 12  censored survival times")

relapse  <- c(0.5, 0.9, 1.1, 2.4, 3.8, 4.0, 5.2)
censored <- rep(6, 20 - length(relapse))
exposure <- sum(relapse) + sum(censored)
say("events d = %d   total exposure = %.1f weeks", length(relapse), exposure)

surv_priors <- list("Gamma(0.01, 0.01)" = c(0.01, 0.01), "Gamma(0.1, 0.1)" = c(0.1, 0.1),
                    "Gamma(1, 1)" = c(1, 1), "Gamma(2, 10)" = c(2, 10))
for (nm in names(surv_priors)) {
  a <- surv_priors[[nm]][1] + length(relapse); b <- surv_priors[[nm]][2] + exposure
  say("%-18s -> Gamma(%.2f, %.2f)  mean = %.4f  95%% CI = %s", nm, a, b, a / b,
      fmt_ci(ci95(qgamma, a, b), 4))
}
say("naive (ignore censored): mean of 7 times -> rate %.3f", 1 / mean(relapse))


# ---- Exercise 13: mixture-of-experts prior ------------------------------------

section("Ex 13  alarm sensitivity, two experts")

beta_from_moments <- function(mean, sd) {
  total <- mean * (1 - mean) / sd^2 - 1
  c(a = mean * total, b = (1 - mean) * total)
}
expert1 <- beta_from_moments(0.95, 0.05)
expert2 <- beta_from_moments(0.80, 0.20)
say("expert 1: Beta(%.2f, %.2f)   expert 2: Beta(%.2f, %.2f)",
    expert1["a"], expert1["b"], expert2["a"], expert2["b"])

y13 <- 5; n13 <- 5
marginal <- function(e) beta(e["a"] + y13, e["b"] + n13 - y13) / beta(e["a"], e["b"])
weights <- 0.5 * c(marginal(expert1), marginal(expert2))
weights <- weights / sum(weights)
say("posterior mixture weights: expert 1 = %.4f, expert 2 = %.4f", weights[1], weights[2])

theta13 <- seq(0.3, 1, length.out = 1401)
dens_prior <- function(e) dbeta(theta13, e["a"], e["b"])
dens_post  <- function(e) dbeta(theta13, e["a"] + y13, e["b"] + n13 - y13)

post_means <- c(uniform = 6 / 7,
                expert1 = unname((expert1["a"] + y13) / (sum(expert1) + n13)),
                expert2 = unname((expert2["a"] + y13) / (sum(expert2) + n13)))
post_means["mixture"] <- sum(weights * post_means[c("expert1", "expert2")])
print(round(post_means, 4))

save_fig("ex13_priors", function() {
  mix <- 0.5 * dens_prior(expert1) + 0.5 * dens_prior(expert2)
  plot(theta13, dens_prior(expert1), type = "l", lwd = 2.5, col = PALETTE["navy"], ylim = c(0, 12),
       xlab = expression("alarm probability " * theta), ylab = "Prior density")
  lines(theta13, dens_prior(expert2), lwd = 2.5, col = PALETTE["gold"])
  lines(theta13, mix, lwd = 2.5, lty = 2, col = PALETTE["red"])
  legend("topleft", c("expert 1: Beta(17.1, 0.9)", "expert 2: Beta(2.4, 0.6)", "50/50 mixture"),
         col = PALETTE[c("navy", "gold", "red")], lty = c(1, 1, 2), lwd = 2.5, bty = "n")
})

save_fig("ex13_posteriors", function() {
  mix <- weights[1] * dens_post(expert1) + weights[2] * dens_post(expert2)
  plot(theta13, dens_post(expert1), type = "l", lwd = 2.5, col = PALETTE["navy"], ylim = c(0, 25),
       xlab = expression("alarm probability " * theta), ylab = "Posterior density")
  lines(theta13, dens_post(expert2), lwd = 2.5, col = PALETTE["gold"])
  lines(theta13, dbeta(theta13, 6, 1), lwd = 2.5, col = PALETTE["green"])
  lines(theta13, mix, lwd = 2.5, lty = 2, col = PALETTE["red"])
  legend("topleft", c("uniform prior", "expert 1", "expert 2", "mixture"),
         col = PALETTE[c("green", "navy", "gold", "red")], lty = c(1, 1, 1, 2), lwd = 2.5, bty = "n")
})


# ---- Exercise 14: priors under reparameterisation -----------------------------

section("Ex 14  probability vs odds")

S <- 1e6
odds <- function(theta) theta / (1 - theta)
keep <- function(g) g[g <= 100]

save_fig("ex14_reparam", width = 7.5, height = 5.2, function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 2.2, 0.8))

  g <- keep(odds(runif(S)))
  hist(g[g < 10], breaks = 60, freq = FALSE, col = PALETTE["navy"], border = NA,
       main = expression("(a) " * theta %~% U(0, 1) %=>% gamma), xlab = expression(gamma ~ "(shown < 10)"))
  curve(1 / (1 + x)^2 / punif(10 / 11), add = TRUE, lwd = 2, col = PALETTE["red"])

  g <- keep(odds(rbeta(S, 0.5, 0.5)))
  hist(g[g < 10], breaks = 60, freq = FALSE, col = PALETTE["navy"], border = NA,
       main = expression("(b) " * theta %~% Beta(0.5, 0.5) %=>% gamma), xlab = expression(gamma ~ "(shown < 10)"))
  curve(1 / (pi * sqrt(x) * (1 + x)) / pbeta(10 / 11, 0.5, 0.5), add = TRUE, lwd = 2, col = PALETTE["red"])

  t <- runif(S, 0, 100); t <- t / (1 + t)
  hist(t, breaks = 60, freq = FALSE, col = PALETTE["navy"], border = NA,
       main = expression("(c) " * gamma %~% U(0, 100) %=>% theta), xlab = expression(theta))

  t <- rexp(S); t <- t / (1 + t)
  hist(t, breaks = 60, freq = FALSE, col = PALETTE["navy"], border = NA,
       main = expression("(d) " * gamma %~% Gamma(1, 1) %=>% theta), xlab = expression(theta))
  curve(exp(-x / (1 - x)) / (1 - x)^2, add = TRUE, lwd = 2, col = PALETTE["red"])
})
say("(c) P(theta > 0.9) = %.3f   (d) P(theta > 0.9) = %.4f",
    mean(runif(S, 0, 100) > 9), pexp(9, lower.tail = FALSE))


# ---- Exercise 18: signal-to-noise mixture prior -------------------------------

section("Ex 18  clinical-trial mixture prior")

mix_w  <- c(0.32, 0.31, 0.30, 0.07)
mix_sd <- c(1.17, 1.74, 2.38, 5.73)

# Posterior of delta | z for a zero-mean normal mixture prior and z | delta ~ N(delta, 1).
mixture_posterior <- function(z, w = mix_w, tau = mix_sd) {
  v    <- tau^2 / (1 + tau^2)
  post <- w * dnorm(z, 0, sqrt(1 + tau^2))
  list(w = post / sum(post), mean = v * z, sd = sqrt(v))
}
p_positive <- function(z, ...) {
  pp <- mixture_posterior(z, ...)
  sum(pp$w * pnorm(pp$mean / pp$sd))
}

pp2 <- mixture_posterior(2)
say("z = 2: weights %s", paste(sprintf("%.3f", pp2$w), collapse = " "))
say("       component means %s", paste(sprintf("%.3f", pp2$mean), collapse = " "))
say("       posterior mean = %.3f   P(delta > 0) = %.4f", sum(pp2$w * pp2$mean), p_positive(2))
say("flat-ish N(0, 100^2) prior: posterior mean = %.4f   P(delta > 0) = %.4f",
    2 * 1e4 / (1 + 1e4), p_positive(2, w = 1, tau = 100))

delta <- seq(-3, 7, length.out = 1001)
mix_density <- function(pp) rowSums(sapply(seq_along(pp$w),
  function(k) pp$w[k] * dnorm(delta, pp$mean[k], pp$sd[k])))

save_fig("ex18_posterior", function() {
  prior <- rowSums(sapply(1:4, function(k) mix_w[k] * dnorm(delta, 0, mix_sd[k])))
  plot(delta, mix_density(pp2), type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(delta), ylab = "Density", ylim = c(0, 0.45))
  lines(delta, prior, lwd = 2, lty = 3, col = PALETTE["grey"])
  lines(delta, dnorm(delta, 2 * 1e4 / (1 + 1e4), sqrt(1e4 / (1 + 1e4))), lwd = 2.5, lty = 2,
        col = PALETTE["red"])
  abline(v = 2, lty = 3, col = PALETTE["grey"])
  legend("topright", c("posterior, mixture prior", "posterior, N(0, 100^2) prior", "mixture prior"),
         col = PALETTE[c("navy", "red", "grey")], lty = c(1, 2, 3), lwd = 2.5, bty = "n")
})

z_grid <- seq(0, 5, by = 0.1)
p_mix  <- sapply(z_grid, p_positive)
p_flat <- sapply(z_grid, p_positive, w = 1, tau = 100)
save_fig("ex18_pz", function() {
  plot(z_grid, p_mix, type = "l", lwd = 2.5, col = PALETTE["navy"], ylim = c(0.5, 1),
       xlab = "observed z", ylab = expression(Prob(delta > 0 ~ "|" ~ z)))
  lines(z_grid, p_flat, lwd = 2.5, lty = 2, col = PALETTE["red"])
  legend("bottomright", c("mixture prior", "N(0, 100^2) prior"), lwd = 2.5, lty = 1:2,
         col = PALETTE[c("navy", "red")], bty = "n")
})
print(data.frame(z = c(1, 1.96, 2, 3, 4), mixture = round(sapply(c(1, 1.96, 2, 3, 4), p_positive), 4),
                 flat = round(sapply(c(1, 1.96, 2, 3, 4), p_positive, w = 1, tau = 100), 4)),
      row.names = FALSE)


# ---- Exercises 19 and 20: Jeffreys priors -------------------------------------

section("Ex 19-20  Jeffreys priors")

save_fig("ex19_20_jeffreys", width = 7.5, height = 3.4, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.4, 2.2, 0.8))
  lam <- seq(0.02, 5, length.out = 500)
  plot(lam, 1 / sqrt(lam), type = "l", lwd = 2.5, col = PALETTE["navy"], main = "Poisson",
       xlab = expression(lambda), ylab = expression(pi(lambda) %prop% lambda^{-1/2}))
  plot(lam, 1 / lam, type = "l", lwd = 2.5, col = PALETTE["navy"], main = "Exponential",
       xlab = expression(lambda), ylab = expression(pi(lambda) %prop% lambda^{-1}), ylim = c(0, 10))
})


# ---- Exercise 21: default analysis of Y = 10 ----------------------------------

section("Ex 21  Poisson with Y = 10, default priors")

default_priors <- list("Jeffreys  Gamma(0.5, 0)" = c(0.5, 0), "flat      Gamma(1, 0)" = c(1, 0),
                       "vague     Gamma(0.01, 0.01)" = c(0.01, 0.01))
for (nm in names(default_priors)) {
  a <- default_priors[[nm]][1] + 10; b <- default_priors[[nm]][2] + 1
  say("%-28s -> Gamma(%.2f, %.2f)  mean = %.3f  sd = %.3f  95%% CI = %s",
      nm, a, b, a / b, sqrt(a) / b, fmt_ci(ci95(qgamma, a, b), 2))
}

save_fig("ex21_posterior", width = 6, height = 3.8, function() {
  lam <- seq(0, 25, length.out = 800)
  plot(lam, dgamma(lam, 10.5, 1), type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(lambda), ylab = "Posterior density")
  ci <- ci95(qgamma, 10.5, 1)
  abline(v = ci, lty = 2, col = PALETTE["red"])
})


# ---- Exercise 22: HOV lane survey, Jeffreys vs empirical Bayes ----------------

section("Ex 22  HOV lane approval")

hov <- data.frame(county = 1:10,
                  approve    = c(12, 90, 80, 5, 63, 15, 67, 22, 56, 33),
                  disapprove = c(50, 150, 63, 10, 63, 8, 56, 19, 63, 19))
hov$n <- hov$approve + hov$disapprove
hov$p_hat <- hov$approve / hov$n

m22 <- mean(hov$p_hat); v22 <- var(hov$p_hat)
eb  <- beta_from_moments(m22, sqrt(v22))
say("mean(p_hat) = %.4f   var(p_hat) = %.5f   ->  Beta(a = %.3f, b = %.3f), a + b = %.2f",
    m22, v22, eb["a"], eb["b"], sum(eb))

interval <- function(a0, b0) {
  A <- a0 + hov$approve; B <- b0 + hov$disapprove
  cbind(mean = A / (A + B), lo = qbeta(0.025, A, B), hi = qbeta(0.975, A, B))
}
jef <- interval(0.5, 0.5)
emp <- interval(eb["a"], eb["b"])
hov_table <- data.frame(county = hov$county, n = hov$n, p_hat = round(hov$p_hat, 3),
                        jef_mean = round(jef[, "mean"], 3), jef_lo = round(jef[, "lo"], 3),
                        jef_hi = round(jef[, "hi"], 3), eb_mean = round(emp[, "mean"], 3),
                        eb_lo = round(emp[, "lo"], 3), eb_hi = round(emp[, "hi"], 3))
print(hov_table, row.names = FALSE)
say("average width: Jeffreys %.3f   empirical Bayes %.3f",
    mean(jef[, "hi"] - jef[, "lo"]), mean(emp[, "hi"] - emp[, "lo"]))

save_fig("ex22_intervals", width = 6.5, height = 4.4, function() {
  y <- hov$county
  plot(NA, xlim = c(0, 1), ylim = c(10.6, 0.4), yaxt = "n", xlab = "approval probability",
       ylab = "county")
  axis(2, at = y, labels = y)
  abline(v = m22, lty = 3, col = PALETTE["grey"])
  segments(jef[, "lo"], y - 0.15, jef[, "hi"], y - 0.15, lwd = 3, col = PALETTE["navy"])
  segments(emp[, "lo"], y + 0.15, emp[, "hi"], y + 0.15, lwd = 3, col = PALETTE["red"])
  points(hov$p_hat, y, pch = 4, cex = 1.1)
  legend("topright", c("Jeffreys 95%", "empirical Bayes 95%", "sample proportion"),
         col = c(PALETTE[c("navy", "red")], "black"), lwd = c(3, 3, NA), pch = c(NA, NA, 4),
         bty = "n")
})
