# =============================================================================
#  Chapter 3 -- Computational approaches
#  Reproduces every number and figure in ch03/ch03.tex with hand-written
#  samplers (Gibbs, Metropolis, grids).  The JAGS versions of the same models
#  live in ch03/code/ch03_jags.R and are run by GitHub Actions, not locally.
#
#  Run from the repository root:   Rscript ch03/code/ch03.R
#  Needs base R + MASS (for the galaxies data).
# =============================================================================

source(file.path("common", "helpers.R"))
set.seed(2026)
FIG_DIR <- file.path("ch03", "figures")


# ---- Exercise 2: MAP vs posterior mean, known variances ------------------------

section("Ex 2  MAP with known variances")

y2 <- c(12, 10, 22); s2 <- c(3, 3, 10)
w2 <- 1 / s2^2
map2 <- sum(w2 * y2) / sum(w2)
log_post2 <- function(mu) sapply(mu, function(m) sum(dnorm(y2, m, s2, log = TRUE)))

mu_grid <- seq(0, 25, length.out = 5001)
dens2   <- exp(log_post2(mu_grid) - max(log_post2(mu_grid)))
dens2   <- dens2 / sum(dens2 * diff(mu_grid)[1])
mean2   <- sum(mu_grid * dens2) * diff(mu_grid)[1]
say("MAP = %.4f   numerical posterior mean = %.4f   posterior sd = %.4f",
    map2, mean2, 1 / sqrt(sum(w2)))

save_fig("ex02_map", width = 6, height = 3.8, function() {
  plot(mu_grid, dens2, type = "l", lwd = 2.5, col = PALETTE["navy"], xlim = c(5, 18),
       xlab = expression(mu), ylab = "Posterior density")
  abline(v = map2, lwd = 2, col = PALETTE["red"])
  abline(v = mean2, lwd = 2, lty = 2, col = PALETTE["gold"])
  legend("topright", c("posterior", "MAP", "posterior mean"), bty = "n",
         col = PALETTE[c("navy", "red", "gold")], lty = c(1, 1, 2), lwd = 2)
})


# ---- Exercise 3: Poisson with exposures, grid vs Bayesian CLT -----------------

section("Ex 3  Poisson exposure model")

N3 <- c(50, 100); y3 <- c(12, 25)
map3 <- min(sum(y3) / sum(N3), 20)
clt_sd3 <- sqrt(map3 / sum(N3))            # 1 / sqrt(observed information)
lam3 <- seq(0.1, 0.42, length.out = 2001)
lp3  <- sum(y3) * log(lam3) - sum(N3) * lam3
dens3 <- exp(lp3 - max(lp3)); dens3 <- dens3 / sum(dens3 * diff(lam3)[1])
say("MAP = %.4f  (grid argmax %.4f)   CLT sd = %.4f   exact Gamma(38, 150) sd = %.4f",
    map3, lam3[which.max(dens3)], clt_sd3, sqrt(38) / 150)

save_fig("ex03_clt", width = 6, height = 3.8, function() {
  plot(lam3, dens3, type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(lambda), ylab = "Posterior density")
  lines(lam3, dnorm(lam3, map3, clt_sd3), lwd = 2, lty = 2, col = PALETTE["red"])
  abline(v = map3, lty = 3, col = PALETTE["grey"])
  legend("topright", c("exact (grid)", "Bayesian CLT"), lty = 1:2, lwd = 2,
         col = PALETTE[c("navy", "red")], bty = "n")
})


# ---- Exercise 4: Gibbs sampler for a hierarchical variance model --------------

section("Ex 4  Gibbs sampler, InvGamma(a, b) variances")

gibbs4 <- function(y, a, iters = 20000) {
  n <- length(y)
  sigma2 <- y^2 + 1; b <- 1                          # initial values
  draws <- matrix(NA, iters, n + 1, dimnames = list(NULL, c(paste0("s2_", 1:n), "b")))
  for (it in seq_len(iters)) {
    sigma2 <- 1 / rgamma(n, a + 0.5, b + y^2 / 2)    # sigma_i^2 | b, y
    b      <- rgamma(1, n * a + 1, 1 + sum(1 / sigma2))   # b | sigma^2
    draws[it, ] <- c(sigma2, b)
  }
  draws
}

y4 <- 1:10
fits4 <- list(a10 = gibbs4(y4, a = 10), a1 = gibbs4(y4, a = 1))
burn <- 1:2000
for (nm in names(fits4)) {
  d <- fits4[[nm]][-burn, ]
  say("%s: b mean = %.3f (ESS %.0f)   sigma2_1 mean = %.3f (ESS %.0f)   sigma2_10 mean = %.2f (ESS %.0f)",
      nm, mean(d[, "b"]), ess(d[, "b"]), mean(d[, 1]), ess(d[, 1]), mean(d[, 10]), ess(d[, 10]))
}

save_fig("ex04_marginals", width = 7.5, height = 5.4, function() {
  par(mfrow = c(3, 4), mar = c(3.4, 2.6, 1.8, 0.4), mgp = c(1.8, 0.5, 0))
  d <- fits4$a10[-burn, ]
  for (j in 1:11) {
    dj <- density(d[, j])
    plot(dj, main = colnames(d)[j], xlab = "", ylab = "", lwd = 2, col = PALETTE["navy"],
         zero.line = FALSE)
  }
  plot.new()
})

save_fig("ex04_traces", width = 7.5, height = 3.4, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.2, 2, 0.6))
  plot(fits4$a10[1:5000, "b"], type = "l", col = PALETTE["navy"], lwd = 0.5,
       xlab = "iteration", ylab = "b", main = "a = 10")
  plot(fits4$a1[1:5000, "b"], type = "l", col = PALETTE["red"], lwd = 0.5,
       xlab = "iteration", ylab = "b", main = "a = 1")
})


# ---- Exercise 5: Gibbs sampler for a two-group normal model -------------------

section("Ex 5  two-group normal model")

n5 <- 50; m5 <- 50
y5 <- c(rnorm(n5, 10, 2), rnorm(m5, 10 + 1, 2))
grp <- rep(0:1, c(n5, m5))
write.csv(data.frame(y = y5, group = grp), file.path("ch03", "code", "ex05_data.csv"), row.names = FALSE)

gibbs5 <- function(y, grp, iters = 20000) {
  N <- length(y); m <- sum(grp)
  mu <- mean(y); delta <- 0; sigma2 <- var(y)
  draws <- matrix(NA, iters, 3, dimnames = list(NULL, c("mu", "delta", "sigma")))
  for (it in seq_len(iters)) {
    prec  <- N / sigma2 + 1 / 100^2
    mu    <- rnorm(1, sum(y - delta * grp) / sigma2 / prec, 1 / sqrt(prec))
    prec  <- m / sigma2 + 1 / 100^2
    delta <- rnorm(1, sum((y - mu)[grp == 1]) / sigma2 / prec, 1 / sqrt(prec))
    sse   <- sum((y - mu - delta * grp)^2)
    sigma2 <- 1 / rgamma(1, 0.01 + N / 2, 0.01 + sse / 2)
    draws[it, ] <- c(mu, delta, sqrt(sigma2))
  }
  draws[-(1:1000), ]
}
fit5 <- gibbs5(y5, grp)
print(t(apply(fit5, 2, post_summary)))
say("sample means: group 1 = %.3f, group 2 - group 1 = %.3f, pooled sd = %.3f",
    mean(y5[grp == 0]), mean(y5[grp == 1]) - mean(y5[grp == 0]),
    sqrt(sum(tapply(y5, grp, function(v) sum((v - mean(v))^2))) / 98))

save_fig("ex05_marginals", width = 7.5, height = 2.8, function() {
  par(mfrow = c(1, 3), mar = c(4, 3, 1.2, 0.5))
  truth <- c(mu = 10, delta = 1, sigma = 2)
  labels <- list(expression(mu), expression(delta), expression(sigma))
  for (j in 1:3) {
    plot(density(fit5[, j]), main = "", xlab = labels[[j]], ylab = "", lwd = 2,
         col = PALETTE["navy"], zero.line = FALSE)
    abline(v = truth[j], lty = 2, col = PALETTE["red"], lwd = 2)
  }
})


# ---- Exercise 6: Metropolis-within-Gibbs, beta prior with random scale --------

section("Ex 6  Metropolis-within-Gibbs")

log_post_b <- function(b, theta) {
  e <- exp(b)
  dnorm(b, 0, 1, log = TRUE) + lgamma(2 * e) - 2 * lgamma(e) + (e - 1) * log(theta * (1 - theta))
}

mwg6 <- function(Y, n, iters = 20000, step = 1) {
  theta <- Y / n; b <- 0; accepted <- 0
  draws <- matrix(NA, iters, 2, dimnames = list(NULL, c("theta", "b")))
  for (it in seq_len(iters)) {
    theta <- rbeta(1, Y + exp(b), n - Y + exp(b))           # Gibbs step
    cand  <- rnorm(1, b, step)                               # Metropolis step
    if (log(runif(1)) < log_post_b(cand, theta) - log_post_b(b, theta)) {
      b <- cand; accepted <- accepted + 1
    }
    draws[it, ] <- c(theta, b)
  }
  attr(draws, "acceptance") <- accepted / iters
  draws
}
fit6 <- mwg6(Y = 10, n = 100)
say("acceptance rate = %.3f", attr(fit6, "acceptance"))
fit6 <- fit6[-(1:1000), ]
print(t(apply(fit6, 2, post_summary)))
say("ESS: theta %.0f   b %.0f", ess(fit6[, "theta"]), ess(fit6[, "b"]))

save_fig("ex06_posterior", width = 7.5, height = 2.8, function() {
  par(mfrow = c(1, 3), mar = c(4, 3.4, 1.2, 0.5))
  plot(density(fit6[, "theta"]), main = "", xlab = expression(theta), lwd = 2,
       col = PALETTE["navy"], zero.line = FALSE)
  plot(density(fit6[, "b"]), main = "", xlab = "b", lwd = 2, col = PALETTE["navy"], zero.line = FALSE)
  curve(dnorm(x), add = TRUE, lty = 2, col = PALETTE["grey"])
  plot(fit6[1:3000, "b"], type = "l", lwd = 0.5, col = PALETTE["navy"], xlab = "iteration", ylab = "b")
})


# ---- Exercise 7: NBA clutch shooting, shrinkage toward the overall rate -------

section("Ex 7  NBA hierarchical beta model")

nba <- data.frame(
  player  = c("Westbrook", "Harden", "Leonard", "James", "Thomas", "Curry",
              "Antetokounmpo", "Wall", "Davis", "Durant"),
  q = c(0.845, 0.847, 0.880, 0.674, 0.909, 0.898, 0.770, 0.801, 0.802, 0.875),
  y = c(64, 72, 55, 27, 75, 24, 28, 66, 40, 13),
  n = c(75, 95, 63, 39, 83, 26, 41, 82, 54, 16))

log_post_m <- function(m, theta, q) {
  e <- exp(m)
  dnorm(m, 0, sqrt(10), log = TRUE) +
    sum(dbeta(theta, e * q, e * (1 - q), log = TRUE))
}

mcmc7 <- function(dat, iters = 30000, step = 0.6) {
  m <- 0; accepted <- 0
  draws <- matrix(NA, iters, nrow(dat) + 1, dimnames = list(NULL, c(dat$player, "m")))
  for (it in seq_len(iters)) {
    e     <- exp(m)
    theta <- rbeta(nrow(dat), dat$y + e * dat$q, dat$n - dat$y + e * (1 - dat$q))   # Gibbs
    cand  <- rnorm(1, m, step)                                                     # Metropolis
    if (log(runif(1)) < log_post_m(cand, theta, dat$q) - log_post_m(m, theta, dat$q)) {
      m <- cand; accepted <- accepted + 1
    }
    draws[it, ] <- c(theta, m)
  }
  attr(draws, "acceptance") <- accepted / iters
  draws
}
fit7 <- mcmc7(nba)
say("acceptance rate for m = %.3f", attr(fit7, "acceptance"))
fit7 <- fit7[-(1:2000), ]
tab7 <- t(apply(fit7, 2, function(x) round(c(mean = mean(x), quantile(x, c(0.025, 0.975))), 3)))
print(cbind(tab7, ESS = round(apply(fit7, 2, ess))))


# ---- Exercise 8: Metropolis for the gamma shape parameter ---------------------

section("Ex 8  Metropolis for Gamma(theta, 1) shape")

log_post8 <- function(theta) ifelse(theta > 0 & theta < 10, dgamma(5, theta, 1, log = TRUE), -Inf)

metropolis8 <- function(iters = 50000, step = 2.5, start = 5) {
  theta <- start; draws <- numeric(iters); accepted <- 0
  for (it in seq_len(iters)) {
    cand <- rnorm(1, theta, step)
    if (log(runif(1)) < log_post8(cand) - log_post8(theta)) {
      theta <- cand; accepted <- accepted + 1
    }
    draws[it] <- theta
  }
  list(draws = draws, acceptance = accepted / iters)
}

for (step in c(0.2, 1, 2.5, 8)) {
  r <- metropolis8(iters = 10000, step = step)
  say("step sd %.1f: acceptance %.3f   ESS %.0f / 10000", step, r$acceptance, ess(r$draws))
}
run8 <- metropolis8()
draws8 <- run8$draws[-(1:1000)]
theta8 <- seq(0.001, 10, by = 0.01)
curve8 <- dgamma(5, theta8, 1); curve8 <- curve8 / sum(curve8 * 0.01)
say("final: acceptance %.3f   mean %.3f (grid %.3f)   95%% CI %s",
    run8$acceptance, mean(draws8), sum(theta8 * curve8) * 0.01,
    fmt_ci(quantile(draws8, c(0.025, 0.975)), 2))

save_fig("ex08_metropolis", width = 7.5, height = 3.2, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.2, 1.2, 0.6))
  hist(draws8, breaks = 60, freq = FALSE, col = adjustcolor(PALETTE["navy"], 0.6), border = "white",
       main = "", xlab = expression(theta), ylab = "Posterior density", xlim = c(0, 10))
  lines(theta8, curve8, lwd = 2.5, col = PALETTE["red"])
  plot(draws8[1:3000], type = "l", lwd = 0.5, col = PALETTE["navy"], xlab = "iteration",
       ylab = expression(theta))
})


# ---- Exercise 9: Student-t model for the galaxies data ------------------------

section("Ex 9  t model for galaxies")

galaxies <- MASS::galaxies
n9 <- length(galaxies)

# Gibbs via the scale-mixture representation Y_i ~ N(mu, sigma^2 / lambda_i),
# lambda_i ~ Gamma(k/2, k/2); k is updated by Metropolis on (1, 30).
gibbs9 <- function(y, iters = 20000, mu = median(y), tau = 1 / var(y), k = 5) {
  n <- length(y); lambda <- rep(1, n)
  log_post_k <- function(k) if (k <= 1 || k >= 30) -Inf else sum(dgamma(lambda, k / 2, k / 2, log = TRUE))
  draws <- matrix(NA, iters, 3, dimnames = list(NULL, c("mu", "sigma", "k")))
  for (it in seq_len(iters)) {
    prec   <- tau * sum(lambda) + 1 / 10000^2
    mu     <- rnorm(1, tau * sum(lambda * y) / prec, 1 / sqrt(prec))
    tau    <- rgamma(1, 0.01 + n / 2, 0.01 + sum(lambda * (y - mu)^2) / 2)
    lambda <- rgamma(n, (k + 1) / 2, (k + tau * (y - mu)^2) / 2)
    cand   <- rnorm(1, k, 3)
    if (log(runif(1)) < log_post_k(cand) - log_post_k(k)) k <- cand
    draws[it, ] <- c(mu, 1 / sqrt(tau), k)
  }
  draws
}

chains9 <- list(gibbs9(galaxies, mu = 15000, tau = 1 / 5000^2, k = 25),
                gibbs9(galaxies, mu = 25000, tau = 1 / 1000^2, k = 3))
keep9 <- -(1:2000)
post9 <- do.call(rbind, lapply(chains9, function(ch) ch[keep9, ]))
print(t(apply(post9, 2, post_summary, digits = 1)))
for (p in colnames(post9))
  say("%-5s  R-hat = %.3f   ESS = %.0f", p,
      rhat(sapply(chains9, function(ch) ch[keep9, p])), sum(sapply(chains9, function(ch) ess(ch[keep9, p]))))

save_fig("ex09_traces", width = 7.5, height = 3, function() {
  par(mfrow = c(1, 3), mar = c(4, 4.6, 1.2, 0.5))
  for (p in colnames(post9)) trace_plot(sapply(chains9, function(ch) ch[, p]), ylab = p)
})

save_fig("ex09_fit", width = 6.5, height = 3.8, function() {
  pm <- colMeans(post9)
  hist(galaxies, breaks = 25, freq = FALSE, col = adjustcolor(PALETTE["navy"], 0.5), border = "white",
       main = "", xlab = "velocity (km/s)", ylab = "Density")
  x <- seq(5000, 40000, length.out = 800)
  lines(x, dt((x - pm["mu"]) / pm["sigma"], pm["k"]) / pm["sigma"], lwd = 2.5, col = PALETTE["red"])
})


# ---- Exercise 10: Laplace model for galaxies on a grid -----------------------

section("Ex 10  Laplace model for galaxies")

mu10    <- seq(19000, 23000, length.out = 300)
sigma10 <- seq(2000, 5500, length.out = 300)
abs_dev <- sapply(mu10, function(m) sum(abs(galaxies - m)))
log_post10 <- outer(abs_dev, sigma10, function(d, s) -n9 * log(s) - d / s)   # flat priors
post10 <- exp(log_post10 - max(log_post10)); post10 <- post10 / sum(post10)
marg_mu <- rowSums(post10); marg_sigma <- colSums(post10)
mu_hat <- sum(mu10 * marg_mu); sigma_hat <- sum(sigma10 * marg_sigma)
say("posterior mean mu = %.1f   sigma = %.1f   (sample median %.1f)", mu_hat, sigma_hat, median(galaxies))

cells <- sample(length(post10), 1e5, replace = TRUE, prob = post10)
mu_s <- mu10[row(post10)[cells]]; sigma_s <- sigma10[col(post10)[cells]]
y_star <- mu_s + sigma_s * (rexp(1e5) - rexp(1e5))         # Laplace draws
say("PPD: mean %.1f  var %.4g    plug-in: mean %.1f  var %.4g   (ratio %.3f)",
    mean(y_star), var(y_star), mu_hat, 2 * sigma_hat^2, var(y_star) / (2 * sigma_hat^2))

dlaplace <- function(x, m, s) exp(-abs(x - m) / s) / (2 * s)

save_fig("ex10_posterior", width = 7.5, height = 2.9, function() {
  par(mfrow = c(1, 3), mar = c(4, 4.4, 1.2, 0.5))
  contour(mu10, sigma10, post10, drawlabels = FALSE, col = PALETTE["navy"],
          xlab = expression(mu), ylab = expression(sigma))
  plot(mu10, marg_mu / diff(mu10)[1], type = "l", lwd = 2, col = PALETTE["navy"],
       xlab = expression(mu), ylab = "density")
  plot(sigma10, marg_sigma / diff(sigma10)[1], type = "l", lwd = 2, col = PALETTE["navy"],
       xlab = expression(sigma), ylab = "density")
})

save_fig("ex10_fit", width = 6.5, height = 3.8, function() {
  hist(galaxies, breaks = 25, freq = FALSE, col = adjustcolor(PALETTE["navy"], 0.5), border = "white",
       main = "", xlab = "velocity (km/s)", ylab = "Density", ylim = c(0, 1.6e-4))
  x <- seq(5000, 40000, length.out = 800)
  lines(x, dlaplace(x, mu_hat, sigma_hat), lwd = 2.5, col = PALETTE["red"])
  lines(density(y_star), lwd = 2, lty = 2, col = PALETTE["gold"])
  legend("topright", c("plug-in Laplace", "posterior predictive"), lty = 1:2, lwd = 2,
         col = PALETTE[c("red", "gold")], bty = "n")
})


# ---- Exercise 11: Reggie Jackson again, with Metropolis + diagnostics ---------

section("Ex 11  Reggie Jackson with diagnostics")

log_post11 <- function(log_lam, hr, games) {        # Uniform(0, 10) prior, sampled on log scale
  lam <- exp(log_lam)
  if (lam >= 10) return(-Inf)
  hr * log_lam - games * lam + log_lam               # + log-Jacobian
}
metro11 <- function(hr, games, start, iters = 20000, step) {
  x <- log(start); draws <- numeric(iters); acc <- 0
  for (it in seq_len(iters)) {
    cand <- rnorm(1, x, step)
    if (log(runif(1)) < log_post11(cand, hr, games) - log_post11(x, hr, games)) { x <- cand; acc <- acc + 1 }
    draws[it] <- exp(x)
  }
  structure(draws, acceptance = acc / iters)
}
starts <- c(0.05, 2)
reg <- sapply(starts, function(s) metro11(563, 2820, s, step = 0.1))
ws  <- sapply(starts, function(s) metro11(10, 27, s, step = 0.7))
keep11 <- -(1:2000)
ratio11 <- ws[keep11, ] / reg[keep11, ]
for (nm in c("regular", "world series", "ratio")) {
  ch <- switch(nm, regular = reg[keep11, ], `world series` = ws[keep11, ], ratio = ratio11)
  s <- post_summary(as.vector(ch))
  say("%-13s mean %.4f  sd %.4f  95%% CI (%.4f, %.4f)   R-hat %.4f   ESS %.0f",
      nm, s[1], s[2], s[3], s[4], rhat(ch), ess(ch[, 1]) + ess(ch[, 2]))
}
say("P(ratio > 1) = %.4f", mean(ratio11 > 1))

save_fig("ex11_traces", width = 7.5, height = 2.9, function() {
  par(mfrow = c(1, 3), mar = c(4, 4.4, 1.6, 0.5))
  trace_plot(reg[1:4000, ], ylab = expression(lambda[R]), main = "regular season")
  trace_plot(ws[1:4000, ], ylab = expression(lambda[WS]), main = "World Series")
  trace_plot((ws / reg)[1:4000, ], ylab = "ratio", main = "ratio")
})


# ---- Exercises 12-13: Poisson regression with a hand-written Metropolis ------

section("Ex 12-13  Poisson log-linear trend")

y12 <- c(64, 13, 33, 18, 30, 20); t12 <- 1:6
log_post12 <- function(ab) {
  eta <- ab[1] + ab[2] * t12
  sum(y12 * eta - exp(eta)) + sum(dnorm(ab, 0, 10, log = TRUE))
}

# One-at-a-time random-walk Metropolis; step sizes tuned for ~40% acceptance.
metro12 <- function(start, steps = c(0.12, 0.035), iters = 30000) {
  ab <- start; lp <- log_post12(ab)
  draws <- matrix(NA, iters, 2, dimnames = list(NULL, c("alpha", "beta")))
  acc <- c(0, 0)
  for (it in seq_len(iters)) {
    for (j in 1:2) {
      cand <- ab; cand[j] <- rnorm(1, ab[j], steps[j])
      lp_cand <- log_post12(cand)
      if (log(runif(1)) < lp_cand - lp) { ab <- cand; lp <- lp_cand; acc[j] <- acc[j] + 1 }
    }
    draws[it, ] <- ab
  }
  structure(draws, acceptance = acc / iters)
}
chains12 <- list(metro12(c(0, 0)), metro12(c(5, -1)), metro12(c(2, 1)))
say("acceptance: alpha %.3f   beta %.3f", attr(chains12[[1]], "acceptance")[1], attr(chains12[[1]], "acceptance")[2])
keep12 <- -(1:3000)
post12 <- do.call(rbind, lapply(chains12, function(ch) ch[keep12, ]))
print(t(apply(post12, 2, post_summary)))
for (p in c("alpha", "beta"))
  say("%-5s R-hat %.4f   ESS %.0f", p, rhat(sapply(chains12, function(ch) ch[keep12, p])),
      sum(sapply(chains12, function(ch) ess(ch[keep12, p]))))
say("P(beta < 0) = %.4f   posterior median yearly change = %.1f%%",
    mean(post12[, "beta"] < 0), 100 * (exp(median(post12[, "beta"])) - 1))
say("cor(alpha, beta) = %.3f", cor(post12[, 1], post12[, 2]))
glm12 <- coef(glm(y12 ~ t12, family = poisson))
say("glm MLE: alpha %.4f  beta %.4f", glm12[1], glm12[2])

save_fig("ex12_traces", width = 7.5, height = 3, function() {
  par(mfrow = c(1, 3), mar = c(4, 4.4, 1.2, 0.5))
  trace_plot(sapply(chains12, function(ch) ch[1:6000, "alpha"]), ylab = expression(alpha))
  trace_plot(sapply(chains12, function(ch) ch[1:6000, "beta"]), ylab = expression(beta))
  plot(density(post12[, "beta"]), main = "", xlab = expression(beta), lwd = 2,
       col = PALETTE["navy"], zero.line = FALSE)
  abline(v = 0, lty = 2, col = PALETTE["red"])
})


# ---- Exercise 14: three-arm clinical trial ------------------------------------

section("Ex 14  three-arm trial")

arms <- data.frame(arm = c("placebo", "low", "high"), pos = c(52, 60, 54), neg = c(48, 40, 46))
draws14 <- sapply(seq_len(nrow(arms)), function(i) rbeta(1e6, 1 + arms$pos[i], 1 + arms$neg[i]))
colnames(draws14) <- arms$arm
print(t(apply(draws14, 2, post_summary)))
best <- table(factor(colnames(draws14)[max.col(draws14)], levels = arms$arm)) / nrow(draws14)
print(round(best, 4))


# ---- Exercise 15: normal mixture with vague proper priors ---------------------

section("Ex 15  two-component normal mixture")

set.seed(27695)
theta_true <- 4; n15 <- 30
B <- rbinom(n15, 1, 0.5)
Y15 <- rnorm(n15, B * theta_true, 1)
set.seed(2026)

nlp <- function(theta, Y) {
  like  <- 0.5 * dnorm(Y, 0, 1) + 0.5 * dnorm(Y, theta, 1)
  prior <- dnorm(theta, 0, 10)
  -sum(log(like)) - sum(log(prior))
}
map15 <- suppressWarnings(optim(mean(Y15), nlp, Y = Y15, hessian = TRUE))
map_est <- map15$par; map_sd <- 1 / sqrt(map15$hessian[1, 1])
say("MAP = %.4f   asymptotic sd = %.4f", map_est, map_sd)

theta15 <- seq(-2, 8, length.out = 4001)
loglik15 <- sapply(theta15, function(th) sum(log(0.5 * dnorm(Y15) + 0.5 * dnorm(Y15, th))))
post_k <- sapply(0:3, function(k) {
  lp <- loglik15 + dnorm(theta15, 0, sqrt(10^k), log = TRUE)
  d <- exp(lp - max(lp)); d / sum(d * diff(theta15)[1])
})
for (k in 0:3) say("prior N(0, 10^%d): posterior mean %.4f  sd %.4f", k,
                   sum(theta15 * post_k[, k + 1]) * diff(theta15)[1],
                   sqrt(sum(theta15^2 * post_k[, k + 1]) * diff(theta15)[1] -
                          (sum(theta15 * post_k[, k + 1]) * diff(theta15)[1])^2))

save_fig("ex15_data", width = 7.5, height = 3, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.2, 1.2, 0.5))
  hist(Y15, breaks = 15, col = adjustcolor(PALETTE["navy"], 0.6), border = "white", main = "",
       xlab = "y", ylab = "count")
  y <- seq(-3, 10, length.out = 500)
  plot(NA, xlim = c(-3, 10), ylim = c(0, 0.22), xlab = "y", ylab = expression(f(y ~ "|" ~ theta)))
  for (k in 1:3) lines(y, 0.5 * dnorm(y) + 0.5 * dnorm(y, c(2, 4, 6)[k]), lwd = 2, col = PALETTE[k])
  legend("topright", paste("theta =", c(2, 4, 6)), col = PALETTE[1:3], lwd = 2, bty = "n")
})

save_fig("ex15_posteriors", width = 6.5, height = 4, function() {
  plot(NA, xlim = c(2, 6), ylim = c(0, max(post_k) * 1.05), xlab = expression(theta),
       ylab = "Posterior density")
  for (k in 0:3) lines(theta15, post_k[, k + 1], lwd = 2, col = PALETTE[k + 1])
  lines(theta15, dnorm(theta15, map_est, map_sd), lwd = 2, lty = 2, col = "black")
  legend("topright", c(paste0("prior N(0, 10^", 0:3, ")"), "asymptotic normal"),
         col = c(PALETTE[1:4], "black"), lty = c(1, 1, 1, 1, 2), lwd = 2, bty = "n")
})

# (e) Gibbs sampler via the latent-allocation representation.
gibbs15 <- function(Y, iters = 20000, theta = mean(Y)) {
  draws <- numeric(iters)
  for (it in seq_len(iters)) {
    p1 <- dnorm(Y, theta, 1); p0 <- dnorm(Y, 0, 1)
    B  <- rbinom(length(Y), 1, p1 / (p0 + p1))
    prec  <- sum(B) + 1 / 100
    theta <- rnorm(1, sum(Y[B == 1]) / prec, 1 / sqrt(prec))
    draws[it] <- theta
  }
  draws[-(1:1000)]
}
g15 <- gibbs15(Y15)
say("Gibbs (latent B): mean %.4f  sd %.4f   ESS %.0f", mean(g15), sd(g15), ess(g15))


# ---- Exercise 16: shifted exponential, demonstration of the sampler ----------

section("Ex 16  shifted exponential")

save_fig("ex16_density", width = 6, height = 3.6, function() {
  y <- seq(0, 5, length.out = 1001)
  plot(y, ifelse(y >= 2, 3 * exp(-3 * (y - 2)), 0), type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = "y", ylab = expression(f(y ~ "|" ~ alpha == 2 ~ "," ~ beta == 3)))
})

# Demonstration on simulated data: alpha = 2, beta = 3, n = 40, priors U(0, 5) and Gamma(1, 1).
y16 <- 2 + rexp(40, 3)
gibbs16 <- function(y, c = 5, a = 1, b = 1, iters = 20000) {
  n <- length(y); upper <- min(c, min(y))
  alpha <- upper / 2; beta <- 1
  draws <- matrix(NA, iters, 2, dimnames = list(NULL, c("alpha", "beta")))
  for (it in seq_len(iters)) {
    rate  <- n * beta                         # alpha | beta has density ∝ exp(rate * alpha) on (0, upper)
    u     <- runif(1)                         # inverse-CDF draw from that truncated exponential
    alpha <- upper + log(u + (1 - u) * exp(-rate * upper)) / rate
    beta  <- rgamma(1, a + n, b + sum(y - alpha))
    draws[it, ] <- c(alpha, beta)
  }
  draws[-(1:1000), ]
}
fit16 <- gibbs16(y16)
print(t(apply(fit16, 2, post_summary)))
say("min(y) = %.4f", min(y16))


# ---- Exercise 17: nonlinear T. rex growth curve -------------------------------

section("Ex 17  nonlinear growth curve")

mass <- c(29.9, 1761, 1807, 2984, 3230, 5040, 5654)
age  <- c(2, 15, 14, 16, 18, 22, 28)
n17  <- length(mass)

rtruncnorm <- function(m, s, lo, hi) {
  u <- runif(1, pnorm(lo, m, s), pnorm(hi, m, s))
  qnorm(u, m, s)
}

mcmc17 <- function(iters = 40000, theta3 = 1, step = 0.08) {
  theta1 <- 0; sigma2 <- 1e6
  x <- age^theta3; theta2 <- min(19999, max(1, sum(x * mass) / sum(x^2)))
  log_post3 <- function(t3, t1, t2, s2) {
    sum(dnorm(mass, t1 + t2 * age^t3, sqrt(s2), log = TRUE)) + dnorm(t3, 0, 1, log = TRUE)
  }
  draws <- matrix(NA, iters, 4, dimnames = list(NULL, c("theta1", "theta2", "theta3", "sigma")))
  acc <- 0
  for (it in seq_len(iters)) {
    x <- age^theta3
    prec   <- n17 / sigma2 + 1 / 100^2
    theta1 <- rnorm(1, sum(mass - theta2 * x) / sigma2 / prec, 1 / sqrt(prec))
    prec   <- sum(x^2) / sigma2
    theta2 <- rtruncnorm(sum(x * (mass - theta1)) / sigma2 / prec, 1 / sqrt(prec), 0, 20000)
    sigma2 <- 1 / rgamma(1, 0.01 + n17 / 2, 0.01 + sum((mass - theta1 - theta2 * x)^2) / 2)
    cand <- rnorm(1, theta3, step)
    if (log(runif(1)) < log_post3(cand, theta1, theta2, sigma2) - log_post3(theta3, theta1, theta2, sigma2)) {
      theta3 <- cand; acc <- acc + 1
    }
    draws[it, ] <- c(theta1, theta2, theta3, sqrt(sigma2))
  }
  structure(draws, acceptance = acc / iters)
}
chains17 <- lapply(c(0.5, 1, 1.5, 2.5), function(t3) mcmc17(theta3 = t3))
say("acceptance for theta3: %s", paste(sprintf("%.2f", sapply(chains17, attr, "acceptance")), collapse = " "))
keep17 <- -(1:5000)
post17 <- do.call(rbind, lapply(chains17, function(ch) ch[keep17, ]))
print(t(apply(post17, 2, post_summary, digits = 3)))
for (p in colnames(post17))
  say("%-6s R-hat %.3f   total ESS %.0f   lag-1 acf %.3f", p,
      rhat(sapply(chains17, function(ch) ch[keep17, p])),
      sum(sapply(chains17, function(ch) ess(ch[keep17, p]))),
      acf(chains17[[1]][keep17, p], plot = FALSE)$acf[2])
say("cor(theta2, theta3) = %.3f", cor(post17[, "theta2"], post17[, "theta3"]))

save_fig("ex17_fit", width = 6, height = 3.8, function() {
  ages <- seq(0, 30, length.out = 200)
  idx  <- sample(nrow(post17), 4000)
  curves <- sapply(ages, function(a) post17[idx, "theta1"] + post17[idx, "theta2"] * a^post17[idx, "theta3"])
  band <- apply(curves, 2, quantile, c(0.025, 0.975))
  plot(age, mass, pch = 19, col = PALETTE["navy"], cex = 1.2, xlim = c(0, 30), ylim = c(-1000, 8000),
       xlab = "age (years)", ylab = "mass (kg)")
  polygon(c(ages, rev(ages)), c(band[1, ], rev(band[2, ])), col = adjustcolor(PALETTE["gold"], 0.25), border = NA)
  lines(ages, colMeans(curves), lwd = 2.5, col = PALETTE["red"])
})

save_fig("ex17_traces", width = 7.5, height = 4.6, function() {
  par(mfrow = c(2, 2), mar = c(4, 4.6, 1.2, 0.5))
  for (p in colnames(post17)) trace_plot(sapply(chains17, function(ch) ch[1:20000, p]), ylab = p)
})


# ---- Exercise 18: binomial with unknown n -------------------------------------

section("Ex 18  unknown binomial n")

gibbs18 <- function(Y, lambda, a, b, start_p, iters = 20000) {
  p <- start_p
  draws <- matrix(NA, iters, 3, dimnames = list(NULL, c("n", "p", "theta")))
  for (it in seq_len(iters)) {
    n <- Y + rpois(1, lambda * (1 - p))            # n | p, Y
    p <- rbeta(1, a + Y, b + n - Y)                # p | n, Y
    draws[it, ] <- c(n, p, n * p)
  }
  draws
}
fits18 <- list()
for (ab in c(1, 10)) {
  chains <- lapply(c(0.05, 0.5, 0.95), function(p0) gibbs18(10, 10, ab, ab, p0))
  fits18[[paste0("a=b=", ab)]] <- chains
  post <- do.call(rbind, lapply(chains, function(ch) ch[-(1:1000), ]))
  say("a = b = %d", ab)
  for (p in colnames(post)) {
    s <- post_summary(post[, p], 3)
    say("   %-5s mean %7.3f  sd %6.3f  95%% CI (%.3f, %.3f)   R-hat %.4f  ESS %6.0f  lag-1 acf %.3f",
        p, s[1], s[2], s[3], s[4], rhat(sapply(chains, function(ch) ch[-(1:1000), p])),
        sum(sapply(chains, function(ch) ess(ch[-(1:1000), p]))),
        acf(chains[[1]][-(1:1000), p], plot = FALSE)$acf[2])
  }
}

save_fig("ex18_traces", width = 7.5, height = 4.4, function() {
  par(mfrow = c(2, 3), mar = c(4, 4.2, 1.8, 0.5))
  for (nm in names(fits18)) for (p in c("n", "p", "theta"))
    trace_plot(sapply(fits18[[nm]], function(ch) ch[1:1500, p]), ylab = p, main = nm)
})


# ---- Exercise 19: an over-parameterised hierarchical model --------------------

section("Ex 19  over-parameterised hierarchical model")

gibbs19 <- function(n, eps, iters = 20000, step = 0.5) {
  y <- 1:n
  mu <- y; sigma2 <- rep(1, n); th <- c(1, 1, 1)
  log_post_th2 <- function(t2) {              # sampled on the log scale
    if (!is.finite(t2) || t2 <= 0) return(-Inf)
    dgamma(t2, eps, eps, log = TRUE) + sum(t2 * log(th[3]) - lgamma(t2) - (t2 + 1) * log(sigma2)) + log(t2)
  }
  draws <- matrix(NA, iters, 3, dimnames = list(NULL, paste0("theta", 1:3)))
  for (it in seq_len(iters)) {
    prec   <- 1 / sigma2 + th[1]
    mu     <- rnorm(n, (y / sigma2) / prec, 1 / sqrt(prec))
    sigma2 <- 1 / rgamma(n, th[2] + 0.5, th[3] + (y - mu)^2 / 2)
    th[1]  <- rgamma(1, eps + n / 2, eps + sum(mu^2) / 2)
    th[3]  <- rgamma(1, eps + n * th[2], eps + sum(1 / sigma2))
    cand   <- th[2] * exp(rnorm(1, 0, step))
    if (log(runif(1)) < log_post_th2(cand) - log_post_th2(th[2])) th[2] <- cand
    draws[it, ] <- th
  }
  draws[-(1:2000), ]
}
ess19 <- t(sapply(list(c(5, 0.1), c(5, 10), c(25, 0.1), c(25, 10)), function(s) {
  d <- gibbs19(s[1], s[2])
  c(n = s[1], eps = s[2], round(apply(d, 2, ess)), round(apply(d, 2, median), 3))
}))
colnames(ess19) <- c("n", "eps", "ESS_t1", "ESS_t2", "ESS_t3", "med_t1", "med_t2", "med_t3")
print(ess19)
