# =============================================================================
#  Chapter 7 -- Case studies using hierarchical modeling
#  Reproduces every number and figure in ch07/ch07.tex with hand-written samplers.
#  The JAGS versions (Exercises 3b, 4 and 6) are in ch07/code/ch07_jags.R and are
#  run by GitHub Actions; this script writes the small data files they need.
#
#  Run from the repository root:   Rscript ch07/code/ch07.R
#  Needs base R + survival (ships with R, for the lung data).
# =============================================================================

source(file.path("common", "helpers.R"))
set.seed(2026)
FIG_DIR  <- file.path("ch07", "figures")
CODE_DIR <- file.path("ch07", "code")
rinvwishart <- function(df, scale) solve(rWishart(1, df, solve(scale))[, , 1])


# ---- Exercise 1: incompatible full conditionals --------------------------------------

section("Ex 1  incompatible conditionals")

# Gibbs sampler that alternates X | Y ~ N(bY, 1) and Y | X ~ N(aX, 1).
gibbs_pair <- function(a, b, iters = 2e5, x = 0) {
  out <- matrix(NA, iters, 2, dimnames = list(NULL, c("x", "y")))
  for (it in seq_len(iters)) {
    y <- rnorm(1, a * x, 1); x <- rnorm(1, b * y, 1)
    out[it, ] <- c(x, y)
  }
  out
}
g1 <- gibbs_pair(1, 0.5)[-(1:1000), ]
fit_xy <- coef(lm(x ~ y, data = as.data.frame(g1))); fit_yx <- coef(lm(y ~ x, data = as.data.frame(g1)))
say("a = 1, b = 0.5:  slope of X on Y in the output = %.3f (b = 0.5),  slope of Y on X = %.3f (a = 1)",
    fit_xy[2], fit_yx[2])
say("theory for the sampler's stationary law: slope of Y on X = b(1+a^2)/(1+b^2) = %.3f (not a = 1)", 0.5 * 2 / 1.25)
g2 <- gibbs_pair(2, 2, iters = 30)
say("a = b = 2: |x| after 5, 10, 20, 30 sweeps: %s",
    paste(sprintf("%.3g", abs(g2[c(5, 10, 20, 30), "x"])), collapse = ", "))
g3 <- gibbs_pair(0.5, 0.5)[-(1:1000), ]
say("a = b = 0.5 (compatible): slopes %.3f and %.3f, Var(X) = %.3f (theory 1/(1-ab) = %.3f)",
    coef(lm(x ~ y, data = as.data.frame(g3)))[2], coef(lm(y ~ x, data = as.data.frame(g3)))[2],
    var(g3[, "x"]), 1 / (1 - 0.25))


# ---- Exercise 3: meta-analysis of corticosteroid trials --------------------------------

section("Ex 3  Cochrane meta-analysis")

coch <- data.frame(name = c("Auckland", "Block", "Doran", "Gamsu", "Morrison", "Papageorgiou", "Tauesch"),
                   ev_trt = c(36, 1, 4, 14, 3, 1, 8), n_trt = c(532, 69, 81, 131, 67, 71, 56),
                   ev_ctrl = c(60, 5, 11, 20, 7, 7, 10), n_ctrl = c(538, 61, 63, 137, 59, 75, 71))
write.csv(coch, file.path(CODE_DIR, "cochrane.csv"), row.names = FALSE)
K3 <- nrow(coch)

# (a) pooled: theta_0, theta_1 ~ Uniform(0, 1)
S <- 1e6
th0 <- rbeta(S, 1 + sum(coch$ev_ctrl), 1 + sum(coch$n_ctrl - coch$ev_ctrl))
th1 <- rbeta(S, 1 + sum(coch$ev_trt), 1 + sum(coch$n_trt - coch$ev_trt))
say("(a) pooled rates: control %.4f, treatment %.4f   P(theta1 < theta0) = %.5f",
    mean(th0), mean(th1), mean(th1 < th0))
say("    odds ratio: median %.3f  95%% CI %s", median(th1 / (1 - th1) / (th0 / (1 - th0))),
    fmt_ci(quantile(th1 / (1 - th1) / (th0 / (1 - th0)), c(0.025, 0.975))))

# (b) bivariate random effects on the logits: alpha_i ~ N2(mu, Sigma),
#     mu ~ N(0, 10^2 I), Sigma ~ InvWishart(3, I).
Yc <- cbind(coch$ev_ctrl, coch$ev_trt); Nc <- cbind(coch$n_ctrl, coch$n_trt)
loglik_trial <- function(alpha) rowSums(dbinom(Yc, Nc, plogis(alpha), log = TRUE))
dmvn_rows <- function(A, mu, Sig_inv) -0.5 * rowSums((sweep(A, 2, mu) %*% Sig_inv) * sweep(A, 2, mu))

meta_mcmc <- function(iters = 40000, burn = 5000) {
  alpha <- qlogis((Yc + 0.5) / (Nc + 1)); mu <- colMeans(alpha); Sig <- diag(2)
  out <- matrix(NA, iters, 2 + 3 + 2 * K3,
                dimnames = list(NULL, c("mu0", "mu1", "s00", "s11", "s01", paste0("a0_", 1:K3), paste0("a1_", 1:K3))))
  acc <- 0
  for (it in seq_len(iters)) {
    Sig_inv <- solve(Sig)
    cand <- alpha + matrix(rnorm(2 * K3, 0, 0.35), K3)              # MH for every trial at once
    lr <- loglik_trial(cand) + dmvn_rows(cand, mu, Sig_inv) - loglik_trial(alpha) - dmvn_rows(alpha, mu, Sig_inv)
    take <- log(runif(K3)) < lr
    alpha[take, ] <- cand[take, ]; acc <- acc + mean(take)
    V <- solve(K3 * Sig_inv + diag(1 / 100, 2))                     # mu | alpha, Sigma
    mu <- drop(V %*% Sig_inv %*% colSums(alpha) + t(chol(V)) %*% rnorm(2))
    D <- sweep(alpha, 2, mu)
    Sig <- rinvwishart(3 + K3, diag(2) + crossprod(D))              # Sigma | alpha, mu
    out[it, ] <- c(mu, Sig[1, 1], Sig[2, 2], Sig[1, 2], alpha[, 1], alpha[, 2])
  }
  structure(out[-seq_len(burn), ], acceptance = acc / iters)
}
chains3 <- lapply(1:3, function(k) meta_mcmc())
post3 <- do.call(rbind, chains3)
say("(b) acceptance for alpha: %.2f", attr(chains3[[1]], "acceptance"))
lor <- post3[, "mu1"] - post3[, "mu0"]
say("    mean log odds ratio mu1 - mu0: mean %.3f  95%% CI %s   P(< 0) = %.4f   OR median %.3f",
    mean(lor), fmt_ci(quantile(lor, c(0.025, 0.975))), mean(lor < 0), exp(median(lor)))
say("    R-hat(mu1 - mu0) = %.4f   ESS = %.0f", rhat(sapply(chains3, function(ch) ch[, "mu1"] - ch[, "mu0"])),
    sum(sapply(chains3, function(ch) ess(ch[, "mu1"] - ch[, "mu0"]))))
trial_lor <- post3[, paste0("a1_", 1:K3)] - post3[, paste0("a0_", 1:K3)]
raw_lor <- log((coch$ev_trt + 0.5) / (coch$n_trt - coch$ev_trt + 0.5)) -
           log((coch$ev_ctrl + 0.5) / (coch$n_ctrl - coch$ev_ctrl + 0.5))
print(data.frame(trial = coch$name, raw_logOR = round(raw_lor, 2), post_mean = round(colMeans(trial_lor), 2),
                 lo = round(apply(trial_lor, 2, quantile, 0.025), 2), hi = round(apply(trial_lor, 2, quantile, 0.975), 2),
                 P_neg = round(colMeans(trial_lor < 0), 3)), row.names = FALSE)
say("    SD of baseline logits sqrt(s00): %.2f   SD of treatment logits sqrt(s11): %.2f",
    mean(sqrt(post3[, "s00"])), mean(sqrt(post3[, "s11"])))

# WAIC for both models (pointwise over the 14 trial-arm counts)
ll_a <- cbind(sapply(1:K3, function(i) dbinom(coch$ev_ctrl[i], coch$n_ctrl[i], th0[1:20000], log = TRUE)),
              sapply(1:K3, function(i) dbinom(coch$ev_trt[i], coch$n_trt[i], th1[1:20000], log = TRUE)))
idx3 <- sample(nrow(post3), 20000)
ll_b <- cbind(sapply(1:K3, function(i) dbinom(coch$ev_ctrl[i], coch$n_ctrl[i], plogis(post3[idx3, paste0("a0_", i)]), log = TRUE)),
              sapply(1:K3, function(i) dbinom(coch$ev_trt[i], coch$n_trt[i], plogis(post3[idx3, paste0("a1_", i)]), log = TRUE)))
print(round(rbind(pooled = waic(ll_a), random_effects = waic(ll_b)), 2))

save_fig("ex03_forest", width = 6.5, height = 4, function() {
  q <- apply(trial_lor, 2, quantile, c(0.025, 0.5, 0.975)); y <- K3:1
  par(mar = c(4.2, 7.5, 1, 1))
  plot(NA, xlim = c(-3.5, 2), ylim = c(-0.5, K3 + 0.5), yaxt = "n", xlab = "log odds ratio (treatment vs control)", ylab = "")
  axis(2, at = c(y, -0.2), labels = c(coch$name, "overall"), las = 1)
  abline(v = 0, lty = 2, col = PALETTE["grey"])
  segments(q[1, ], y, q[3, ], y, lwd = 3, col = PALETTE["navy"]); points(q[2, ], y, pch = 19, col = PALETTE["navy"])
  points(raw_lor, y + 0.25, pch = 4, col = PALETTE["red"])
  qo <- quantile(lor, c(0.025, 0.5, 0.975))
  segments(qo[1], -0.2, qo[3], -0.2, lwd = 5, col = PALETTE["gold"]); points(qo[2], -0.2, pch = 18, cex = 1.8, col = PALETTE["gold"])
  legend("bottomright", c("trial posterior (RE model)", "raw estimate", "overall effect"), bty = "n", cex = 0.8,
         col = c(PALETTE["navy"], PALETTE["red"], PALETTE["gold"]), pch = c(19, 4, 18), lwd = c(3, NA, 5))
})


# ---- Exercise 4: hierarchical AR(1) model for marathon pacing ------------------------

section("Ex 4  marathon pacing model")

mar <- new.env(); load(file.path("data", "BostonMarathon2016.Rdata"), envir = mar)
Ysp <- mar$Marathon_female$SPEED_mile
Ysp <- Ysp[complete.cases(Ysp), ]
write.csv(Ysp, file.path(CODE_DIR, "marathon_female_complete.csv"), row.names = FALSE)
n4 <- nrow(Ysp); J4 <- ncol(Ysp)
say("runners %d, miles %d, mean speed %.3f min/mile", n4, J4, mean(Ysp))

marathon_gibbs <- function(iters = 25000, burn = 5000, theta = c(6, 1, 0.5, 0.1, 2, 0.1), common_rho = FALSE) {
  mu <- rowMeans(Ysp); rho <- rep(0.5, n4); s2 <- rep(0.1, n4); s20 <- 0.5
  Yc <- Ysp[, -1]; Yl <- Ysp[, -J4]
  log_post_t5 <- function(t5, t6) {               # theta5 sampled on the log scale
    if (t5 <= 0) return(-Inf)
    dgamma(t5, 0.1, 0.1, log = TRUE) + sum(t5 * log(t6) - lgamma(t5) - (t5 + 1) * log(s2)) + log(t5)
  }
  out <- matrix(NA, iters, 6 + 1 + 3 * n4,
                dimnames = list(NULL, c(paste0("theta", 1:6), "sigma0", paste0("mu", 1:n4), paste0("rho", 1:n4), paste0("s2_", 1:n4))))
  for (it in seq_len(iters)) {
    # mu_i: first mile N(mu_i, s20) and miles 2..26 linear in mu_i with coefficient (1 - rho_i)
    w <- 1 - rho
    prec <- 1 / s20 + (J4 - 1) * w^2 / s2 + 1 / theta[2]
    lin  <- Ysp[, 1] / s20 + w * rowSums(Yc - rho * Yl) / s2 + theta[1] / theta[2]
    mu   <- rnorm(n4, lin / prec, 1 / sqrt(prec))
    # rho_i: regression of e_ij on e_{i,j-1}
    E <- Yc - mu; El <- Yl - mu
    if (common_rho) {                              # simplified model: one rho for every runner
      prec <- sum(rowSums(El^2) / s2) + 1 / 10^2
      rho  <- rep(rnorm(1, sum(rowSums(E * El) / s2) / prec, 1 / sqrt(prec)), n4)
    } else {
      prec <- rowSums(El^2) / s2 + 1 / theta[4]
      rho  <- rnorm(n4, (rowSums(E * El) / s2 + theta[3] / theta[4]) / prec, 1 / sqrt(prec))
    }
    s2   <- 1 / rgamma(n4, theta[5] + (J4 - 1) / 2, theta[6] + rowSums((E - rho * El)^2) / 2)
    s20  <- 1 / rgamma(1, 0.1 + n4 / 2, 0.1 + sum((Ysp[, 1] - mu)^2) / 2)
    # hyperparameters
    p1 <- n4 / theta[2] + 1 / 100^2;  theta[1] <- rnorm(1, sum(mu) / theta[2] / p1, 1 / sqrt(p1))
    theta[2] <- 1 / rgamma(1, 0.1 + n4 / 2, 0.1 + sum((mu - theta[1])^2) / 2)
    if (common_rho) { theta[3] <- rho[1]; theta[4] <- 0 } else {
      p3 <- n4 / theta[4] + 1 / 10^2;   theta[3] <- rnorm(1, sum(rho) / theta[4] / p3, 1 / sqrt(p3))
      theta[4] <- 1 / rgamma(1, 0.1 + n4 / 2, 0.1 + sum((rho - theta[3])^2) / 2)
    }
    theta[6] <- rgamma(1, 0.1 + n4 * theta[5], 0.1 + sum(1 / s2))
    cand <- theta[5] * exp(rnorm(1, 0, 0.3))
    if (log(runif(1)) < log_post_t5(cand, theta[6]) - log_post_t5(theta[5], theta[6])) theta[5] <- cand
    out[it, ] <- c(theta, sqrt(s20), mu, rho, s2)
  }
  out[-seq_len(burn), ]
}
chains4 <- lapply(list(c(6, 1, 0.5, 0.1, 2, 0.1), c(5, 0.1, 0.9, 0.01, 0.5, 0.01), c(7, 3, 0, 0.5, 5, 1)),
                  function(th0) marathon_gibbs(theta = th0))
post4 <- do.call(rbind, chains4)
thetas <- c(paste0("theta", 1:6), "sigma0")
diag4 <- t(sapply(thetas, function(p) {
  ch <- sapply(chains4, function(c_) c_[, p])
  c(post_summary(as.vector(ch), 4), Rhat = rhat(ch), ESS = sum(apply(ch, 2, ess)))
}))
print(round(diag4, 4))
say("posterior means of rho_i: range %.2f to %.2f;  sd of mu_i means %.3f",
    min(colMeans(post4[, paste0("rho", 1:n4)])), max(colMeans(post4[, paste0("rho", 1:n4)])),
    sd(colMeans(post4[, paste0("mu", 1:n4)])))
# prior vs posterior spread for (d): prior SDs implied by the hyperpriors
prior_sd <- c(theta1 = 100, theta3 = 10)
say("theta1: prior sd %.0f -> posterior sd %.3f;  theta3: prior sd %.0f -> posterior sd %.3f",
    prior_sd[1], sd(post4[, "theta1"]), prior_sd[2], sd(post4[, "theta3"]))
say("InvGamma(0.1,0.1) / Gamma(0.1,0.1) priors have infinite / very large spread; posterior 95%% widths: theta2 %.3f, theta4 %.3f, theta5 %.2f, theta6 %.3f",
    diff(quantile(post4[, "theta2"], c(.025, .975))), diff(quantile(post4[, "theta4"], c(.025, .975))),
    diff(quantile(post4[, "theta5"], c(.025, .975))), diff(quantile(post4[, "theta6"], c(.025, .975))))

save_fig("ex04_traces", width = 7.5, height = 5, function() {
  par(mfrow = c(2, 3), mar = c(4, 4.4, 1.6, 0.5))
  for (p in paste0("theta", 1:6)) trace_plot(sapply(chains4, function(c_) c_[seq(1, 20000, by = 10), p]), ylab = p)
})

save_fig("ex04_runners", width = 7.5, height = 3.4, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.4, 1.2, 0.6))
  matplot(t(Ysp[1:8, ]), type = "l", lty = 1, col = adjustcolor(PALETTE, 0.8), lwd = 1.5,
          xlab = "mile", ylab = "pace (minutes per mile)")
  plot(colMeans(post4[, paste0("mu", 1:n4)]), colMeans(post4[, paste0("rho", 1:n4)]), pch = 19,
       col = PALETTE["navy"], xlab = expression("posterior mean of " * mu[i]),
       ylab = expression("posterior mean of " * rho[i]))
})

# (e) simplified model with common rho: compare by WAIC using the one-step-ahead likelihood
waic_marathon <- function(post) {
  idx <- sample(nrow(post), 4000)
  ll <- matrix(NA, length(idx), n4 * J4)
  for (k in seq_along(idx)) {
    s <- post[idx[k], ]
    mu <- s[paste0("mu", 1:n4)]; rho <- s[paste0("rho", 1:n4)]; s2 <- s[paste0("s2_", 1:n4)]
    m <- mu + rho * (Ysp[, -J4] - mu)
    ll[k, ] <- c(dnorm(Ysp[, 1], mu, s["sigma0"], log = TRUE), dnorm(Ysp[, -1], m, sqrt(s2), log = TRUE))
  }
  waic(ll)
}
post4s <- marathon_gibbs(common_rho = TRUE)
say("common-rho model: rho = %.4f (95%% CI %s)", mean(post4s[, "theta3"]), fmt_ci(quantile(post4s[, "theta3"], c(0.025, 0.975)), 4))
print(round(rbind(full = waic_marathon(post4), common_rho = waic_marathon(post4s)), 2))


# ---- Exercise 5: exponential survival for the lung data --------------------------------

section("Ex 5  lung cancer survival")

lung <- survival::lung
d5 <- as.numeric(lung$status == 2)                      # 1 = death observed, 0 = censored
t5 <- lung$time
say("n = %d   deaths = %d   censored = %d   total follow-up = %d days", nrow(lung), sum(d5), sum(1 - d5), sum(t5))
priors5 <- list("Gamma(0.01, 0.01)" = c(0.01, 0.01), "Gamma(0.1, 0.1)" = c(0.1, 0.1), "Gamma(1, 1)" = c(1, 1),
                "Gamma(10, 1000) informative" = c(10, 1000))
for (nm in names(priors5)) {
  a <- priors5[[nm]][1] + sum(d5); b <- priors5[[nm]][2] + sum(t5)
  say("%-28s -> Gamma(%.2f, %.2f): mean %.6f/day (median survival %.0f days)  95%% CI %s",
      nm, a, b, a / b, log(2) * b / a, fmt_ci(ci95(qgamma, a, b), 6))
}
sex <- factor(lung$sex, labels = c("male", "female"))
by_sex <- sapply(levels(sex), function(s) c(deaths = sum(d5[sex == s]), time = sum(t5[sex == s]), n = sum(sex == s)))
print(by_sex)
lam <- sapply(levels(sex), function(s) rgamma(1e6, 0.01 + by_sex["deaths", s], 0.01 + by_sex["time", s]))
say("rate male %.5f  female %.5f   P(lambda_male > lambda_female) = %.4f   hazard ratio M/F median %.3f  95%% CI %s",
    mean(lam[, "male"]), mean(lam[, "female"]), mean(lam[, "male"] > lam[, "female"]),
    median(lam[, "male"] / lam[, "female"]), fmt_ci(quantile(lam[, "male"] / lam[, "female"], c(0.025, 0.975))))
# Bayes factor: separate rates vs common rate, Gamma(a, b) priors (closed form)
log_m <- function(d, tt, a = 1, b = 0.001) a * log(b) - lgamma(a) + lgamma(a + d) - (a + d) * log(b + tt)
for (ab in list(c(1, 0.001), c(1, 100), c(1, 420), c(2, 840))) {
  lbf <- log_m(by_sex["deaths", "male"], by_sex["time", "male"], ab[1], ab[2]) +
         log_m(by_sex["deaths", "female"], by_sex["time", "female"], ab[1], ab[2]) - log_m(sum(d5), sum(t5), ab[1], ab[2])
  say("BF(separate : common) with Gamma(%g, %g) priors (prior mean rate %.4g/day) = %.3g", ab[1], ab[2], ab[1] / ab[2], exp(lbf))
}

save_fig("ex05_sex", width = 6.5, height = 3.8, function() {
  x <- seq(0.001, 0.008, length.out = 600)
  plot(x * 1000, dgamma(x, 0.01 + by_sex["deaths", "male"], 0.01 + by_sex["time", "male"]) / 1000, type = "l", lwd = 2.5,
       col = PALETTE["navy"], xlab = "death rate per 1,000 patient-days", ylab = "Posterior density",
       ylim = c(0, 2.1))
  lines(x * 1000, dgamma(x, 0.01 + by_sex["deaths", "female"], 0.01 + by_sex["time", "female"]) / 1000, lwd = 2.5,
        col = PALETTE["red"])
  legend("topright", c("men", "women"), lwd = 2.5, col = PALETTE[c("navy", "red")], bty = "n")
})


# ---- Exercise 6: spatial model for village malaria prevalence --------------------------

section("Ex 6  spatial malaria model")

gam <- read.csv(file.path("data", "gambia.csv"))
key <- paste(gam$x, gam$y)
vill <- data.frame(x = tapply(gam$x, key, `[`, 1) / 1000, y = tapply(gam$y, key, `[`, 1) / 1000,   # km
                   n = tapply(gam$pos, key, length), Y = tapply(gam$pos, key, sum))
write.csv(vill, file.path(CODE_DIR, "gambia_villages.csv"), row.names = FALSE)
V <- nrow(vill); Dmat <- as.matrix(dist(vill[, c("x", "y")])); dstar <- max(Dmat)
say("villages %d   children per village %d-%d   max distance %.1f km", V, min(vill$n), max(vill$n), dstar)

spatial_mcmc <- function(iters = 25000, burn = 5000, rho = dstar / 4) {
  alpha <- qlogis((vill$Y + 0.5) / (vill$n + 1)); mu <- mean(alpha); s2 <- var(alpha)
  chol_R <- function(rho) chol(exp(-Dmat / rho) + diag(1e-8, V))
  Rc <- chol_R(rho)
  loglik <- function(a) sum(dbinom(vill$Y, vill$n, plogis(a), log = TRUE))
  quad <- function(a, Rc) sum(backsolve(Rc, a - mu, transpose = TRUE)^2)
  out <- matrix(NA, iters, 3 + V, dimnames = list(NULL, c("mu", "sigma", "rho", paste0("a", 1:V))))
  acc_rho <- 0
  for (it in seq_len(iters)) {
    # alpha | mu, sigma, rho, Y : elliptical slice sampling (no tuning needed)
    nu <- mu + sqrt(s2) * drop(crossprod(Rc, rnorm(V)))
    logy <- loglik(alpha) + log(runif(1))
    ang <- runif(1, 0, 2 * pi); lo <- ang - 2 * pi; hi <- ang
    repeat {
      prop <- mu + (alpha - mu) * cos(ang) + (nu - mu) * sin(ang)
      if (loglik(prop) > logy) { alpha <- prop; break }
      if (ang < 0) lo <- ang else hi <- ang
      ang <- runif(1, lo, hi)
    }
    # mu | rest (normal) and sigma^2 | rest (inverse gamma)
    ones <- backsolve(Rc, rep(1, V), transpose = TRUE); za <- backsolve(Rc, alpha, transpose = TRUE)
    prec <- sum(ones^2) / s2 + 1 / 100
    mu <- rnorm(1, sum(ones * za) / s2 / prec, 1 / sqrt(prec))
    # rho | alpha, mu with sigma^2 integrated out (InvGamma conjugacy), Metropolis on logit(rho / d*)
    log_marg_rho <- function(Rc_, r) -sum(log(diag(Rc_))) - (0.1 + V / 2) * log(0.1 + quad(alpha, Rc_) / 2) +
      log(r) + log(dstar - r)                                    # + log-Jacobian of the logit transform
    eta_c <- qlogis(rho / dstar) + rnorm(1, 0, 0.6); cand <- dstar * plogis(eta_c)
    Rc_c <- chol_R(cand)
    if (log(runif(1)) < log_marg_rho(Rc_c, cand) - log_marg_rho(Rc, rho)) { rho <- cand; Rc <- Rc_c; acc_rho <- acc_rho + 1 }
    s2 <- 1 / rgamma(1, 0.1 + V / 2, 0.1 + quad(alpha, Rc) / 2)  # sigma^2 | alpha, mu, rho
    out[it, ] <- c(mu, sqrt(s2), rho, alpha)
  }
  structure(out[-seq_len(burn), ], acc_rho = acc_rho / iters)
}
chains6 <- lapply(c(5, 40, 120), function(r0) spatial_mcmc(rho = r0))
post6 <- do.call(rbind, chains6)
say("rho acceptance: %s", paste(sprintf("%.2f", sapply(chains6, attr, "acc_rho")), collapse = " "))
diag6 <- t(sapply(c("mu", "sigma", "rho"), function(p) {
  ch <- sapply(chains6, function(c_) c_[, p])
  c(post_summary(as.vector(ch), 3), Rhat = rhat(ch), ESS = sum(apply(ch, 2, ess)))
}))
print(round(diag6, 3))
pv <- plogis(post6[, paste0("a", 1:V)])
vill$p_mean <- colMeans(pv); vill$p_sd <- apply(pv, 2, sd); vill$p_gt_half <- colMeans(pv > 0.5)
p_rhat <- sapply(1:V, function(v) rhat(sapply(chains6, function(c_) plogis(c_[, paste0("a", v)]))))
p_ess  <- sapply(1:V, function(v) sum(sapply(chains6, function(c_) ess(plogis(c_[, paste0("a", v)])))))
say("village probabilities p_v: R-hat %.3f-%.3f   ESS min %.0f, median %.0f", min(p_rhat), max(p_rhat), min(p_ess), median(p_ess))
say("logit scale: max R-hat %.3f, min ESS %.0f (villages with no positive children have long left tails)",
    max(sapply(1:V, function(v) rhat(sapply(chains6, function(c_) c_[, paste0("a", v)])))),
    min(sapply(1:V, function(v) sum(sapply(chains6, function(c_) ess(c_[, paste0("a", v)]))))))
say("villages with P(p > 0.5) > 0.9: %d   correlation(sample prop, posterior mean) = %.3f",
    sum(vill$p_gt_half > 0.9), cor(vill$Y / vill$n, vill$p_mean))
say("shrinkage example: village with n = %d, sample proportion %.2f -> posterior mean %.2f",
    vill$n[which.min(vill$n)], (vill$Y / vill$n)[which.min(vill$n)], vill$p_mean[which.min(vill$n)])

save_fig("ex06_maps", width = 7.5, height = 6.6, function() {
  par(mfrow = c(3, 2), mar = c(3.6, 3.6, 2.2, 5.2), mgp = c(2.1, 0.5, 0))
  map_panel <- function(z, main, lim = range(z), digits = 2) {
    ramp <- colorRampPalette(c("#F4F6F7", PALETTE["gold"], PALETTE["red"]))(100)
    col_idx <- pmax(1, pmin(100, round(1 + 99 * (z - lim[1]) / diff(lim))))
    plot(vill$x, vill$y, pch = 21, bg = ramp[col_idx], cex = 1.7, col = "grey30", main = main,
         xlab = "easting (km)", ylab = "northing (km)", asp = 1)
    legend("topright", inset = c(-0.24, 0), xpd = TRUE, bty = "n", cex = 0.85, pch = 21, pt.cex = 1.4,
           pt.bg = ramp[c(100, 50, 1)], legend = formatC(c(lim[2], mean(lim), lim[1]), format = "f", digits = digits))
  }
  map_panel(vill$n, "sample size n_v", digits = 0)
  map_panel(vill$Y / vill$n, "sample proportion Y_v / n_v", c(0, 1))
  map_panel(vill$p_mean, "posterior mean of p_v", c(0, 1))
  map_panel(vill$p_sd, "posterior SD of p_v")
  map_panel(vill$p_gt_half, "P(p_v > 0.5 | data)", c(0, 1))
  plot(density(post6[, "rho"]), main = "posterior of rho (km)", xlab = expression(rho), lwd = 2,
       col = PALETTE["navy"], zero.line = FALSE)
})
say("villages with P(p > 0.5) > 0.9 by easting: %s km", paste(round(sort(vill$x[vill$p_gt_half > 0.9])), collapse = ", "))
