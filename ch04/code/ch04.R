# =============================================================================
#  Chapter 4 -- Linear models
#  Reproduces every number and figure in ch04/ch04.tex.
#  The JAGS version of Exercise 9 is in ch04/code/ch04_jags.R (run on GitHub).
#
#  Run from the repository root:   Rscript ch04/code/ch04.R
#  Needs base R + MASS (Boston, galaxies, glm.nb) + splines (ships with R).
#  Data files are in data/ (see data/README.md for sources).
# =============================================================================

source(file.path("common", "helpers.R"))
set.seed(2026)
FIG_DIR <- file.path("ch04", "figures")
load_data <- function(file) { e <- new.env(); load(file.path("data", file), envir = e); e }


# ---- Exercise 1: two-sample comparison of weight change -----------------------

section("Ex 1  placebo vs treatment")

placebo   <- c(2.0, -3.1, -1.0, 0.2, 0.3, 0.4)
treatment <- c(-3.5, -1.6, -4.6, -0.9, -5.1, 0.1)
n1 <- length(placebo); n2 <- length(treatment)

# (i) Equal variances, Jeffreys prior p(mu1, mu2, sigma^2) ∝ 1/sigma^2:
#     delta = mu2 - mu1 | data ~ t_{n1+n2-2}(ybar2 - ybar1, s_p^2 (1/n1 + 1/n2)).
sp2   <- ((n1 - 1) * var(placebo) + (n2 - 1) * var(treatment)) / (n1 + n2 - 2)
d_hat <- mean(treatment) - mean(placebo)
d_se  <- sqrt(sp2 * (1 / n1 + 1 / n2))
df    <- n1 + n2 - 2
say("means: placebo %.3f  treatment %.3f   diff %.3f   pooled sd %.3f   se %.3f",
    mean(placebo), mean(treatment), d_hat, sqrt(sp2), d_se)
say("equal variances: 95%% CI %s   P(delta < 0) = %.4f",
    fmt_ci(d_hat + qt(c(0.025, 0.975), df) * d_se), pt(-d_hat / d_se, df))

# (ii) Unequal variances: mu_g | data ~ t_{n_g - 1}(ybar_g, s_g^2 / n_g), independently.
S <- 1e6
mu_draw <- function(y) mean(y) + sd(y) / sqrt(length(y)) * rt(S, length(y) - 1)
delta_uneq <- mu_draw(treatment) - mu_draw(placebo)
say("unequal variances: mean %.3f  95%% CI %s   P(delta < 0) = %.4f",
    mean(delta_uneq), fmt_ci(quantile(delta_uneq, c(0.025, 0.975))), mean(delta_uneq < 0))

# (iii) Sensitivity: informative N(0, tau^2) prior on delta (sceptical), via Gibbs.
gibbs_two_sample <- function(tau, iters = 40000) {
  y <- c(placebo, treatment); g <- rep(0:1, c(n1, n2))
  mu <- mean(placebo); delta <- 0; s2 <- sp2
  out <- numeric(iters)
  for (it in seq_len(iters)) {
    prec <- length(y) / s2 + 1 / 100^2
    mu <- rnorm(1, sum(y - delta * g) / s2 / prec, 1 / sqrt(prec))
    prec <- n2 / s2 + 1 / tau^2
    delta <- rnorm(1, sum((y - mu)[g == 1]) / s2 / prec, 1 / sqrt(prec))
    s2 <- 1 / rgamma(1, 0.01 + length(y) / 2, 0.01 + sum((y - mu - delta * g)^2) / 2)
    out[it] <- delta
  }
  out[-(1:2000)]
}
sens <- t(sapply(c(100, 10, 2, 1, 0.5), function(tau) {
  d <- gibbs_two_sample(tau)
  c(tau = tau, mean = mean(d), lo = unname(quantile(d, 0.025)), hi = unname(quantile(d, 0.975)),
    p_neg = mean(d < 0))
}))
print(round(sens, 3))
say("classical two-sample t-test p-value (one-sided): %.4f",
    t.test(treatment, placebo, alternative = "less", var.equal = TRUE)$p.value)

save_fig("ex01_delta", width = 6.5, height = 3.8, function() {
  x <- seq(-7, 3, length.out = 800)
  plot(x, dt((x - d_hat) / d_se, df) / d_se, type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(delta == mu[T] - mu[P] ~ "(pounds)"), ylab = "Posterior density")
  lines(density(delta_uneq, from = -7, to = 3), lwd = 2, lty = 2, col = PALETTE["gold"])
  abline(v = 0, lty = 3, col = PALETTE["red"])
  legend("topleft", c("equal variances", "unequal variances"), lty = 1:2, lwd = 2,
         col = PALETTE[c("navy", "gold")], bty = "n")
})


# ---- Exercise 2: Boston housing ---------------------------------------------------

section("Ex 2  Boston housing regression")

boston <- MASS::Boston
covs <- setdiff(names(boston), "medv")
Xb <- cbind(intercept = 1, scale(boston[, covs]))           # standardised covariates
yb <- boston$medv

fit2 <- gibbs_lm(yb, Xb, prior_sd = 100, iters = 12000, burn = 2000)
ols  <- lm(medv ~ ., data = data.frame(medv = yb, Xb[, -1]))
tab2 <- cbind(summary_table(fit2[, 1:ncol(Xb)]),
              ols_est = round(coef(ols), 3), ols_se = round(summary(ols)$coef[, 2], 3),
              ESS = round(apply(fit2[, 1:ncol(Xb)], 2, ess)))
print(tab2)
say("sigma: Bayes %.3f   OLS %.3f", mean(fit2[, "sigma"]), summary(ols)$sigma)

# (c) Bayesian LASSO (Park & Casella 2008), lambda^2 ~ Gamma(1, 0.1).
bayes_lasso <- function(y, X, iters = 12000, burn = 2000) {
  n <- nrow(X); p <- ncol(X)
  yc <- y - mean(y); XtX <- crossprod(X); Xty <- crossprod(X, yc)
  beta <- rep(0, p); s2 <- var(yc); inv_tau2 <- rep(1, p); lambda2 <- 1
  draws <- matrix(NA, iters, p + 3, dimnames = list(NULL, c("intercept", colnames(X), "sigma", "lambda")))
  for (it in seq_len(iters)) {
    A    <- XtX + diag(inv_tau2, p)
    R    <- chol(A)
    m    <- backsolve(R, forwardsolve(t(R), Xty))
    beta <- drop(m + sqrt(s2) * backsolve(R, rnorm(p)))
    sse  <- sum((yc - X %*% beta)^2)
    s2   <- 1 / rgamma(1, (n - 1 + p) / 2, (sse + sum(beta^2 * inv_tau2)) / 2)
    inv_tau2 <- rinvgauss(p, sqrt(lambda2 * s2 / beta^2), lambda2)
    lambda2  <- rgamma(1, p + 1, 0.1 + sum(1 / inv_tau2) / 2)
    alpha <- rnorm(1, mean(y), sqrt(s2 / n))
    draws[it, ] <- c(alpha, beta, sqrt(s2), sqrt(lambda2))
  }
  draws[-seq_len(burn), ]
}
fit2c <- bayes_lasso(yb, Xb[, -1])
tab2c <- cbind(flat = round(colMeans(fit2[, 1:ncol(Xb)]), 3), summary_table(fit2c[, 1:ncol(Xb)]))
print(tab2c)
say("lambda posterior mean %.3f", mean(fit2c[, "lambda"]))

save_fig("ex02_coefficients", width = 7, height = 4.6, function() {
  k <- ncol(Xb) - 1; at <- seq_len(k)
  par(mar = c(4.2, 6, 1, 1))
  q_flat  <- apply(fit2[, 2:ncol(Xb)], 2, quantile, c(0.025, 0.5, 0.975))
  q_lasso <- apply(fit2c[, 2:ncol(Xb)], 2, quantile, c(0.025, 0.5, 0.975))
  plot(NA, xlim = range(q_flat, q_lasso), ylim = c(0.5, k + 0.5), yaxt = "n",
       xlab = "coefficient (change in medv per SD of covariate)", ylab = "")
  axis(2, at = at, labels = covs, las = 1)
  abline(v = 0, lty = 3, col = PALETTE["grey"])
  segments(q_flat[1, ], at + 0.15, q_flat[3, ], at + 0.15, lwd = 3, col = PALETTE["navy"])
  points(q_flat[2, ], at + 0.15, pch = 19, col = PALETTE["navy"])
  segments(q_lasso[1, ], at - 0.15, q_lasso[3, ], at - 0.15, lwd = 3, col = PALETTE["red"])
  points(q_lasso[2, ], at - 0.15, pch = 19, col = PALETTE["red"])
  legend("bottomright", c("Gaussian prior", "double-exponential prior"), lwd = 3,
         col = PALETTE[c("navy", "red")], bty = "n")
})

# (d) fit on the first 500 rows, predict the last 6.
train <- 1:500; test <- 501:506
fit2d <- gibbs_lm(yb[train], Xb[train, ], prior_sd = 100, iters = 12000, burn = 2000)
beta_d <- fit2d[, 1:ncol(Xb)]
ppd <- beta_d %*% t(Xb[test, ]) + fit2d[, "sigma"] * matrix(rnorm(nrow(fit2d) * 6), ncol = 6)
ppd_tab <- data.frame(row = test, actual = yb[test], mean = colMeans(ppd),
                      lo = apply(ppd, 2, quantile, 0.025), hi = apply(ppd, 2, quantile, 0.975))
print(round(ppd_tab, 2))

save_fig("ex02_ppd", width = 7.5, height = 3.4, function() {
  par(mfrow = c(1, 6), mar = c(4, 2, 2, 0.5))
  for (j in 1:6) {
    d <- density(ppd[, j])
    plot(d, main = paste("obs", test[j]), xlab = "medv", ylab = "", yaxt = "n", lwd = 2,
         col = PALETTE["navy"], zero.line = FALSE, xlim = range(d$x, yb[test[j]]))
    abline(v = yb[test[j]], col = PALETTE["red"], lwd = 2)
  }
})


# ---- Exercise 3: 2016 election ------------------------------------------------------

section("Ex 3  2016 election swing")

el <- load_data("election_2008_2016.RData")
keep <- !el$junk
y3 <- 100 * (el$pct_gop_16 - el$pct_gop_12)[keep]
X3raw <- el$X[keep, ]
X3 <- cbind(intercept = 1, scale(X3raw))
county <- paste0(el$all_dat$area_name, ", ", el$all_dat$state_abbreviation)[keep]
state  <- el$all_dat$state_abbreviation[keep]
labels3 <- el$names$description[match(colnames(X3raw), el$names$column_name)]
say("counties used: %d   mean swing %.2f points   sd %.2f", length(y3), mean(y3), sd(y3))

fit3 <- gibbs_lm(y3, X3, prior_sd = 100, iters = 6000, burn = 1000)
tab3 <- summary_table(fit3)
print(cbind(tab3, ESS = round(apply(fit3, 2, ess))))

beta3 <- colMeans(fit3[, 1:ncol(X3)])
resid3 <- drop(y3 - X3 %*% beta3)
say("residual sd %.3f  skewness %.3f  excess kurtosis %.3f", sd(resid3),
    mean((resid3 - mean(resid3))^3) / sd(resid3)^3, mean((resid3 - mean(resid3))^4) / sd(resid3)^4 - 3)
ord <- order(resid3)
cat("Most negative residuals:\n"); print(data.frame(county = county[head(ord, 5)], resid = round(resid3[head(ord, 5)], 2)))
cat("Most positive residuals:\n"); print(data.frame(county = county[tail(ord, 5)], resid = round(resid3[tail(ord, 5)], 2)))
by_state <- sort(tapply(resid3, state, mean))
say("state-average residual range: %s %.2f  ...  %s %.2f",
    names(by_state)[1], by_state[1], names(by_state)[length(by_state)], tail(by_state, 1))

save_fig("ex03_residuals", width = 7.5, height = 3.4, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.2, 1.2, 0.8))
  hist(resid3, breaks = 60, freq = FALSE, col = adjustcolor(PALETTE["navy"], 0.6), border = "white",
       main = "", xlab = "residual (percentage points)")
  curve(dnorm(x, 0, sd(resid3)), add = TRUE, lwd = 2, col = PALETTE["red"])
  qqnorm(resid3, pch = 16, cex = 0.4, col = PALETTE["navy"], main = "")
  qqline(resid3, col = PALETTE["red"], lwd = 2)
})

# (c) state random effects, hierarchically centred for good mixing:
#     Y_i = X_i beta + alpha_{s_i} + e_i,  alpha_l ~ N(mu, tau^2)   (mu plays the intercept's role,
#     so X has no intercept column).  Reported "intercept" = mu.
gibbs_re <- function(y, X, group, iters = 6000, burn = 1000) {
  g <- as.integer(factor(group)); L <- max(g); n_l <- tabulate(g, L)
  n <- length(y); p <- ncol(X)
  XtX <- crossprod(X)
  beta <- qr.solve(X, y - mean(y)); mu <- mean(y); alpha <- rep(mu, L); s2 <- var(y); t2 <- 1
  draws <- matrix(NA, iters, p + 3 + L,
                  dimnames = list(NULL, c("intercept", colnames(X), "sigma", "tau", levels(factor(group)))))
  for (it in seq_len(iters)) {
    r    <- y - alpha[g]                                        # beta | alpha
    prec <- XtX / s2 + diag(1 / 100^2, p)
    R    <- chol(prec)
    beta <- drop(backsolve(R, forwardsolve(t(R), crossprod(X, r) / s2)) + backsolve(R, rnorm(p)))
    e    <- drop(y - X %*% beta)
    prec_a <- n_l / s2 + 1 / t2                                 # alpha_l | beta, mu
    alpha  <- rnorm(L, (rowsum(e, g)[, 1] / s2 + mu / t2) / prec_a, 1 / sqrt(prec_a))
    prec_m <- L / t2 + 1 / 100^2                                # mu | alpha
    mu <- rnorm(1, sum(alpha) / t2 / prec_m, 1 / sqrt(prec_m))
    s2 <- 1 / rgamma(1, 0.01 + n / 2, 0.01 + sum((e - alpha[g])^2) / 2)
    t2 <- 1 / rgamma(1, 0.01 + L / 2, 0.01 + sum((alpha - mu)^2) / 2)
    draws[it, ] <- c(mu, beta, sqrt(s2), sqrt(t2), alpha - mu)  # report effects as deviations
  }
  draws[-seq_len(burn), ]
}
fit3c <- gibbs_re(y3, X3[, -1], state)
p3 <- ncol(X3)
tab3c <- cbind(no_RE_mean = round(colMeans(fit3[, 1:p3]), 3), no_RE_sd = round(apply(fit3[, 1:p3], 2, sd), 3),
               RE_mean = round(colMeans(fit3c[, 1:p3]), 3), RE_sd = round(apply(fit3c[, 1:p3], 2, sd), 3))
print(tab3c)
say("sigma: no RE %.3f  with RE %.3f   tau %.3f", mean(fit3[, "sigma"]), mean(fit3c[, "sigma"]), mean(fit3c[, "tau"]))
alpha_mean <- sort(colMeans(fit3c[, -(1:(p3 + 2))]))
cat("Lowest state effects:\n");  print(round(head(alpha_mean, 5), 2))
cat("Highest state effects:\n"); print(round(tail(alpha_mean, 5), 2))
say("ESS (of 5000): tau %.0f   state effects min %.0f median %.0f", ess(fit3c[, "tau"]),
    min(apply(fit3c[, -(1:(p3 + 2))], 2, ess)), median(apply(fit3c[, -(1:(p3 + 2))], 2, ess)))
say("ESS (of 5000) for coefficients: min %.0f", min(apply(fit3c[, 1:p3], 2, ess)))

save_fig("ex03_state_effects", width = 7.5, height = 4, function() {
  a <- fit3c[, -(1:(p3 + 2))]
  o <- order(colMeans(a)); q <- apply(a[, o], 2, quantile, c(0.025, 0.5, 0.975))
  par(mar = c(3, 4.2, 1, 0.5))
  plot(seq_along(o), q[2, ], ylim = range(q), pch = 19, cex = 0.6, col = PALETTE["navy"], xaxt = "n",
       xlab = "", ylab = "state effect (points)")
  segments(seq_along(o), q[1, ], seq_along(o), q[3, ], col = PALETTE["navy"])
  axis(1, at = seq_along(o), labels = colnames(a)[o], las = 2, cex.axis = 0.55)
  abline(h = 0, lty = 3, col = PALETTE["grey"])
})


# ---- Exercise 4: gun laws and homicides -------------------------------------------

section("Ex 4  gun laws and homicides")

gun <- load_data("guns.Rdata")
laws <- rowSums(gun$X)
Z <- gun$Z
colnames(Z) <- c("rate2009", "ownership", "unemployment", "nonfirearm", "export")
X4 <- cbind(intercept = 1, scale(Z), laws = (laws - mean(laws)) / sd(laws))
y4 <- gun$Y; off4 <- log(gun$N)
say("laws: mean %.2f  sd %.2f  range %d-%d", mean(laws), sd(laws), min(laws), max(laws))

log_post_pois <- function(b) {
  eta <- off4 + drop(X4 %*% b)
  sum(y4 * eta - exp(eta)) + sum(dnorm(b, 0, 10, log = TRUE))
}
glm4 <- glm(y4 ~ X4 - 1 + offset(off4), family = poisson)
init4 <- setNames(coef(glm4), colnames(X4))
chains4 <- lapply(1:3, function(k) metropolis_block(log_post_pois, init4 + rnorm(7, 0, 0.1),
                                                    (2.4^2 / 7) * vcov(glm4), iters = 25000, burn = 5000))
say("Poisson acceptance: %s", paste(sprintf("%.2f", sapply(chains4, attr, "acceptance")), collapse = " "))
post4 <- do.call(rbind, chains4)
rh4 <- sapply(colnames(X4), function(p) rhat(sapply(chains4, function(ch) ch[, p])))
print(cbind(summary_table(post4, 4), Rhat = round(rh4, 3),
            ESS = round(rowSums(sapply(chains4, function(ch) apply(ch, 2, ess))))))
say("laws effect: rate ratio per extra law = %.4f (95%% CI %s)",
    mean(exp(post4[, "laws"] / sd(laws))), fmt_ci(exp(quantile(post4[, "laws"] / sd(laws), c(0.025, 0.975))), 4))
say("Pearson dispersion (Poisson glm): %.2f", sum(residuals(glm4, "pearson")^2) / glm4$df.residual)

# (b) negative binomial: Y ~ NB(mean N*lambda, size r), log r sampled too.
nb4 <- MASS::glm.nb(y4 ~ X4 - 1 + offset(off4))
log_post_nb <- function(th) {
  b <- th[1:7]; r <- exp(th[8])
  mu <- exp(off4 + drop(X4 %*% b))
  sum(dnbinom(y4, size = r, mu = mu, log = TRUE)) + sum(dnorm(b, 0, 10, log = TRUE)) +
    dgamma(r, 0.1, 0.1, log = TRUE) + th[8]
}
init_nb <- c(setNames(coef(nb4), colnames(X4)), log_r = log(nb4$theta))
cov_nb <- matrix(0, 8, 8)
cov_nb[1:7, 1:7] <- vcov(nb4)
cov_nb[8, 8] <- (nb4$SE.theta / nb4$theta)^2             # delta-method variance of log(r)
cov_nb <- (2.4^2 / 8) * cov_nb
chains4b <- lapply(1:3, function(k) metropolis_block(log_post_nb, init_nb + rnorm(8, 0, 0.05), cov_nb,
                                                     iters = 25000, burn = 5000))
post4b <- do.call(rbind, chains4b)
say("NB acceptance: %s", paste(sprintf("%.2f", sapply(chains4b, attr, "acceptance")), collapse = " "))
rh4b <- sapply(colnames(post4b), function(p) rhat(sapply(chains4b, function(ch) ch[, p])))
print(cbind(poisson_mean = c(round(colMeans(post4), 4), NA), poisson_sd = c(round(apply(post4, 2, sd), 4), NA),
            summary_table(post4b, 4), Rhat = round(rh4b, 3)))
say("NB size r: posterior mean %.1f", mean(exp(post4b[, "log_r"])))

# (c) posterior predictive national totals with laws = 0 and laws = 25 (Poisson model).
ppd_total <- function(n_laws, post) {
  Xs <- X4; Xs[, "laws"] <- (n_laws - mean(laws)) / sd(laws)
  lam <- exp(sweep(post %*% t(Xs), 2, off4, "+"))
  rowSums(matrix(rpois(length(lam), lam), nrow(lam)))
}
idx <- sample(nrow(post4), 20000)
tot_obs <- sum(y4)
tot0  <- ppd_total(0, post4[idx, ]); tot25 <- ppd_total(25, post4[idx, ])
tot_now <- { lam <- exp(sweep(post4[idx, ] %*% t(X4), 2, off4, "+")); rowSums(matrix(rpois(length(lam), lam), nrow(lam))) }
say("observed total %d   PPD current laws mean %.0f", tot_obs, mean(tot_now))
say("laws = 0 : mean %.0f  95%% PI %s", mean(tot0), fmt_ci(quantile(tot0, c(0.025, 0.975)), 0))
say("laws = 25: mean %.0f  95%% PI %s", mean(tot25), fmt_ci(quantile(tot25, c(0.025, 0.975)), 0))

save_fig("ex04_ppd", width = 6.5, height = 3.8, function() {
  r <- range(tot0, tot25, tot_now)
  hist(tot0, breaks = 60, freq = FALSE, col = adjustcolor(PALETTE["red"], 0.5), border = NA, xlim = r,
       main = "", xlab = "nationwide homicides (posterior predictive)", ylab = "Density")
  hist(tot25, breaks = 60, freq = FALSE, col = adjustcolor(PALETTE["green"], 0.5), border = NA, add = TRUE)
  hist(tot_now, breaks = 60, freq = FALSE, col = adjustcolor(PALETTE["grey"], 0.5), border = NA, add = TRUE)
  abline(v = tot_obs, lwd = 2, lty = 2)
  legend("topright", c("0 laws in every state", "current laws", "25 laws in every state", "observed"),
         fill = c(adjustcolor(PALETTE[c("red", "grey", "green")], 0.5), NA), border = NA,
         lty = c(NA, NA, NA, 2), bty = "n")
})


# ---- Exercise 5: Titanic logistic regression ------------------------------------

section("Ex 5  Titanic survival")

tit <- read.csv(file.path("data", "titanic_train.csv"))
tit <- tit[!is.na(tit$Age), ]
X5 <- cbind(intercept = 1, age10 = (tit$Age - 30) / 10, female = as.numeric(tit$Sex == "female"),
            class2 = as.numeric(tit$Pclass == 2), class3 = as.numeric(tit$Pclass == 3))
y5 <- tit$Survived
say("complete cases: %d of 891   survival rate %.3f", nrow(tit), mean(y5))
glm5 <- glm(y5 ~ X5 - 1, family = binomial)
log_post_logit <- function(b, X, y, prior_sd = 10) {
  eta <- drop(X %*% b)
  sum(y * eta - log1p(exp(eta))) + sum(dnorm(b, 0, prior_sd, log = TRUE))
}
chains5 <- lapply(1:3, function(k) metropolis_block(function(b) log_post_logit(b, X5, y5),
                                                    setNames(coef(glm5), colnames(X5)) + rnorm(5, 0, 0.2),
                                                    (2.4^2 / 5) * vcov(glm5), iters = 20000, burn = 4000))
post5 <- do.call(rbind, chains5)
say("acceptance %.2f", attr(chains5[[1]], "acceptance"))
tab5 <- summary_table(post5)
odds5 <- t(apply(exp(post5), 2, function(x) round(c(OR = median(x), quantile(x, c(0.025, 0.975))), 3)))
rh5 <- sapply(colnames(X5), function(p) rhat(sapply(chains5, function(ch) ch[, p])))
print(cbind(tab5, odds5, Rhat = round(rh5, 3), glm = round(coef(glm5), 3)))
# Survival probabilities for a 30-year-old in each class, by sex.
grid5 <- expand.grid(female = 0:1, class = 1:3)
prob5 <- apply(grid5, 1, function(r) {
  x <- c(1, 0, r["female"], r["class"] == 2, r["class"] == 3)
  quantile(plogis(post5 %*% x), c(0.5, 0.025, 0.975))
})
print(cbind(grid5, round(t(prob5), 3)))


# ---- Exercise 6: gamma regression for T. rex mass -------------------------------

section("Ex 6  gamma regression for T. rex")

mass <- c(29.9, 1761, 1807, 2984, 3230, 5040, 5654)
age  <- c(2, 15, 14, 16, 18, 22, 28)
glm6 <- glm(mass ~ age, family = Gamma(link = "log"))
log_post6 <- function(th) {
  mu <- exp(th[1] + th[2] * age); nu <- exp(th[3])
  sum(dgamma(mass, shape = nu, rate = nu / mu, log = TRUE)) +
    sum(dnorm(th[1:2], 0, 10, log = TRUE)) + dgamma(nu, 0.1, 0.1, log = TRUE) + th[3]
}
nu0 <- 1 / summary(glm6)$dispersion
init6 <- c(beta0 = coef(glm6)[[1]], beta1 = coef(glm6)[[2]], log_nu = log(nu0))
cov6 <- diag(c(0.3, 0.015, 0.5)^2)
chains6 <- lapply(1:3, function(k) metropolis_block(log_post6, init6 + rnorm(3, 0, 0.1), cov6,
                                                    iters = 60000, burn = 10000))
post6 <- do.call(rbind, chains6)
post6 <- cbind(post6, nu = exp(post6[, "log_nu"]), growth_pct = 100 * (exp(post6[, "beta1"]) - 1))
say("acceptance %s", paste(sprintf("%.2f", sapply(chains6, attr, "acceptance")), collapse = " "))
rh6 <- sapply(c("beta0", "beta1", "log_nu"), function(p) rhat(sapply(chains6, function(ch) ch[, p])))
print(summary_table(post6, 4)); print(round(rh6, 4))
say("glm: beta0 %.3f beta1 %.4f shape %.2f", coef(glm6)[1], coef(glm6)[2], nu0)

save_fig("ex06_fit", width = 6, height = 3.8, function() {
  ages <- seq(0, 30, length.out = 200)
  curves <- exp(outer(post6[, "beta0"], rep(1, 200)) + outer(post6[, "beta1"], ages))
  band <- apply(curves, 2, quantile, c(0.025, 0.5, 0.975))
  plot(age, mass, pch = 19, cex = 1.2, col = PALETTE["navy"], xlim = c(0, 30), ylim = c(0, 12000),
       xlab = "age (years)", ylab = "mass (kg)")
  polygon(c(ages, rev(ages)), c(band[1, ], rev(band[3, ])), col = adjustcolor(PALETTE["gold"], 0.25), border = NA)
  lines(ages, colMeans(curves), lwd = 2.5, col = PALETTE["red"])
})


# ---- Exercise 8: Gambia malaria, logistic and random-effects logistic ------------

section("Ex 8  Gambia malaria")

gam <- read.csv(file.path("data", "gambia.csv"))
loc <- as.integer(factor(paste(gam$x, gam$y)))
L8 <- max(loc)
X8 <- cbind(intercept = 1, age = (gam$age / 365 - mean(gam$age / 365)) / sd(gam$age / 365),
            netuse = gam$netuse, treated = gam$treated,
            green = (gam$green - mean(gam$green)) / sd(gam$green), phc = gam$phc,
            x = (gam$x - mean(gam$x)) / sd(gam$x), y = (gam$y - mean(gam$y)) / sd(gam$y))
y8 <- gam$pos
say("children %d  locations %d  prevalence %.3f", nrow(gam), L8, mean(y8))

glm8 <- glm(y8 ~ X8 - 1, family = binomial)
chains8 <- lapply(1:3, function(k) metropolis_block(function(b) log_post_logit(b, X8, y8),
                                                    setNames(coef(glm8), colnames(X8)) + rnorm(8, 0, 0.1),
                                                    (2.4^2 / 8) * vcov(glm8), iters = 15000, burn = 3000))
post8 <- do.call(rbind, chains8)
rh8 <- sapply(colnames(X8), function(p) rhat(sapply(chains8, function(ch) ch[, p])))
say("acceptance %.2f", attr(chains8[[1]], "acceptance"))
print(cbind(summary_table(post8), Rhat = round(rh8, 3)))

# (b) random location effects, hierarchically centred: alpha_l ~ N(mu, tau^2), X without intercept.
#     Metropolis for beta (block) and for all 65 alpha_l at once; Gibbs for mu and tau^2.
X8r <- X8[, -1]
re_logit <- function(iters = 25000, burn = 5000, init = coef(glm8)[-1]) {
  pb <- ncol(X8r)
  beta <- init; mu <- coef(glm8)[1]; alpha <- rep(mu, L8); t2 <- 0.5
  prop_b <- t(chol((2.4^2 / pb) * vcov(glm8)[-1, -1]))    # conditional on alpha, spread is close to the GLM's
  eta_fixed <- drop(X8r %*% beta)
  loglik_obs <- function(eta) y8 * eta - log1p(exp(eta))
  out <- matrix(NA, iters, pb + 2 + L8, dimnames = list(NULL, c("intercept", colnames(X8r), "tau", paste0("a", 1:L8))))
  acc <- c(beta = 0, alpha = 0)
  for (it in seq_len(iters)) {
    cand <- beta + drop(prop_b %*% rnorm(pb))
    eta_c <- drop(X8r %*% cand)
    lr <- sum(loglik_obs(eta_c + alpha[loc])) - sum(loglik_obs(eta_fixed + alpha[loc])) +
      sum(dnorm(cand, 0, 10, log = TRUE)) - sum(dnorm(beta, 0, 10, log = TRUE))
    if (log(runif(1)) < lr) { beta <- cand; eta_fixed <- eta_c; acc["beta"] <- acc["beta"] + 1 }
    a_c <- alpha + rnorm(L8, 0, 0.45)
    ll_now <- rowsum(loglik_obs(eta_fixed + alpha[loc]), loc)[, 1]
    ll_c   <- rowsum(loglik_obs(eta_fixed + a_c[loc]), loc)[, 1]
    lr_a <- ll_c - ll_now + dnorm(a_c, mu, sqrt(t2), log = TRUE) - dnorm(alpha, mu, sqrt(t2), log = TRUE)
    take <- log(runif(L8)) < lr_a
    alpha[take] <- a_c[take]; acc["alpha"] <- acc["alpha"] + mean(take)
    prec_m <- L8 / t2 + 1 / 100
    mu <- rnorm(1, sum(alpha) / t2 / prec_m, 1 / sqrt(prec_m))
    t2 <- 1 / rgamma(1, 0.1 + L8 / 2, 0.1 + sum((alpha - mu)^2) / 2)
    out[it, ] <- c(mu, beta, sqrt(t2), alpha - mu)
  }
  structure(out[-seq_len(burn), ], acceptance = acc / iters)
}
chains8b <- lapply(1:3, function(k) re_logit(init = coef(glm8)[-1] + rnorm(7, 0, 0.2)))
post8b <- do.call(rbind, chains8b)
say("RE acceptance: beta %.2f  alpha %.2f", attr(chains8b[[1]], "acceptance")[1], attr(chains8b[[1]], "acceptance")[2])
show8 <- c(colnames(X8), "tau")
rh8b <- sapply(show8, function(p) rhat(sapply(chains8b, function(ch) ch[, p])))
ess8b <- sapply(show8, function(p) sum(sapply(chains8b, function(ch) ess(ch[, p]))))
print(cbind(no_RE = c(round(colMeans(post8), 3), tau = NA),
            summary_table(post8b[, show8]), Rhat = round(rh8b, 3), ESS = round(ess8b)))
alpha8 <- colMeans(post8b[, paste0("a", 1:L8)])
coords <- unique(data.frame(loc = loc, x = gam$x, y = gam$y))
coords <- coords[order(coords$loc), ]

save_fig("ex08_map", width = 6.5, height = 3.6, function() {
  par(mar = c(4.2, 4.4, 1, 6))
  ramp <- colorRampPalette(c(PALETTE["navy"], "#F4F6F7", PALETTE["red"]))(101)
  lim <- max(abs(alpha8))
  col_idx <- round(50 + 50 * alpha8 / lim) + 1
  plot(coords$x / 1000, coords$y / 1000, pch = 21, bg = ramp[col_idx], cex = 1.8, col = "grey30",
       xlab = "easting (km)", ylab = "northing (km)", asp = 1)
  legend("topright", inset = c(-0.27, 0), xpd = TRUE, bty = "n", title = "village effect",
         legend = sprintf("%+.1f", c(lim, lim / 2, 0, -lim / 2, -lim)),
         pt.bg = ramp[c(101, 76, 51, 26, 1)], pch = 21, pt.cex = 1.6)
})
say("village effects range %.2f to %.2f", min(alpha8), max(alpha8))


# ---- Exercise 9: AR(1) model for the popularity of "Sophia" -------------------

section("Ex 9  AR(1) model for Sophia")

soph <- read.csv(file.path("data", "sophia.csv"))
Y9 <- log(soph$prop / (1 - soph$prop)); t9 <- soph$year - 1950; T9 <- length(Y9)

rtrunc_norm <- function(m, s, lo, hi) qnorm(runif(1, pnorm(lo, m, s), pnorm(hi, m, s)), m, s)

gibbs_ar1 <- function(iters = 30000, burn = 5000, rho = 0) {
  y <- Y9[-1]; ylag <- Y9[-T9]; tt <- t9[-1]; tlag <- t9[-T9]
  ab <- c(0, 0); s2 <- 1
  out <- matrix(NA, iters, 4, dimnames = list(NULL, c("alpha", "beta", "rho", "sigma")))
  for (it in seq_len(iters)) {
    # (alpha, beta) | rho, sigma: linear regression of y - rho*ylag on (1 - rho, t - rho*tlag)
    Z  <- cbind(1 - rho, tt - rho * tlag); w <- y - rho * ylag
    prec <- crossprod(Z) / s2 + diag(1 / 100^2, 2)
    V <- solve(prec); ab <- drop(V %*% crossprod(Z, w) / s2 + t(chol(V)) %*% rnorm(2))
    # rho | rest: regression of e_t on e_{t-1}, truncated to (-1, 1)
    e <- y - ab[1] - ab[2] * tt; elag <- ylag - ab[1] - ab[2] * tlag
    rho <- rtrunc_norm(sum(e * elag) / sum(elag^2), sqrt(s2 / sum(elag^2)), -1, 1)
    res <- e - rho * elag
    s2 <- 1 / rgamma(1, 0.1 + length(y) / 2, 0.1 + sum(res^2) / 2)
    out[it, ] <- c(ab, rho, sqrt(s2))
  }
  out[-seq_len(burn), ]
}
chains9 <- lapply(c(-0.5, 0, 0.9), function(r0) gibbs_ar1(rho = r0))
post9 <- do.call(rbind, chains9)
rh9 <- sapply(colnames(post9), function(p) rhat(sapply(chains9, function(ch) ch[, p])))
ess9 <- sapply(colnames(post9), function(p) sum(sapply(chains9, function(ch) ess(ch[, p]))))
print(cbind(summary_table(post9, 4), Rhat = round(rh9, 4), ESS = round(ess9)))

# (c) PPD for 2020 (t = 70): simulate 2018, 2019, 2020 forward from 2017.
idx9 <- sample(nrow(post9), 30000)
pp <- post9[idx9, ]
y_prev <- rep(Y9[T9], nrow(pp)); t_prev <- t9[T9]
for (tt in (t_prev + 1):70) {
  mu_t <- pp[, "alpha"] + pp[, "beta"] * tt; mu_prev <- pp[, "alpha"] + pp[, "beta"] * (tt - 1)
  y_prev <- mu_t + pp[, "rho"] * (y_prev - mu_prev) + rnorm(nrow(pp), 0, pp[, "sigma"])
}
y2020 <- y_prev
say("last observed (2017): log-odds %.3f  proportion %.5f", Y9[T9], soph$prop[T9])
say("PPD 2020 log-odds: mean %.3f  95%% PI %s", mean(y2020), fmt_ci(quantile(y2020, c(0.025, 0.975))))
say("PPD 2020 proportion: median %.5f  95%% PI %s", median(plogis(y2020)),
    fmt_ci(quantile(plogis(y2020), c(0.025, 0.975)), 5))

save_fig("ex09_series", width = 7.5, height = 3.4, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.4, 1.2, 0.6))
  plot(soph$year, Y9, type = "b", pch = 19, cex = 0.5, col = PALETTE["navy"],
       xlab = "year", ylab = "log-odds of the name Sophia", xlim = c(1951, 2020))
  plot(density(y2020), main = "", lwd = 2.5, col = PALETTE["navy"], zero.line = FALSE,
       xlab = "log-odds in 2020 (posterior predictive)")
  abline(v = Y9[T9], lty = 2, col = PALETTE["red"])
})
write.csv(data.frame(t = t9, Y = Y9), file.path("ch04", "code", "ex09_sophia_logodds.csv"), row.names = FALSE)


# ---- Exercise 10: B-spline regression for global temperature ------------------

section("Ex 10  B-spline temperature curve")

tmp <- load_data("global_mean_temp.RData")
yr <- tmp$year; anom <- tmp$anomaly
B <- splines::bs(yr, df = 15, intercept = TRUE)
colnames(B) <- paste0("B", 1:15)
fit10 <- gibbs_lm(anom, B, prior_sd = 100, iters = 12000, burn = 2000)
g_draws <- fit10[, 1:15] %*% t(B)                       # posterior draws of g(X_i)
g_band <- apply(g_draws, 2, quantile, c(0.025, 0.5, 0.975))
delta_trend <- g_draws[, yr == 2020] - g_draws[, yr == 1980]
say("years %d-%d   sigma %.3f", min(yr), max(yr), mean(fit10[, "sigma"]))
say("Delta = g(2020) - g(1980): mean %.3f  95%% CI %s  P(Delta > 0) = %.4f",
    mean(delta_trend), fmt_ci(quantile(delta_trend, c(0.025, 0.975))), mean(delta_trend > 0))
d_i <- g_draws[, -1] - g_draws[, -ncol(g_draws)]
d_band <- apply(d_i, 2, quantile, c(0.025, 0.5, 0.975))
d_pos  <- colMeans(d_i > 0)
say("years with P(delta_i > 0) > 0.95: %d of %d; first such year %d",
    sum(d_pos > 0.95), length(d_pos), yr[-1][which(d_pos > 0.95)[1]])
say("largest posterior mean warming rate %.4f C/yr in %d", max(colMeans(d_i)), yr[-1][which.max(colMeans(d_i))])
say("posterior mean rate in 2021: %.4f C/yr (95%% CI %s)", mean(d_i[, ncol(d_i)]),
    fmt_ci(quantile(d_i[, ncol(d_i)], c(0.025, 0.975)), 4))

save_fig("ex10_curve", width = 7, height = 3.8, function() {
  plot(yr, anom, pch = 19, cex = 0.5, col = PALETTE["grey"], xlab = "year", ylab = "temperature anomaly (C)")
  polygon(c(yr, rev(yr)), c(g_band[1, ], rev(g_band[3, ])), col = adjustcolor(PALETTE["red"], 0.25), border = NA)
  lines(yr, colMeans(g_draws), lwd = 2.5, col = PALETTE["red"])
})
save_fig("ex10_rate", width = 7, height = 3.6, function() {
  yy <- yr[-1]
  plot(yy, colMeans(d_i), type = "n", ylim = range(d_band), xlab = "year",
       ylab = expression(delta[i] == g(X[i]) - g(X[i - 1]) ~ "(C / year)"))
  polygon(c(yy, rev(yy)), c(d_band[1, ], rev(d_band[3, ])), col = adjustcolor(PALETTE["navy"], 0.25), border = NA)
  lines(yy, colMeans(d_i), lwd = 2.5, col = PALETTE["navy"])
  abline(h = 0, lty = 2, col = PALETTE["red"])
})


# ---- Exercise 13: three-component normal mixture for galaxies -----------------

section("Ex 13  normal mixture for galaxies")

gal <- MASS::galaxies / 1000                              # work in 1,000 km/s
K <- 3; n13 <- length(gal)
gibbs_mix <- function(iters = 20000, burn = 5000) {
  mu <- quantile(gal, c(0.1, 0.5, 0.9)); s2 <- rep(4, K); w <- rep(1 / K, K)
  grid <- seq(5, 40, by = 0.1)
  dens <- matrix(NA, iters - burn, length(grid))
  pars <- matrix(NA, iters - burn, 3 * K)
  for (it in seq_len(iters)) {
    lp <- sapply(1:K, function(k) log(w[k]) + dnorm(gal, mu[k], sqrt(s2[k]), log = TRUE))
    pr <- exp(lp - apply(lp, 1, max)); pr <- pr / rowSums(pr)
    u  <- runif(n13)
    z  <- 1 + (u > pr[, 1]) + (u > pr[, 1] + pr[, 2])             # draw each label from its row of pr
    nk <- tabulate(z, K)
    w  <- { g <- rgamma(K, 1 + nk); g / sum(g) }                  # Dirichlet(1 + n_k)
    for (k in 1:K) {
      yk <- gal[z == k]
      prec <- nk[k] / s2[k] + 1 / 100
      mu[k] <- rnorm(1, (sum(yk) / s2[k] + 20 / 100) / prec, 1 / sqrt(prec))   # prior N(20, 10^2)
      s2[k] <- 1 / rgamma(1, 1 + nk[k] / 2, 1 + sum((yk - mu[k])^2) / 2)      # prior InvGamma(1, 1)
    }
    if (it > burn) {
      dens[it - burn, ] <- rowSums(sapply(1:K, function(k) w[k] * dnorm(grid, mu[k], sqrt(s2[k])))) / 1000
      o <- order(mu); pars[it - burn, ] <- c(mu[o], sqrt(s2[o]), w[o])
    }
  }
  list(grid = grid * 1000, dens = dens, pars = pars)
}
mix <- gibbs_mix()
colnames(mix$pars) <- c(paste0("mu", 1:3), paste0("sd", 1:3), paste0("w", 1:3))
print(summary_table(mix$pars, 3))
band13 <- apply(mix$dens, 2, quantile, c(0.025, 0.5, 0.975))
say("density matrix: %d draws x %d grid points", nrow(mix$dens), ncol(mix$dens))

save_fig("ex13_mixture", width = 7, height = 3.8, function() {
  hist(MASS::galaxies, breaks = 25, freq = FALSE, col = adjustcolor(PALETTE["navy"], 0.4), border = "white",
       main = "", xlab = "velocity (km/s)", ylab = "Density", ylim = c(0, max(band13) * 1.05))
  polygon(c(mix$grid, rev(mix$grid)), c(band13[1, ], rev(band13[3, ])),
          col = adjustcolor(PALETTE["red"], 0.25), border = NA)
  lines(mix$grid, band13[2, ], lwd = 2.5, col = PALETTE["red"])
})
