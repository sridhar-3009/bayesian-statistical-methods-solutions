# =============================================================================
#  Chapter 6 -- Model selection and diagnostics
#  Reproduces every number and figure in ch06/ch06.tex.
#
#  Run from the repository root:   Rscript ch06/code/ch06.R
#  Needs base R + MASS (galaxies) + splines (ships with R).
# =============================================================================

source(file.path("common", "helpers.R"))
set.seed(2026)
FIG_DIR <- file.path("ch06", "figures")
skewness <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurtosis <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3      # excess kurtosis
log_mean_exp <- function(l) max(l) + log(mean(exp(l - max(l))))


# ---- Exercises 1-2: airquality, cross-validation and posterior predictive checks ------

section("Ex 1  airquality: 5-fold cross-validation")

aq <- na.omit(airquality[, c("Ozone", "Solar.R", "Temp", "Wind")])
n1 <- nrow(aq)
design <- function(d, vars, centre, scale_) cbind(intercept = 1, sweep(sweep(as.matrix(d[, vars, drop = FALSE]), 2, centre[vars]), 2, scale_[vars], "/"))
centre1 <- colMeans(aq); scale1 <- apply(aq, 2, sd)
models1 <- list(M1 = "Solar.R", M2 = c("Solar.R", "Temp", "Wind"))
say("complete cases: %d   mean ozone %.1f   sd %.1f", n1, mean(aq$Ozone), sd(aq$Ozone))

folds <- sample(rep(1:5, length.out = n1))
cv1 <- sapply(models1, function(vars) {
  sq_err <- numeric(n1); lpd <- numeric(n1)
  for (f in 1:5) {
    tr <- folds != f; te <- folds == f
    Xtr <- design(aq[tr, ], vars, centre1, scale1); Xte <- design(aq[te, ], vars, centre1, scale1)
    fit <- gibbs_lm(aq$Ozone[tr], Xtr, prior_sd = 100, a = 0.1, b = 0.1, iters = 4000, burn = 1000)
    mu <- fit[, 1:ncol(Xtr)] %*% t(Xte)                          # S x n_test predictive means
    sq_err[te] <- (aq$Ozone[te] - colMeans(mu))^2
    lpd[te] <- sapply(seq_len(sum(te)), function(i)
      log_mean_exp(dnorm(aq$Ozone[te][i], mu[, i], fit[, "sigma"], log = TRUE)))
  }
  c(RMSE = sqrt(mean(sq_err)), MAE = mean(sqrt(sq_err)), log_score = sum(lpd))
})
print(round(cv1, 3))

fit_full <- lapply(models1, function(vars)
  gibbs_lm(aq$Ozone, design(aq, vars, centre1, scale1), prior_sd = 100, a = 0.1, b = 0.1, iters = 11000, burn = 1000))
for (m in names(fit_full)) { cat(m, "\n"); print(summary_table(fit_full[[m]])) }

section("Ex 2  posterior predictive checks for M2")

X2 <- design(aq, models1$M2, centre1, scale1)
f2 <- fit_full$M2
S2 <- nrow(f2)
mu2 <- f2[, 1:4] %*% t(X2)
yrep <- mu2 + f2[, "sigma"] * matrix(rnorm(S2 * n1), S2)
stats_fun <- list(mean = mean, sd = sd, min = min, max = max, skewness = skewness,
                  `share < 0` = function(y) mean(y < 0))
ppc2 <- t(sapply(names(stats_fun), function(s) {
  obs <- stats_fun[[s]](aq$Ozone); rep_s <- apply(yrep, 1, stats_fun[[s]])
  c(observed = obs, rep_mean = mean(rep_s), rep_lo = unname(quantile(rep_s, 0.025)),
    rep_hi = unname(quantile(rep_s, 0.975)), ppp = ppp(rep_s, obs))
}))
print(round(ppc2, 3))
fitted2 <- colMeans(mu2); resid2 <- aq$Ozone - fitted2
say("cor(fitted, resid^2) = %.3f   share of fitted values < 0: %.3f", cor(fitted2, resid2^2), mean(fitted2 < 0))

save_fig("ex02_ppc", width = 7.5, height = 5, function() {
  par(mfrow = c(2, 3), mar = c(4, 3.5, 2, 0.5))
  for (s in c("sd", "max", "skewness", "share < 0")) {
    rep_s <- apply(yrep, 1, stats_fun[[s]]); obs <- stats_fun[[s]](aq$Ozone)
    hist(rep_s, breaks = 40, col = adjustcolor(PALETTE["navy"], 0.5), border = "white", main = s,
         xlab = "replicated", ylab = "", xlim = range(rep_s, obs))
    abline(v = obs, lwd = 2.5, col = PALETTE["red"])
  }
  plot(fitted2, resid2, pch = 19, cex = 0.6, col = PALETTE["navy"], xlab = "fitted", ylab = "residual",
       main = "residuals vs fitted")
  abline(h = 0, lty = 2, col = PALETTE["grey"])
  plot(aq$Temp, aq$Ozone, pch = 19, cex = 0.6, col = PALETTE["navy"], xlab = "Temp", ylab = "Ozone",
       main = "ozone vs temperature")
  lines(lowess(aq$Temp, aq$Ozone), lwd = 2, col = PALETTE["red"])
})


# ---- Exercise 3: DIC and WAIC for Mr. October ----------------------------------------

section("Ex 3  DIC and WAIC, one rate or two")

Y3 <- c(563, 10); N3 <- c(2820, 27)
rtrunc_gamma <- function(n, shape, rate, upper) qgamma(runif(n, 0, pgamma(upper, shape, rate)), shape, rate)
S3 <- 50000
ic3 <- do.call(rbind, lapply(c(1, 10), function(c3) {
  lam_sep <- sapply(1:2, function(j) rtrunc_gamma(S3, Y3[j] + 1, N3[j], c3))       # M1
  lam_com <- rtrunc_gamma(S3, sum(Y3) + 1, sum(N3), c3)                              # M2
  ll1 <- cbind(dpois(Y3[1], N3[1] * lam_sep[, 1], log = TRUE), dpois(Y3[2], N3[2] * lam_sep[, 2], log = TRUE))
  ll2 <- cbind(dpois(Y3[1], N3[1] * lam_com, log = TRUE), dpois(Y3[2], N3[2] * lam_com, log = TRUE))
  d1 <- dic(-2 * rowSums(ll1), -2 * sum(dpois(Y3, N3 * colMeans(lam_sep), log = TRUE)))
  d2 <- dic(-2 * rowSums(ll2), -2 * sum(dpois(Y3, N3 * mean(lam_com), log = TRUE)))
  rbind(c(c = c3, model = 1, d1, waic(ll1)), c(c = c3, model = 2, d2, waic(ll2)))
}))
print(round(ic3, 3))


# ---- Exercise 4: collinear genetic markers (small simulation for intuition) -----------

section("Ex 4  collinear markers, p = 3")

n4 <- 1000; rho4 <- 0.99
Sigma4 <- rho4^abs(outer(1:3, 1:3, "-"))
X4 <- matrix(rnorm(n4 * 3), n4) %*% chol(Sigma4)
y4 <- drop(X4 %*% c(0.3, 0.3, 0.3)) + rnorm(n4)              # all three markers matter
ols4 <- summary(lm(y4 ~ X4))
print(round(ols4$coefficients, 3))
say("cor(beta1_hat, beta2_hat) = %.3f   sum of slopes %.3f (se %.3f)",
    cov2cor(vcov(lm(y4 ~ X4)))[2, 3], sum(coef(lm(y4 ~ X4))[-1]),
    sqrt(sum(vcov(lm(y4 ~ X4))[-1, -1])))
# exact posterior over the 8 subsets, same prior as Exercise 5.9 (W = 1, gamma_j ~ Bernoulli(0.5))
yc4 <- y4 - mean(y4); Xc4 <- scale(X4, scale = FALSE)
lm4 <- function(g) {
  k <- sum(g); yy <- sum(yc4^2)
  if (k == 0) return(-(n4 - 1) / 2 * log(yy))
  Xg <- Xc4[, g, drop = FALSE]; A <- crossprod(Xg) + diag(k); R <- chol(A)
  v <- forwardsolve(t(R), crossprod(Xg, yc4))
  -sum(log(diag(R))) - (n4 - 1) / 2 * log(yy - sum(v^2))
}
mods4 <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 3)))
l4 <- apply(mods4, 1, lm4); pm4 <- exp(l4 - max(l4)); pm4 <- pm4 / sum(pm4)
print(data.frame(model = apply(mods4, 1, function(g) if (!any(g)) "none" else paste0("X", which(g), collapse = "+")),
                 prob = round(pm4, 3)), row.names = FALSE)
say("inclusion probabilities: %s", paste(sprintf("%.3f", colSums(mods4 * pm4)), collapse = "  "))


# ---- Exercises 5-6: Gambia, logit vs probit, and PPCs ---------------------------------

section("Ex 5  Gambia: logit vs probit")

gam <- read.csv(file.path("data", "gambia.csv"))
vil <- as.integer(factor(paste(gam$x, gam$y)))
X5 <- cbind(intercept = 1, age = as.numeric(scale(gam$age / 365)), netuse = gam$netuse,
            treated = gam$treated, green = as.numeric(scale(gam$green)), phc = gam$phc)
y5 <- gam$pos
link_fun <- list(logit = plogis, probit = pnorm)
fit5 <- lapply(names(link_fun), function(lk) {
  F <- link_fun[[lk]]
  log_post <- function(b) {
    p <- F(drop(X5 %*% b))
    sum(dbinom(y5, 1, p, log = TRUE)) + sum(dnorm(b, 0, 10, log = TRUE))
  }
  g <- glm(y5 ~ X5 - 1, family = binomial(link = lk))
  draws <- metropolis_block(log_post, setNames(coef(g), colnames(X5)), (2.4^2 / 6) * vcov(g),
                            iters = 22000, burn = 2000)
  P <- F(draws %*% t(X5))                                          # S x n success probabilities
  ll <- sweep(log(P), 2, y5, "*") + sweep(log(1 - P), 2, 1 - y5, "*")   # S x n pointwise log-lik
  dev_mean <- -2 * sum(dbinom(y5, 1, F(drop(X5 %*% colMeans(draws))), log = TRUE))
  list(draws = draws, P = P, ic = c(dic(-2 * rowSums(ll), dev_mean), waic(ll)),
       acc = attr(draws, "acceptance"))
})
names(fit5) <- names(link_fun)
print(round(t(sapply(fit5, `[[`, "ic")), 2))
print(round(cbind(logit = colMeans(fit5$logit$draws), probit = colMeans(fit5$probit$draws),
                  ratio = colMeans(fit5$logit$draws) / colMeans(fit5$probit$draws)), 3))
say("acceptance: logit %.2f  probit %.2f", fit5$logit$acc, fit5$probit$acc)

section("Ex 6  Gambia: PPC for the logistic model")

P6 <- fit5$logit$P; S6 <- nrow(P6)
idx6 <- sample(S6, 4000)
yrep6 <- matrix(rbinom(length(idx6) * ncol(P6), 1, P6[idx6, ]), length(idx6))
vil_prev <- function(y) tapply(y, vil, mean)
stat6 <- list(
  overall = function(y) mean(y),
  `sd of village prevalence` = function(y) sd(vil_prev(y)),
  `max village prevalence` = function(y) max(vil_prev(y)),
  `villages with 0 positives` = function(y) sum(vil_prev(y) == 0))
ppc6 <- t(sapply(names(stat6), function(s) {
  obs <- stat6[[s]](y5); rep_s <- apply(yrep6, 1, stat6[[s]])
  c(observed = obs, rep_mean = mean(rep_s), rep_lo = unname(quantile(rep_s, 0.025)),
    rep_hi = unname(quantile(rep_s, 0.975)), ppp = ppp(rep_s, obs))
}))
print(round(ppc6, 3))
age_grp <- cut(gam$age / 365, c(0, 1.5, 2.5, 3.5, 5))
obs_age <- tapply(y5, age_grp, mean)
rep_age <- apply(yrep6, 1, function(y) tapply(y, age_grp, mean))
print(round(cbind(observed = obs_age, rep_lo = apply(rep_age, 1, quantile, 0.025),
                  rep_hi = apply(rep_age, 1, quantile, 0.975)), 3))

save_fig("ex06_ppc", width = 7.5, height = 3, function() {
  par(mfrow = c(1, 3), mar = c(4, 3.5, 2, 0.5))
  for (s in c("overall", "sd of village prevalence", "villages with 0 positives")) {
    rep_s <- apply(yrep6, 1, stat6[[s]]); obs <- stat6[[s]](y5)
    hist(rep_s, breaks = 30, col = adjustcolor(PALETTE["navy"], 0.5), border = "white", main = s,
         xlab = "replicated", ylab = "", xlim = range(rep_s, obs))
    abline(v = obs, lwd = 2.5, col = PALETTE["red"])
  }
})


# ---- Exercises 7-8: NBA clutch shooting, LOO cross-validation and PPC ----------------

section("Ex 7  NBA: leave-one-out cross-validation")

nba <- data.frame(
  player = c("Westbrook", "Harden", "Leonard", "James", "Thomas", "Curry",
             "Antetokounmpo", "Wall", "Davis", "Durant"),
  q = c(0.845, 0.847, 0.880, 0.674, 0.909, 0.898, 0.770, 0.801, 0.802, 0.875),
  y = c(64, 72, 55, 27, 75, 24, 28, 66, 40, 13),
  n = c(75, 95, 63, 39, 83, 26, 41, 82, 54, 16))
nba$lq <- qlogis(nba$q)

fit_nba <- function(d, model, iters = 12000) {
  lp <- function(b) {
    eta <- b[1] + (if (model == 1) 1 else b[2]) * d$lq
    sum(dbinom(d$y, d$n, plogis(eta), log = TRUE)) + sum(dnorm(b, 0, 10, log = TRUE))
  }
  g <- if (model == 1) glm(cbind(y, n - y) ~ 1 + offset(lq), family = binomial, data = d)
       else glm(cbind(y, n - y) ~ lq, family = binomial, data = d)
  init <- setNames(coef(g), if (model == 1) "b1" else c("b1", "b2"))
  metropolis_block(lp, init, (2.4^2 / length(init)) * vcov(g), iters = iters, burn = 2000)
}
loo7 <- sapply(1:2, function(model) {
  sapply(seq_len(nrow(nba)), function(i) {
    dr <- fit_nba(nba[-i, ], model)
    eta <- dr[, 1] + (if (model == 1) 1 else dr[, 2]) * nba$lq[i]
    log_mean_exp(dbinom(nba$y[i], nba$n[i], plogis(eta), log = TRUE))
  })
})
colnames(loo7) <- c("M1", "M2")
print(cbind(nba[, c("player", "y", "n")], round(loo7, 3)))
say("LOO log score: M1 %.3f   M2 %.3f   difference (M2 - M1) %.3f  (se %.3f)",
    sum(loo7[, 1]), sum(loo7[, 2]), sum(loo7[, 2] - loo7[, 1]),
    sqrt(nrow(nba) * var(loo7[, 2] - loo7[, 1])))
full7 <- lapply(1:2, function(m) fit_nba(nba, m, iters = 22000))
print(summary_table(full7[[1]])); print(summary_table(full7[[2]]))

section("Ex 8  NBA: posterior predictive checks for M2")

d8 <- full7[[2]]
P8 <- plogis(d8[, 1] + d8[, 2] %o% nba$lq)                            # S x 10
yrep8 <- matrix(rbinom(length(P8), rep(nba$n, each = nrow(P8)), P8), nrow(P8))
chi2 <- function(y, P) rowSums((matrix(y, nrow(P), length(y), byrow = TRUE) - sweep(P, 2, nba$n, "*"))^2 /
                                 sweep(P * (1 - P), 2, nba$n, "*"))
chi_obs <- chi2(nba$y, P8); chi_rep <- sapply(seq_len(nrow(P8)), function(s) chi2(yrep8[s, ], P8[s, , drop = FALSE]))
say("chi-square discrepancy: observed mean %.2f   PPP = P(rep >= obs) = %.3f", mean(chi_obs), mean(chi_rep >= chi_obs))
ppp_player <- sapply(1:10, function(i) mean(yrep8[, i] <= nba$y[i]))
print(data.frame(player = nba$player, observed = nba$y, rep_mean = round(colMeans(yrep8), 1),
                 rep_lo = apply(yrep8, 2, quantile, 0.025), rep_hi = apply(yrep8, 2, quantile, 0.975),
                 `P(rep <= obs)` = round(ppp_player, 3), check.names = FALSE), row.names = FALSE)

save_fig("ex08_ppc", width = 7.5, height = 3.6, function() {
  par(mar = c(6.5, 4.2, 1, 0.5))
  q <- apply(yrep8 / matrix(nba$n, nrow(yrep8), 10, byrow = TRUE), 2, quantile, c(0.025, 0.5, 0.975))
  plot(1:10, nba$y / nba$n, pch = 19, col = PALETTE["red"], ylim = range(q, nba$y / nba$n), xaxt = "n",
       xlab = "", ylab = "clutch make proportion")
  segments(1:10, q[1, ], 1:10, q[3, ], lwd = 3, col = adjustcolor(PALETTE["navy"], 0.6))
  points(1:10, q[2, ], pch = 3, col = PALETTE["navy"])
  axis(1, at = 1:10, labels = nba$player, las = 2, cex.axis = 0.8)
  legend("bottomleft", c("observed", "95% posterior predictive interval"), pch = c(19, NA),
         lwd = c(NA, 3), col = c(PALETTE["red"], PALETTE["navy"]), bty = "n")
})


# ---- Exercises 9-10: WWWusage autoregressions ------------------------------------------

section("Ex 9  WWWusage: choosing the lag")

www <- as.numeric(WWWusage)
t9 <- 5:100
y9 <- www[t9]
fits9 <- lapply(1:4, function(L) {
  X <- cbind(intercept = 1, sapply(1:L, function(l) www[t9 - l]))
  colnames(X)[-1] <- paste0("lag", 1:L)
  dr <- gibbs_lm(y9, X, prior_sd = 100, a = 0.1, b = 0.1, iters = 11000, burn = 1000)
  mu <- dr[, 1:ncol(X)] %*% t(X)
  ll <- dnorm(matrix(y9, nrow(mu), length(y9), byrow = TRUE), mu, dr[, "sigma"], log = TRUE)
  dev_mean <- -2 * sum(dnorm(y9, drop(X %*% colMeans(dr[, 1:ncol(X)])), mean(dr[, "sigma"]), log = TRUE))
  list(draws = dr, X = X, ic = c(dic(-2 * rowSums(ll), dev_mean), waic(ll)))
})
ic9 <- t(sapply(fits9, `[[`, "ic")); rownames(ic9) <- paste0("L = ", 1:4)
print(round(ic9, 2))
for (L in 1:4) { cat("L =", L, "\n"); print(summary_table(fits9[[L]]$draws, 3)) }

section("Ex 10  WWWusage: posterior predictive checks for L = 2")

f10 <- fits9[[2]]; X10 <- f10$X; d10 <- f10$draws
mu10 <- d10[, 1:3] %*% t(X10)
res_obs <- sweep(-mu10, 2, y9, "+")                                 # residual draws
acf1 <- function(e) { e <- e - mean(e); sum(e[-1] * e[-length(e)]) / sum(e^2) }
yrep10 <- mu10 + d10[, "sigma"] * matrix(rnorm(length(mu10)), nrow(mu10))
res_rep <- yrep10 - mu10
stat10 <- function(e) c(acf1 = acf1(e), acf2 = { e <- e - mean(e); sum(e[-(1:2)] * e[-(length(e) - 0:1)]) / sum(e^2) },
                       kurtosis = kurtosis(e), max_abs = max(abs(e)))
obs10 <- t(apply(res_obs, 1, stat10)); rep10 <- t(apply(res_rep, 1, stat10))
print(round(rbind(observed_mean = colMeans(obs10), replicated_mean = colMeans(rep10),
                  PPP = colMeans(rep10 >= obs10)), 3))
say("sample ACF of differenced series at lag 1: %.3f", acf(diff(www), plot = FALSE)$acf[2])

save_fig("ex10_ppc", width = 7.5, height = 3.2, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.2, 1.5, 0.5))
  plot(1:100, www, type = "l", lwd = 2, col = PALETTE["navy"], xlab = "minute", ylab = "users connected")
  lines(t9, colMeans(mu10), lwd = 1.5, lty = 2, col = PALETTE["red"])
  plot(rep10[, "acf1"], obs10[, "acf1"], pch = 16, cex = 0.3, col = adjustcolor(PALETTE["navy"], 0.3),
       xlab = "lag-1 autocorrelation, replicated", ylab = "lag-1 autocorrelation, observed")
  abline(0, 1, lty = 2, col = PALETTE["red"])
})


# ---- Exercise 11: a prior on c that makes p_D uniform ----------------------------------

section("Ex 11  uniform prior on p_D")

n11 <- 5; p11 <- 10
u <- runif(1e5); c_draw <- u / (n11 * (1 - u))                       # c = u / (n (1 - u))
pD <- p11 * c_draw / (c_draw + 1 / n11)
say("simulated p_D: mean %.3f (theory %.1f)   KS distance to Uniform(0, %d): %.4f",
    mean(pD), p11 / 2, p11, max(abs(ecdf(pD)(seq(0, p11, length.out = 1001)) - seq(0, 1, length.out = 1001))))
say("median of the prior on c = %.3f (= 1/n)", median(c_draw))


# ---- Exercise 12: choosing the number of B-spline basis functions ---------------------

section("Ex 12  temperature curve: choosing J")

tmp <- new.env(); load(file.path("data", "global_mean_temp.RData"), envir = tmp)
yr12 <- tmp$year; y12 <- tmp$anomaly; n12 <- length(y12)
J_grid <- seq(4, 30, by = 2)
basis <- function(J) { B <- splines::bs(yr12, df = J, intercept = TRUE); colnames(B) <- paste0("B", 1:J); B }

ic12 <- t(sapply(J_grid, function(J) {
  B <- basis(J)
  dr <- gibbs_lm(y12, B, prior_sd = 100, a = 0.01, b = 0.001, iters = 5000, burn = 1000)
  mu <- dr[, 1:J] %*% t(B)
  ll <- dnorm(matrix(y12, nrow(mu), n12, byrow = TRUE), mu, dr[, "sigma"], log = TRUE)
  dev_mean <- -2 * sum(dnorm(y12, drop(B %*% colMeans(dr[, 1:J])), mean(dr[, "sigma"]), log = TRUE))
  c(J = J, dic(-2 * rowSums(ll), dev_mean), waic(ll))
}))
folds12 <- sample(rep(1:10, length.out = n12))
cv12 <- sapply(J_grid, function(J) {
  B <- basis(J); sq <- numeric(n12)
  for (f in 1:10) {
    tr <- folds12 != f
    dr <- gibbs_lm(y12[tr], B[tr, ], prior_sd = 100, a = 0.01, b = 0.001, iters = 2500, burn = 500)
    sq[!tr] <- (y12[!tr] - drop(B[!tr, , drop = FALSE] %*% colMeans(dr[, 1:J])))^2
  }
  sqrt(mean(sq))
})
res12 <- cbind(ic12[, c("J", "DIC", "pD", "WAIC", "pW")], CV_RMSE = cv12)
print(round(res12, 4))
J_dic <- J_grid[which.min(ic12[, "DIC"])]; J_waic <- J_grid[which.min(ic12[, "WAIC"])]; J_cv <- J_grid[which.min(cv12)]
say("best J: DIC %d   WAIC %d   10-fold CV %d", J_dic, J_waic, J_cv)

save_fig("ex12_selection", width = 7.5, height = 3.2, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.4, 1.2, 0.6))
  plot(J_grid, ic12[, "DIC"], type = "b", pch = 19, col = PALETTE["navy"], xlab = "J", ylab = "criterion",
       ylim = range(ic12[, c("DIC", "WAIC")]))
  lines(J_grid, ic12[, "WAIC"], type = "b", pch = 17, col = PALETTE["red"])
  legend("topright", c("DIC", "WAIC"), pch = c(19, 17), col = PALETTE[c("navy", "red")], bty = "n")
  plot(J_grid, cv12, type = "b", pch = 19, col = PALETTE["green"], xlab = "J", ylab = "10-fold CV RMSE (C)")
})

# (c) posterior predictive checks for the selected J
J_sel <- J_waic
B12 <- basis(J_sel)
d12 <- gibbs_lm(y12, B12, prior_sd = 100, a = 0.01, b = 0.001, iters = 11000, burn = 1000)
mu12 <- d12[, 1:J_sel] %*% t(B12)
yrep12 <- mu12 + d12[, "sigma"] * matrix(rnorm(length(mu12)), nrow(mu12))
res_o <- sweep(-mu12, 2, y12, "+"); res_r <- yrep12 - mu12
st12 <- function(e) c(acf1 = acf1(e), kurtosis = kurtosis(e), max_abs = max(abs(e)),
                      runs = sum(diff(sign(e)) != 0))
o12 <- t(apply(res_o, 1, st12)); r12 <- t(apply(res_r, 1, st12))
print(round(rbind(observed_mean = colMeans(o12), replicated_mean = colMeans(r12),
                  PPP = colMeans(r12 >= o12)), 3))


# ---- Exercise 13: t vs Gaussian for the galaxies --------------------------------------

section("Ex 13  galaxies: PPC for t and Gaussian models")

gal <- MASS::galaxies; n13 <- length(gal)
t_chains <- gibbs_t(gal, iters = 22000)[-(1:2000), ]
gibbs_normal <- function(y, iters = 22000, burn = 2000) {
  n <- length(y); mu <- mean(y); tau <- 1 / var(y)
  out <- matrix(NA, iters, 2, dimnames = list(NULL, c("mu", "sigma")))
  for (it in seq_len(iters)) {
    prec <- n * tau + 1 / 10000^2
    mu   <- rnorm(1, tau * sum(y) / prec, 1 / sqrt(prec))
    tau  <- rgamma(1, 0.01 + n / 2, 0.01 + sum((y - mu)^2) / 2)
    out[it, ] <- c(mu, 1 / sqrt(tau))
  }
  out[-seq_len(burn), ]
}
n_chains <- gibbs_normal(gal)
idx13 <- sample(nrow(t_chains), 5000)
rep_t <- t(sapply(idx13, function(s) t_chains[s, "mu"] + t_chains[s, "sigma"] * rt(n13, t_chains[s, "k"])))
rep_n <- t(sapply(idx13, function(s) rnorm(n13, n_chains[s, "mu"], n_chains[s, "sigma"])))
st13 <- list(mean = mean, variance = var, skewness = skewness, kurtosis = kurtosis)
ppc13 <- do.call(rbind, lapply(names(st13), function(s) {
  obs <- st13[[s]](gal)
  rt_ <- apply(rep_t, 1, st13[[s]]); rn_ <- apply(rep_n, 1, st13[[s]])
  data.frame(statistic = s, observed = signif(obs, 4),
             t_median = signif(median(rt_), 4), t_lo = signif(quantile(rt_, 0.025), 4),
             t_hi = signif(quantile(rt_, 0.975), 4), t_ppp = round(ppp(rt_, obs), 3),
             normal_median = signif(median(rn_), 4), normal_lo = signif(quantile(rn_, 0.025), 4),
             normal_hi = signif(quantile(rn_, 0.975), 4), normal_ppp = round(ppp(rn_, obs), 3))
}))
print(ppc13, row.names = FALSE)
print(summary_table(n_chains, 1))

save_fig("ex13_ppc", width = 7.5, height = 5, function() {
  par(mfrow = c(2, 4), mar = c(4, 3, 2.2, 0.5))
  for (model in c("t", "Gaussian")) for (s in names(st13)) {
    rs <- apply(if (model == "t") rep_t else rep_n, 1, st13[[s]]); obs <- st13[[s]](gal)
    if (s == "variance") { rs <- log10(rs); obs <- log10(obs) }        # heavy tail under the t model
    lim <- range(quantile(rs, c(0.01, 0.95)), obs)                       # trim the far tail for display
    hist(rs[rs >= lim[1] & rs <= lim[2]], breaks = 40, col = adjustcolor(PALETTE["navy"], 0.5), border = "white",
         main = paste(model, "-", if (s == "variance") "log10 variance" else s), xlab = "", ylab = "", xlim = lim)
    abline(v = obs, lwd = 2.5, col = PALETTE["red"])
  }
})
