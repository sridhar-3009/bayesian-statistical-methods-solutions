# =============================================================================
#  Chapter 8 -- Machine learning
#  Reproduces every number and figure in ch08/ch08.tex.
#
#  Run from the repository root:   Rscript ch08/code/ch08.R
#  Needs base R + MASS (mcycle) + splines + BART (Exercises 7-10) + causaldata (Exercise 10).
# =============================================================================

source(file.path("common", "helpers.R"))
set.seed(2026)
FIG_DIR <- file.path("ch08", "figures")
quiet <- function(expr) { invisible(capture.output(res <- expr)); res }      # silence BART's progress output
log_mean_exp <- function(l) max(l) + log(mean(exp(l - max(l))))


# ---- Exercise 1: shrinkage factor under the horseshoe --------------------------------

section("Ex 1  horseshoe shrinkage factor")

tau <- abs(rcauchy(1e6))                       # half-Cauchy(1)
kappa <- 1 / (1 + tau^2)
say("kappa ~ Beta(1/2, 1/2)?  KS distance %.4f   mean %.4f (theory 0.5)   P(kappa < 0.1) = %.4f (theory %.4f)",
    max(abs(ecdf(kappa)(seq(0, 1, by = 0.001)) - pbeta(seq(0, 1, by = 0.001), 0.5, 0.5))),
    mean(kappa), mean(kappa < 0.1), pbeta(0.1, 0.5, 0.5))
save_fig("ex01_kappa", width = 6, height = 3.6, function() {
  hist(kappa, breaks = 60, freq = FALSE, col = adjustcolor(PALETTE["navy"], 0.5), border = "white", main = "",
       xlab = expression(kappa == 1 / (1 + tau^2)), ylab = "Density", ylim = c(0, 5))
  curve(dbeta(x, 0.5, 0.5), add = TRUE, lwd = 2.5, col = PALETTE["red"], n = 500)
})


# ---- Exercise 2: horseshoe logistic regression for the Gambia data --------------------

section("Ex 2  horseshoe logistic regression")

gam <- read.csv(file.path("data", "gambia.csv"))
Xg <- scale(cbind(age = gam$age, netuse = gam$netuse, treated = gam$treated, green = gam$green, phc = gam$phc))
yg <- gam$pos; pg <- ncol(Xg)
X1 <- cbind(intercept = 1, Xg)
rinvgamma <- function(n, a, b) 1 / rgamma(n, a, b)

# Horseshoe: beta_j ~ N(0, lambda_j^2 tau^2), lambda_j, tau ~ half-Cauchy(1), via the
# inverse-gamma auxiliary representation (Makalic & Schmidt, 2016); beta by block Metropolis.
horseshoe_logit <- function(iters = 25000, burn = 5000) {
  g <- glm(yg ~ Xg, family = binomial)
  b <- coef(g); lam2 <- rep(1, pg); nu <- rep(1, pg); tau2 <- 1; xi <- 1
  L <- t(chol((2.4^2 / (pg + 1)) * vcov(g)))
  log_post <- function(b, prior_var) {
    eta <- drop(X1 %*% b)
    sum(yg * eta - log1p(exp(eta))) + dnorm(b[1], 0, 10, log = TRUE) + sum(dnorm(b[-1], 0, sqrt(prior_var), log = TRUE))
  }
  out <- matrix(NA, iters, 2 * pg + 2, dimnames = list(NULL, c(colnames(X1), paste0("kappa_", colnames(Xg)), "tau")))
  acc <- 0
  for (it in seq_len(iters)) {
    pv <- lam2 * tau2
    cand <- b + drop(L %*% rnorm(pg + 1))
    if (log(runif(1)) < log_post(cand, pv) - log_post(b, pv)) { b <- cand; acc <- acc + 1 }
    bj <- b[-1]
    lam2 <- rinvgamma(pg, 1, 1 / nu + bj^2 / (2 * tau2))
    nu   <- rinvgamma(pg, 1, 1 + 1 / lam2)
    tau2 <- rinvgamma(1, (pg + 1) / 2, 1 / xi + sum(bj^2 / lam2) / 2)
    xi   <- rinvgamma(1, 1, 1 + 1 / tau2)
    out[it, ] <- c(b, 1 / (1 + lam2 * tau2), sqrt(tau2))
  }
  structure(out[-seq_len(burn), ], acceptance = acc / iters)
}
chains2 <- lapply(1:3, function(k) horseshoe_logit())
post2 <- do.call(rbind, chains2)
say("acceptance %.2f", attr(chains2[[1]], "acceptance"))
cols2 <- colnames(X1)
tab2 <- t(sapply(cols2, function(p) {
  ch <- sapply(chains2, function(c_) c_[, p])
  c(median = median(ch), lo = unname(quantile(ch, 0.025)), hi = unname(quantile(ch, 0.975)),
    Rhat = rhat(ch), ESS = sum(apply(ch, 2, ess)))
}))
print(round(tab2, 3))
say("posterior mean shrinkage kappa_j: %s", paste(sprintf("%s %.3f", colnames(Xg), colMeans(post2[, paste0("kappa_", colnames(Xg))])), collapse = "  "))
say("tau: median %.3f", median(post2[, "tau"]))
print(round(coef(glm(yg ~ Xg, family = binomial)), 3))

save_fig("ex02_horseshoe", width = 7.5, height = 3, function() {
  par(mfrow = c(1, 5), mar = c(4, 2, 2, 0.5))
  for (v in colnames(Xg)) {
    d <- density(post2[, v]); plot(d, main = v, xlab = expression(beta), ylab = "", yaxt = "n", lwd = 2,
                                    col = PALETTE["navy"], zero.line = FALSE)
    abline(v = 0, lty = 2, col = PALETTE["grey"])
  }
})


# ---- Exercise 3: spline regression of ozone on solar radiation ------------------------

section("Ex 3  ozone vs solar radiation")

aq <- na.omit(airquality[, c("Ozone", "Solar.R", "Wind", "Temp")])
x3 <- aq$Solar.R; y3 <- aq$Ozone; n3 <- length(y3)
basis3 <- function(J, x = x3) { B <- splines::bs(x, df = J, intercept = TRUE, Boundary.knots = range(x3)); colnames(B) <- paste0("B", 1:J); B }
fit_spline <- function(J, iters = 6000, burn = 1000) gibbs_lm(y3, basis3(J), prior_sd = 100, a = 0.1, b = 0.1, iters = iters, burn = burn)
ic3 <- t(sapply(4:12, function(J) {
  dr <- fit_spline(J); B <- basis3(J); mu <- dr[, 1:J] %*% t(B)
  ll <- dnorm(matrix(y3, nrow(mu), n3, byrow = TRUE), mu, dr[, "sigma"], log = TRUE)
  dev_mean <- -2 * sum(dnorm(y3, drop(B %*% colMeans(dr[, 1:J])), mean(dr[, "sigma"]), log = TRUE))
  c(J = J, dic(-2 * rowSums(ll), dev_mean), waic(ll))
}))
print(round(ic3[, c("J", "DIC", "pD", "WAIC", "pW")], 2))
J3 <- ic3[which.min(ic3[, "WAIC"]), "J"]
say("selected J = %d", J3)
d3 <- fit_spline(J3, iters = 21000)
grid3 <- seq(min(x3), max(x3), length.out = 200)
f_grid <- d3[, 1:J3] %*% t(basis3(J3, grid3))
f_obs  <- d3[, 1:J3] %*% t(basis3(J3))
band <- apply(f_grid, 2, quantile, c(0.025, 0.5, 0.975))
cover_f   <- mean(y3 >= apply(f_obs, 2, quantile, 0.025) & y3 <= apply(f_obs, 2, quantile, 0.975))
yrep3 <- f_obs + d3[, "sigma"] * matrix(rnorm(length(f_obs)), nrow(f_obs))
cover_ppd <- mean(y3 >= apply(yrep3, 2, quantile, 0.025) & y3 <= apply(yrep3, 2, quantile, 0.975))
slope_lin <- gibbs_lm(y3, cbind(1, (x3 - mean(x3)) / sd(x3)), iters = 6000, burn = 1000)[, 2]
say("coverage of the 95%% band for f: %.3f   coverage of 95%% prediction intervals: %.3f", cover_f, cover_ppd)
say("linear slope per SD of Solar.R: mean %.2f  95%% CI %s", mean(slope_lin), fmt_ci(quantile(slope_lin, c(0.025, 0.975)), 2))
say("range of f over the grid: posterior median of max f - min f = %.1f ppb  (95%% CI %s)",
    median(apply(f_grid, 1, function(f) max(f) - min(f))),
    fmt_ci(quantile(apply(f_grid, 1, function(f) max(f) - min(f)), c(0.025, 0.975)), 1))
say("P(f increasing on the whole grid) = %.3f", mean(apply(f_grid, 1, function(f) all(diff(f) > 0))))

save_fig("ex03_spline", width = 6.5, height = 4, function() {
  plot(x3, y3, pch = 19, cex = 0.7, col = PALETTE["grey"], xlab = "solar radiation (Langleys)", ylab = "ozone (ppb)")
  polygon(c(grid3, rev(grid3)), c(band[1, ], rev(band[3, ])), col = adjustcolor(PALETTE["red"], 0.25), border = NA)
  lines(grid3, band[2, ], lwd = 2.5, col = PALETTE["red"])
})


# ---- Exercise 4: heteroskedastic linear model ----------------------------------------

section("Ex 4  heteroskedastic regression")

X4 <- cbind(intercept = 1, scale(aq[, c("Solar.R", "Wind", "Temp")])); y4 <- aq$Ozone; p4 <- ncol(X4)
hetero <- function(iters = 25000, burn = 5000) {
  beta <- qr.solve(X4, y4); alpha <- c(log(var(y4 - X4 %*% beta)), 0, 0, 0)
  L <- t(chol((2.4^2 / p4) * 2 * solve(crossprod(X4))))
  log_lik_alpha <- function(a, r) sum(dnorm(r, 0, exp(drop(X4 %*% a) / 2), log = TRUE)) + sum(dnorm(a, 0, 10, log = TRUE))
  out <- matrix(NA, iters, 2 * p4, dimnames = list(NULL, c(paste0("beta_", colnames(X4)), paste0("alpha_", colnames(X4)))))
  acc <- 0
  for (it in seq_len(iters)) {
    w <- exp(-drop(X4 %*% alpha))                                   # beta | alpha: weighted regression
    prec <- crossprod(X4 * w, X4) + diag(1 / 100^2, p4)
    R <- chol(prec)
    beta <- drop(backsolve(R, forwardsolve(t(R), crossprod(X4 * w, y4))) + backsolve(R, rnorm(p4)))
    r <- drop(y4 - X4 %*% beta)
    cand <- alpha + drop(L %*% rnorm(p4))                           # alpha | beta: block Metropolis
    if (log(runif(1)) < log_lik_alpha(cand, r) - log_lik_alpha(alpha, r)) { alpha <- cand; acc <- acc + 1 }
    out[it, ] <- c(beta, alpha)
  }
  structure(out[-seq_len(burn), ], acceptance = acc / iters)
}
chains4 <- lapply(1:3, function(k) hetero())
post4 <- do.call(rbind, chains4)
say("alpha acceptance %.2f", attr(chains4[[1]], "acceptance"))
print(round(t(sapply(colnames(post4), function(p) {
  ch <- sapply(chains4, function(c_) c_[, p])
  c(post_summary(as.vector(ch), 3), Rhat = rhat(ch), ESS = sum(apply(ch, 2, ess)))
})), 3))
ll_het <- t(apply(post4[sample(nrow(post4), 5000), ], 1, function(s)
  dnorm(y4, X4 %*% s[1:p4], exp(drop(X4 %*% s[(p4 + 1):(2 * p4)]) / 2), log = TRUE)))
hom <- gibbs_lm(y4, X4, prior_sd = 100, a = 0.1, b = 0.1, iters = 11000, burn = 1000)
ll_hom <- t(apply(hom[sample(nrow(hom), 5000), ], 1, function(s) dnorm(y4, X4 %*% s[1:p4], s["sigma"], log = TRUE)))
print(round(rbind(heteroskedastic = waic(ll_het), constant_variance = waic(ll_hom)), 2))


# ---- Exercise 5: additive model with B-splines ---------------------------------------

section("Ex 5  additive model")

vars5 <- c("Solar.R", "Wind", "Temp")
bs5 <- function(v, x = aq[[v]]) {
  kn <- quantile(aq[[v]], seq(1, 5) / 6)                          # five interior knots
  B <- splines::bs(x, knots = kn, Boundary.knots = range(aq[[v]]))
  sweep(B, 2, colMeans(splines::bs(aq[[v]], knots = kn)))         # centre each basis on the data
}
B5 <- do.call(cbind, lapply(vars5, bs5)); grp5 <- rep(vars5, each = ncol(B5) / 3)
X5 <- cbind(intercept = 1, B5); colnames(X5) <- c("intercept", paste0(grp5, "_", sequence(rep(ncol(B5) / 3, 3))))
d5 <- gibbs_lm(aq$Ozone, X5, prior_sd = 100, a = 0.1, b = 0.1, iters = 21000, burn = 1000)
lin5 <- gibbs_lm(aq$Ozone, cbind(1, scale(aq[, vars5])), prior_sd = 100, a = 0.1, b = 0.1, iters = 11000, burn = 1000)
ll_of <- function(dr, X) t(apply(dr[sample(nrow(dr), 5000), ], 1, function(s) dnorm(aq$Ozone, X %*% s[1:ncol(X)], s["sigma"], log = TRUE)))
print(round(rbind(additive = waic(ll_of(d5, X5)), linear = waic(ll_of(lin5, cbind(1, scale(aq[, vars5]))))), 2))
curves5 <- lapply(vars5, function(v) {
  g <- seq(min(aq[[v]]), max(aq[[v]]), length.out = 150)
  F <- d5[, which(grp5 == v) + 1] %*% t(bs5(v, g))
  list(grid = g, band = apply(F, 2, quantile, c(0.025, 0.5, 0.975)), range = apply(F, 1, function(f) diff(range(f))))
})
names(curves5) <- vars5
for (v in vars5) say("%-8s effect range: median %.1f ppb (95%% CI %s)",
                     v, median(curves5[[v]]$range), fmt_ci(quantile(curves5[[v]]$range, c(0.025, 0.975)), 1))
# Is each effect non-linear?  Refit with that variable entering linearly (others smooth) and compare WAIC.
waic_lin <- t(sapply(vars5, function(v) {
  Xv <- cbind(intercept = 1, lin = as.numeric(scale(aq[[v]])), B5[, grp5 != v])
  dv <- gibbs_lm(aq$Ozone, Xv, prior_sd = 100, a = 0.1, b = 0.1, iters = 11000, burn = 1000)
  waic(ll_of(dv, Xv))
}))
rownames(waic_lin) <- paste(vars5, "linear, others smooth")
print(round(rbind(all_smooth = waic(ll_of(d5, X5)), waic_lin), 2))
# Is each variable important?  Drop it entirely and compare WAIC.
waic_drop <- t(sapply(vars5, function(v) {
  Xv <- cbind(intercept = 1, B5[, grp5 != v])
  dv <- gibbs_lm(aq$Ozone, Xv, prior_sd = 100, a = 0.1, b = 0.1, iters = 11000, burn = 1000)
  waic(ll_of(dv, Xv))
}))
rownames(waic_drop) <- paste("without", vars5)
print(round(waic_drop, 2))
save_fig("ex05_gam", width = 7.5, height = 2.9, function() {
  par(mfrow = c(1, 3), mar = c(4, 4.2, 1.2, 0.5))
  for (v in vars5) {
    cv <- curves5[[v]]
    plot(cv$grid, cv$band[2, ], type = "n", ylim = c(-40, 60), xlab = v, ylab = paste0("f(", v, ")"))
    polygon(c(cv$grid, rev(cv$grid)), c(cv$band[1, ], rev(cv$band[3, ])), col = adjustcolor(PALETTE["navy"], 0.25), border = NA)
    lines(cv$grid, cv$band[2, ], lwd = 2.5, col = PALETTE["navy"]); abline(h = 0, lty = 3, col = PALETTE["grey"])
    rug(aq[[v]], col = PALETTE["grey"])
  }
})


# ---- Exercise 7: BART for the motorcycle data -----------------------------------------

section("Ex 7  BART tuning for mcycle")

mc <- MASS::mcycle
xg7 <- matrix(seq(min(mc$times), max(mc$times), length.out = 200))
settings7 <- list("ntree=200, power=2, k=2 (default)" = c(200, 2, 2), "ntree=20" = c(20, 2, 2),
                  "ntree=50" = c(50, 2, 2), "power=0.5" = c(200, 0.5, 2), "power=5" = c(200, 5, 2),
                  "k=1" = c(200, 2, 1), "k=5" = c(200, 2, 5))
fits7 <- lapply(settings7, function(s) quiet(BART::wbart(matrix(mc$times), mc$accel, x.test = xg7,
                                                         ntree = s[1], power = s[2], k = s[3], ndpost = 800, nskip = 300)))
# 5-fold CV RMSE for each setting
folds7 <- sample(rep(1:5, length.out = nrow(mc)))
cv7 <- sapply(settings7, function(s) {
  sq <- numeric(nrow(mc))
  for (f in 1:5) {
    tr <- folds7 != f
    fb <- quiet(BART::wbart(matrix(mc$times[tr]), mc$accel[tr], x.test = matrix(mc$times[!tr]),
                            ntree = s[1], power = s[2], k = s[3], ndpost = 500, nskip = 200))
    sq[!tr] <- (mc$accel[!tr] - fb$yhat.test.mean)^2
  }
  sqrt(mean(sq))
})
print(data.frame(setting = names(settings7), CV_RMSE = round(cv7, 2),
                 sigma = round(sapply(fits7, function(f) mean(tail(f$sigma, 800))), 2)), row.names = FALSE)

save_fig("ex07_bart", width = 7.5, height = 3.4, function() {
  par(mfrow = c(1, 3), mar = c(4, 4.2, 2, 0.5))
  groups <- list(ntree = c(1, 2, 3), power = c(1, 4, 5), k = c(1, 6, 7))
  for (g in names(groups)) {
    plot(mc$times, mc$accel, pch = 19, cex = 0.5, col = PALETTE["grey"], xlab = "time (ms)", ylab = "acceleration (g)",
         main = paste("varying", g))
    for (m in seq_along(groups[[g]])) lines(xg7, fits7[[groups[[g]][m]]]$yhat.test.mean, lwd = 2, col = PALETTE[m], lty = m)
    legend("bottomright", names(settings7)[groups[[g]]], col = PALETTE[1:3], lty = 1:3, lwd = 2, bty = "n", cex = 0.7)
  }
})


# ---- Exercise 8: bike rentals, BART vs linear regression ------------------------------

section("Ex 8  bike rentals")

bike <- read.csv(file.path("data", "day.csv"))
make_X <- function(d) {
  f <- function(v, lv) sapply(lv[-1], function(l) as.numeric(d[[v]] == l))
  X <- cbind(yr = d$yr, holiday = d$holiday,                    # workingday = weekday & !holiday: redundant
             temp = d$temp, hum = d$hum, windspeed = d$windspeed,
             f("season", 1:4), f("mnth", 1:12), f("weekday", 0:6), f("weathersit", 1:3))
  colnames(X) <- c("yr", "holiday", "temp", "hum", "windspeed",
                   paste0("season", 2:4), paste0("mnth", 2:12), paste0("weekday", 1:6), paste0("weather", 2:3))
  X
}
X8 <- make_X(bike); y8 <- bike$cnt
say("predictors after processing: %d   cor(temp, atemp) = %.3f", ncol(X8), cor(bike$temp, bike$atemp))
Xs8 <- cbind(intercept = 1, scale(X8))
lm8 <- gibbs_lm(y8, Xs8, prior_sd = 1e4, a = 0.1, b = 0.1, iters = 6000, burn = 1000)
coef8 <- summary_table(lm8[, -ncol(lm8)], 1)
print(coef8[order(-abs(coef8[, "mean"])), ][1:10, ])
say("MLR residual sd: %.0f", mean(lm8[, "sigma"]))
bart8 <- quiet(BART::wbart(X8, y8, ndpost = 600, nskip = 300))
vc <- colMeans(bart8$varcount); vc <- sort(vc / sum(vc), decreasing = TRUE)
say("BART variable inclusion proportions (top 8): %s", paste(sprintf("%s %.3f", names(vc)[1:8], vc[1:8]), collapse = "  "))
say("BART residual sd: %.0f", mean(tail(bart8$sigma, 600)))

folds8 <- sample(rep(1:5, length.out = nrow(bike)))
cv8 <- sapply(c("MLR", "BART"), function(method) {
  sq <- cov <- numeric(nrow(bike))
  for (f in 1:5) {
    tr <- folds8 != f; te <- !tr
    if (method == "MLR") {
      ctr <- colMeans(X8[tr, ]); sdv <- apply(X8[tr, ], 2, sd); sdv[sdv == 0] <- 1
      Z <- function(X) cbind(1, sweep(sweep(X, 2, ctr), 2, sdv, "/"))
      dr <- gibbs_lm(y8[tr], Z(X8[tr, ]), prior_sd = 1e4, a = 0.1, b = 0.1, iters = 3000, burn = 500)
      mu <- dr[, 1:(ncol(X8) + 1)] %*% t(Z(X8[te, ])); sig <- dr[, "sigma"]
    } else {
      fb <- quiet(BART::wbart(X8[tr, ], y8[tr], x.test = X8[te, ], ndpost = 400, nskip = 200))
      mu <- fb$yhat.test; sig <- tail(fb$sigma, 400)
    }
    yrep <- mu + sig * matrix(rnorm(length(mu)), nrow(mu))
    lo <- apply(yrep, 2, quantile, 0.025); hi <- apply(yrep, 2, quantile, 0.975)
    sq[te] <- (y8[te] - colMeans(mu))^2; cov[te] <- y8[te] >= lo & y8[te] <= hi
  }
  c(RMSE = sqrt(mean(sq)), MSE = mean(sq), coverage = mean(cov))
})
print(round(cv8, 3))


# ---- Exercise 9: individualised treatment rules with BART -----------------------------

section("Ex 9  treatment rules")

set.seed(919)
X9 <- matrix(runif(500), 100, 5)
A9 <- rbinom(100, 1, X9[, 1])
g0 <- 2 + cos(10 * X9[, 3])
g1 <- 1 + 2 * X9[, 2] + cos(10 * X9[, 3])
Y9 <- rnorm(100, (1 - A9) * g0 + A9 * g1, 1)
set.seed(2026)
xtr9 <- cbind(A = A9, X9); colnames(xtr9) <- c("A", paste0("X", 1:5))
x1 <- xtr9; x1[, "A"] <- 1; x0 <- xtr9; x0[, "A"] <- 0
# S-learner: one BART with A as an extra input.
fit9 <- quiet(BART::wbart(xtr9, Y9, x.test = rbind(x1, x0), ndpost = 1500, nskip = 500))
cate_s <- fit9$yhat.test[, 1:100] - fit9$yhat.test[, 101:200]
# T-learner: separate BART models for the treated and the untreated.
colnames(X9) <- paste0("X", 1:5)
f1 <- quiet(BART::wbart(X9[A9 == 1, ], Y9[A9 == 1], x.test = X9, ndpost = 1500, nskip = 500))
f0 <- quiet(BART::wbart(X9[A9 == 0, ], Y9[A9 == 0], x.test = X9, ndpost = 1500, nskip = 500))
cate <- f1$yhat.test - f0$yhat.test
rule_true <- X9[, 2] > 0.5
for (nm in c("S-learner", "T-learner")) {
  cc <- if (nm == "S-learner") cate_s else cate
  say("%s: agreement with optimal rule %.2f   cor(CATE, 2*X2 - 1) = %.3f   range of posterior mean CATE %.2f to %.2f",
      nm, mean((colMeans(cc) > 0) == rule_true), cor(colMeans(cc), 2 * X9[, 2] - 1), min(colMeans(cc)), max(colMeans(cc)))
}
rule_fit <- colMeans(cate) > 0
say("treated: %d of 100   sigma (T-learner arms) %.3f %.3f", sum(A9), mean(tail(f1$sigma, 1500)), mean(tail(f0$sigma, 1500)))
say("in-sample fit of g: cor(fitted g1, true g1) %.3f   cor(fitted g0, true g0) %.3f", cor(colMeans(f1$yhat.test), g1), cor(colMeans(f0$yhat.test), g0))
save_fig("ex09_cate", width = 6.5, height = 3.8, function() {
  q <- apply(cate, 2, quantile, c(0.025, 0.975))
  plot(X9[, 2], colMeans(cate), pch = 19, col = ifelse(rule_fit == rule_true, PALETTE["navy"], PALETTE["red"]),
       ylim = range(q), xlab = expression(X[2]), ylab = expression(hat(g)[1] - hat(g)[0]))
  segments(X9[, 2], q[1, ], X9[, 2], q[2, ], col = adjustcolor(PALETTE["navy"], 0.3))
  abline(-1, 2, lty = 2, lwd = 2, col = PALETTE["green"]); abline(h = 0, v = 0.5, lty = 3, col = PALETTE["grey"])
  legend("topleft", c("rule agrees with optimal", "rule disagrees", "true effect 2X2 - 1"), bty = "n", cex = 0.8,
         pch = c(19, 19, NA), lty = c(NA, NA, 2), col = c(PALETTE["navy"], PALETTE["red"], PALETTE["green"]))
})


# ---- Exercise 10: black legislators and email responses --------------------------------

section("Ex 10  legislator responses")

bp <- as.data.frame(causaldata::black_politicians)
bp <- bp[bp$treat_out == 1, ]
conf10 <- c("totalpop", "medianhhincom", "black_medianhh", "white_medianhh", "blackpercent", "statessquireindex",
            "nonblacknonwhite", "urbanpercent", "leg_senator", "leg_democrat", "south")
y10 <- bp$responded; a10 <- bp$leg_black
say("emails %d   response rate: black legislators %.3f, others %.3f", nrow(bp), mean(y10[a10 == 1]), mean(y10[a10 == 0]))
logit_fit <- function(X) {
  g <- glm(y10 ~ X - 1, family = binomial)
  lp <- function(b) { eta <- drop(X %*% b); sum(y10 * eta - log1p(exp(eta))) + sum(dnorm(b, 0, 10, log = TRUE)) }
  metropolis_block(lp, setNames(coef(g), colnames(X)), (2.4^2 / ncol(X)) * vcov(g), iters = 22000, burn = 2000)
}
Xa <- cbind(intercept = 1, leg_black = a10)
Xb <- cbind(intercept = 1, leg_black = a10, scale(bp[, conf10]))
fa <- logit_fit(Xa); fb10 <- logit_fit(Xb)
for (nm in c("treatment only", "with confounders")) {
  f <- if (nm == "treatment only") fa else fb10
  say("%-17s log OR %.3f  95%% CI %s   OR %.3f   P(OR > 1) = %.3f   acceptance %.2f", nm, mean(f[, "leg_black"]),
      fmt_ci(quantile(f[, "leg_black"], c(0.025, 0.975))), exp(median(f[, "leg_black"])), mean(f[, "leg_black"] > 0),
      attr(f, "acceptance"))
}
print(round(summary_table(fb10, 3), 3))
xb <- as.matrix(cbind(leg_black = a10, bp[, conf10]))
x1b <- xb; x1b[, "leg_black"] <- 1; x0b <- xb; x0b[, "leg_black"] <- 0
lb <- quiet(BART::lbart(xb, y10, x.test = rbind(x1b, x0b), ntree = 50, ndpost = 400, nskip = 200))
nb <- nrow(xb)
mate <- rowMeans(lb$prob.test[, 1:nb] - lb$prob.test[, (nb + 1):(2 * nb)])
say("BART MATE (difference in response probability): mean %.4f  95%% CI %s   P(MATE > 0) = %.3f",
    mean(mate), fmt_ci(quantile(mate, c(0.025, 0.975)), 4), mean(mate > 0))
mate_lin <- apply(fb10, 1, function(b) mean(plogis(cbind(1, 1, scale(bp[, conf10])) %*% b) - plogis(cbind(1, 0, scale(bp[, conf10])) %*% b)))
say("linear-logistic MATE: mean %.4f  95%% CI %s", mean(mate_lin), fmt_ci(quantile(mate_lin, c(0.025, 0.975)), 4))
say("naive difference in response rates: %.4f", mean(y10[a10 == 1]) - mean(y10[a10 == 0]))
say("confounder imbalance (standardised mean difference, black - other): %s",
    paste(sprintf("%s %.2f", conf10, sapply(conf10, function(v) (mean(bp[a10 == 1, v]) - mean(bp[a10 == 0, v])) / sd(bp[[v]]))), collapse = "  "))
save_fig("ex10_mate", width = 6.5, height = 3.6, function() {
  plot(density(mate), lwd = 2.5, col = PALETTE["navy"], main = "", xlab = "MATE (difference in response probability)",
       zero.line = FALSE, xlim = range(mate, mate_lin))
  lines(density(mate_lin), lwd = 2, lty = 2, col = PALETTE["red"])
  abline(v = 0, lty = 3, col = PALETTE["grey"])
  legend("topright", c("BART (lbart)", "linear logistic"), lty = 1:2, lwd = 2, col = PALETTE[c("navy", "red")], bty = "n")
})


# ---- Exercise 11: equation learning for an SIR epidemic --------------------------------

section("Ex 11  SIR equation learning")

# covid19.analytics::generate.SIR.model() needs live data feeds; we simulate an SIR epidemic of the
# same size directly (RK4 integration), with N = 3.5 million and a 90-day window.
Npop <- 3.5e6; beta_true <- 0.40; gamma_true <- 0.25; days <- 90
sir_rhs <- function(s, i) c(-beta_true * s * i / Npop, beta_true * s * i / Npop - gamma_true * i)
state <- c(Npop - 50, 50); traj <- matrix(NA, days, 2)
for (t in 1:days) {
  traj[t, ] <- state
  for (k in 1:10) {                               # 10 RK4 steps per day
    h <- 0.1
    k1 <- sir_rhs(state[1], state[2]); k2 <- sir_rhs(state[1] + h / 2 * k1[1], state[2] + h / 2 * k1[2])
    k3 <- sir_rhs(state[1] + h / 2 * k2[1], state[2] + h / 2 * k2[2]); k4 <- sir_rhs(state[1] + h * k3[1], state[2] + h * k3[2])
    state <- state + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
}
set.seed(1)
S11 <- traj[, 1] * rgamma(90, 100, 100)
I11 <- traj[, 2] * rgamma(90, 100, 100)
set.seed(2026)
# The SIR equation dI/dt = beta S I / N - gamma I is linear in (beta, gamma) after dividing by I:
#   d log I / dt = beta (S / N) - gamma.
# Multiplicative noise is additive on the log scale. Smooth log I and log S (smoothness chosen by
# GCV), differentiate the smooth, then do Bayesian regression on a library of candidate terms.
tt <- 1:days
smI <- smooth.spline(tt, log(I11)); smS <- smooth.spline(tt, log(S11))
say("GCV-chosen smoothness: df %.1f for log I, %.1f for log S", smI$df, smS$df)
dlogI <- predict(smI, tt, deriv = 1)$y
I_s <- exp(predict(smI, tt)$y); S_s <- exp(predict(smS, tt)$y)
lib <- cbind(constant = 1, S_over_N = S_s / Npop, I_over_N = I_s / Npop)
y11 <- 100 * dlogI                      # percent per day, so the InvGamma(0.1, 0.1) prior on sigma^2 is vague
d11_full <- gibbs_lm(y11, lib, prior_sd = 1000, a = 0.1, b = 0.1, iters = 11000, burn = 1000)
print(round(summary_table(d11_full / 100, 5), 5))
d11 <- gibbs_lm(y11, lib[, c("constant", "S_over_N")], prior_sd = 1000, a = 0.1, b = 0.1, iters = 11000, burn = 1000) / 100
beta_hat <- d11[, "S_over_N"]; gamma_hat <- -d11[, "constant"]; R0 <- beta_hat / gamma_hat
say("cor(beta, gamma) in the posterior: %.3f   S/N falls to %.2f by day 90", cor(beta_hat, gamma_hat), min(traj[, 1]) / Npop)
say("beta: mean %.4f  95%% CI %s  (true %.2f)", mean(beta_hat), fmt_ci(quantile(beta_hat, c(0.025, 0.975)), 4), beta_true)
say("gamma: mean %.4f  95%% CI %s  (true %.2f)", mean(gamma_hat), fmt_ci(quantile(gamma_hat, c(0.025, 0.975)), 4), gamma_true)
say("R0 = beta/gamma: mean %.3f  95%% CI %s  (true %.2f)", mean(R0), fmt_ci(quantile(R0, c(0.025, 0.975)), 3), beta_true / gamma_true)
say("peak infections: day %d (%.0f)", which.max(traj[, 2]), max(traj[, 2]))
save_fig("ex11_sir", width = 7.5, height = 3.2, function() {
  par(mfrow = c(1, 2), mar = c(4, 4.6, 1.2, 0.5))
  plot(tt, I11 / 1000, pch = 19, cex = 0.5, col = PALETTE["grey"], xlab = "day", ylab = "infected (thousands)")
  lines(tt, traj[, 2] / 1000, lwd = 2, col = PALETTE["navy"]); lines(tt, I_s / 1000, lwd = 2, lty = 2, col = PALETTE["red"])
  legend("topleft", c("noisy I(t)", "true I(t)", "smoothed"), pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 2),
         col = c(PALETTE["grey"], PALETTE["navy"], PALETTE["red"]), bty = "n", cex = 0.8)
  plot(density(R0), lwd = 2.5, col = PALETTE["navy"], main = "", xlab = expression(R[0] == beta / gamma), zero.line = FALSE)
  abline(v = beta_true / gamma_true, lty = 2, col = PALETTE["red"])
})
