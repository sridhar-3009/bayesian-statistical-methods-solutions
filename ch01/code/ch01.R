# =============================================================================
#  Chapter 1 -- Basics of Bayesian inference
#  Reproduces every number and figure in ch01/ch01.tex.
#
#  Run from the repository root:   Rscript ch01/code/ch01.R
#  Output:  figures -> ch01/figures/*.pdf     numbers -> console
#  Needs base R only (no extra packages).
# =============================================================================

set.seed(2026)
FIG_DIR <- file.path("ch01", "figures")


# ---- helpers ----------------------------------------------------------------

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


# ---- Exercise 8: conditionals of a bivariate t(1) density -------------------

section("Ex 8  conditional of X1 | X2")

cond_density <- function(x1, x2) {
  a2 <- 1 + x2^2
  0.5 * a2 * (a2 + x1^2)^(-1.5)
}

x1_grid <- seq(-10, 10, length.out = 801)
x2_vals <- -3:3

save_fig("ex08_conditionals", function() {
  plot(NA, xlim = range(x1_grid), ylim = c(0, 0.5),
       xlab = expression(x[1]), ylab = expression(f(x[1] ~ "|" ~ x[2])))
  for (k in seq_along(x2_vals))
    lines(x1_grid, cond_density(x1_grid, x2_vals[k]), col = PALETTE[k], lwd = 2)
  legend("topright", legend = paste("x2 =", x2_vals), col = PALETTE, lwd = 2,
         bty = "n", cex = 0.85)
})

for (v in x2_vals)
  say("x2 = %2d   integral = %.4f   peak height = %.4f",
      v, integrate(cond_density, -Inf, Inf, x2 = v)$value, cond_density(0, v))


# ---- Exercise 13: Bayesian spell checker -------------------------------------

section("Ex 13  spell checker")

DICTIONARY <- c("fun", "sun", "sit", "sat", "fan", "for")

n_mismatches <- function(a, b) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])

# Posterior over intended words.  `strict = TRUE` uses theta/25 for a specific
# wrong letter (rule iii); FALSE follows the likelihood written in the book.
spell_check <- function(typed, alpha, theta, strict = FALSE) {
  prior <- ifelse(DICTIONARY == "for", alpha, 1)
  k     <- sapply(DICTIONARY, n_mismatches, b = typed)
  error <- if (strict) theta / 25 else theta
  post  <- prior * error^k * (1 - theta)^(3 - k)
  post / sum(post)
}

cases <- list(a = c(alpha = 2,  theta = 0.10),
              b = c(alpha = 50, theta = 0.10),
              c = c(alpha = 2,  theta = 0.95))

for (strict in c(FALSE, TRUE)) {
  say("%s", if (strict) "\nstrict rule (theta/25 per specific wrong letter):"
            else "likelihood exactly as written in the exercise:")
  for (id in names(cases)) for (typed in c("sun", "the", "foo")) {
    post <- spell_check(typed, cases[[id]]["alpha"], cases[[id]]["theta"], strict)
    say("(%s) %s : %s", id, typed, paste(sprintf("%s %.4f", DICTIONARY, post), collapse = "  "))
  }
}


# ---- Exercise 14: zero-inflated Poisson --------------------------------------

section("Ex 14  zero-inflated Poisson")

n <- 10000; p <- 0.2; lambda <- 4
U <- rbinom(n, 1, p)
V <- rpois(n, lambda)
Y <- U * V

zip_pmf <- function(y) (1 - p) * (y == 0) + p * dpois(y, lambda)

save_fig("ex14_zip", function() {
  support <- 0:max(Y)
  freq    <- table(factor(Y, levels = support)) / n
  centers <- barplot(freq, col = PALETTE["navy"], border = NA, xlab = "y",
                     ylab = "Proportion", names.arg = support)
  points(centers, zip_pmf(support), pch = 19, col = PALETTE["red"])
  legend("topright", c("simulated", "exact PMF"), bty = "n",
         fill = c(PALETTE["navy"], NA), border = NA,
         pch = c(NA, 19), col = c(NA, PALETTE["red"]))
})

say("simulated   P(Y=0) = %.4f   P(Y=4) = %.4f", mean(Y == 0), mean(Y == 4))
say("exact       P(Y=0) = %.6f P(Y=4) = %.6f   P(U=0 | Y=0) = %.6f",
    zip_pmf(0), zip_pmf(4), (1 - p) / zip_pmf(0))


# ---- Exercise 15: occupancy with imperfect detection -------------------------

section("Ex 15  species occupancy")

n_trees     <- 0:60
p_absent    <- 1 / (1 + 0.9^n_trees)          # theta = 0.5, lambda = 0.1
n_needed    <- min(n_trees[p_absent >= 0.95])

say("smallest n with P(X1 = 0 | X2 = 0) >= 0.95 : %d", n_needed)
say("n = 27 -> %.4f    n = 28 -> %.4f", p_absent[n_trees == 27], p_absent[n_trees == 28])

save_fig("ex15_occupancy", width = 6, height = 3.8, function() {
  plot(n_trees, p_absent, type = "b", pch = 19, cex = 0.6, col = PALETTE["navy"],
       ylim = c(0.5, 1), xlab = "n (trees sampled)",
       ylab = expression(Prob(X[1] == 0 ~ "|" ~ X[2] == 0)))
  abline(h = 0.95, lty = 2, col = PALETTE["red"])
  abline(v = n_needed, lty = 3, col = PALETTE["grey"])
})


# ---- Exercise 16: Poisson rate of species discoveries ------------------------

section("Ex 16  bivalve discoveries")

discoveries <- c(64, 13, 33, 18, 30, 20)
shape <- sum(discoveries) + 1                 # posterior Gamma(shape, rate)
rate  <- length(discoveries)

lambda_grid <- seq(15, 45, length.out = 1001)
log_post    <- sum(discoveries) * log(lambda_grid) - rate * lambda_grid
grid_post   <- exp(log_post - max(log_post))
grid_post   <- grid_post / sum(grid_post * diff(lambda_grid)[1])

save_fig("ex16_lambda", width = 6, height = 3.8, function() {
  plot(lambda_grid, grid_post, type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(lambda), ylab = "Posterior density")
  lines(lambda_grid, dgamma(lambda_grid, shape, rate), lty = 2, lwd = 1.5, col = PALETTE["red"])
  legend("topright", c("grid posterior", "Gamma(179, 6)"), lty = 1:2, lwd = 2,
         col = PALETTE[c("navy", "red")], bty = "n")
})

say("sum y = %d   mean = %.3f   sd = %.3f   95%% CI = (%.2f, %.2f)",
    sum(discoveries), shape / rate, sqrt(shape) / rate,
    qgamma(0.025, shape, rate), qgamma(0.975, shape, rate))
say("data sd = %.2f   Gamma mass above 100 = %.1e",
    sd(discoveries), pgamma(100, shape, rate, lower.tail = FALSE))


# ---- Exercise 17: posterior of a correlation ---------------------------------

section("Ex 17  correlation")

pairs <- rbind(c(-3.3, -2.6), c(0.1, -0.2), c(-1.1, -1.5),
               c( 2.7,  1.5), c(2.0,  1.9), c(-0.4, -0.3))
Sxx <- sum(pairs[, 1]^2)
Syy <- sum(pairs[, 2]^2)
Sxy <- sum(pairs[, 1] * pairs[, 2])

rho_grid <- seq(-0.999, 0.999, length.out = 3999)
step     <- diff(rho_grid)[1]
log_lik  <- -3 * log(1 - rho_grid^2) - (Sxx - 2 * rho_grid * Sxy + Syy) / (2 * (1 - rho_grid^2))
rho_post <- exp(log_lik - max(log_lik))
rho_post <- rho_post / sum(rho_post * step)
rho_cdf  <- cumsum(rho_post) * step
quantile_at <- function(q) rho_grid[which.min(abs(rho_cdf - q))]

say("Sxx = %.2f  Syy = %.2f  Sxy = %.2f   mean (x-y)^2 = %.3f",
    Sxx, Syy, Sxy, (Sxx + Syy - 2 * Sxy) / 6)
say("mode = %.3f  mean = %.3f  95%% CI = (%.3f, %.3f)  sample cor = %.3f",
    rho_grid[which.max(rho_post)], sum(rho_grid * rho_post) * step,
    quantile_at(0.025), quantile_at(0.975), cor(pairs[, 1], pairs[, 2]))

save_fig("ex17_scatter", width = 4.2, height = 4.2, function() {
  plot(pairs, pch = 19, cex = 1.3, col = PALETTE["navy"], xlab = "x", ylab = "y",
       xlim = c(-4, 4), ylim = c(-4, 4))
  abline(h = 0, v = 0, col = "grey85")
})

save_fig("ex17_rho", height = 4.2, function() {
  plot(rho_grid, rho_post, type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(rho), ylab = "Posterior density")
  abline(h = 0.5, lty = 2, col = PALETTE["grey"])
  legend("topleft", c("posterior", "Uniform(-1, 1) prior"), lty = 1:2,
         lwd = c(2.5, 1), col = PALETTE[c("navy", "grey")], bty = "n")
})


# ---- Exercise 18: classifying deforested pixels ------------------------------

section("Ex 18  NDVI classifier")

prior_deforested <- 0.1
f_forest     <- function(x) dbeta(x, 25, 10)
f_deforested <- function(x) dbeta(x, 10, 15)
p_deforested <- function(x) {
  num <- prior_deforested * f_deforested(x)
  num / (num + (1 - prior_deforested) * f_forest(x))
}

cutoff <- uniroot(function(x) p_deforested(x) - 0.9, c(0.2, 0.6))$root
ndvi   <- seq(0.001, 0.999, length.out = 999)

save_fig("ex18_densities", function() {
  marginal <- 0.9 * f_forest(ndvi) + 0.1 * f_deforested(ndvi)
  plot(ndvi, f_forest(ndvi), type = "l", lwd = 2.5, col = PALETTE["green"],
       xlab = "NDVI", ylab = "Density")
  lines(ndvi, f_deforested(ndvi), lwd = 2.5, col = PALETTE["gold"])
  lines(ndvi, marginal, lwd = 2.5, lty = 2, col = PALETTE["navy"])
  legend("topleft", c("forest: Beta(25, 10)", "deforested: Beta(10, 15)", "marginal (0.9 / 0.1)"),
         col = PALETTE[c("green", "gold", "navy")], lty = c(1, 1, 2), lwd = 2.5, bty = "n")
})

save_fig("ex18_prob", function() {
  plot(ndvi, p_deforested(ndvi), type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = "NDVI", ylab = "Prob(deforested | NDVI)")
  abline(h = 0.9, lty = 2, col = PALETTE["red"])
  abline(v = cutoff, lty = 3, col = PALETTE["grey"])
})

say("classify as deforested when NDVI < %.4f", cutoff)
say("P(NDVI < cutoff | forest) = %.4f   P(NDVI < cutoff | deforested) = %.4f",
    pbeta(cutoff, 25, 10), pbeta(cutoff, 10, 15))


# ---- Exercise 19: number of customers ----------------------------------------

section("Ex 19  number of customers")

n_max <- 30
customer_posterior <- function(y, theta) {
  w <- dpois(0:n_max, 5) * dbinom(y, 0:n_max, theta)
  w / sum(w)
}

save_fig("ex19_customers", width = 7.5, height = 5, function() {
  par(mfrow = c(2, 3), mar = c(4, 4, 2.2, 0.6))
  for (theta in c(0.2, 0.5)) for (y in c(0, 5, 10)) {
    post <- customer_posterior(y, theta)
    barplot(post, names.arg = 0:n_max, col = PALETTE["navy"], border = NA,
            ylim = c(0, 0.6), xlab = "n", ylab = "Posterior",
            main = bquote(Y == .(y) ~ "," ~ theta == .(theta)))
    mu <- sum(0:n_max * post)
    say("Y = %2d  theta = %.1f   E(n|Y) = %5.2f   sd = %.2f   (theory: %5.2f)",
        y, theta, mu, sqrt(sum((0:n_max)^2 * post) - mu^2), y + 5 * (1 - theta))
  }
})


# ---- Exercise 20: seedling survival ------------------------------------------

section("Ex 20  seedling survival")

a <- 3; b <- 9                                # Beta(1 + 2, 1 + 8)

equal_tailed <- qbeta(c(0.025, 0.975), a, b)
hpd_width    <- function(lo) qbeta(pbeta(lo, a, b) + 0.95, a, b) - lo
hpd_lower    <- optimize(hpd_width, c(0, qbeta(0.05, a, b)))$minimum
hpd          <- c(hpd_lower, qbeta(pbeta(hpd_lower, a, b) + 0.95, a, b))

say("mean = %.4f   sd = %.4f", a / (a + b), sqrt(a * b / ((a + b)^2 * (a + b + 1))))
say("equal-tailed (%.4f, %.4f)  width %.4f", equal_tailed[1], equal_tailed[2], diff(equal_tailed))
say("HPD          (%.4f, %.4f)  width %.4f  (density at ends %.3f, %.3f)",
    hpd[1], hpd[2], diff(hpd), dbeta(hpd[1], a, b), dbeta(hpd[2], a, b))
say("P(at least one of 10 survives) = %.5f", 1 - beta(a, b + 10) / beta(a, b))

save_fig("ex20_intervals", function() {
  theta <- seq(0, 0.8, length.out = 801)
  plot(theta, dbeta(theta, a, b), type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(theta), ylab = "Posterior density")
  segments(equal_tailed[1], 0.25, equal_tailed[2], 0.25, lwd = 4, col = PALETTE["red"])
  segments(hpd[1], 0.55, hpd[2], 0.55, lwd = 4, col = PALETTE["green"])
  legend("topright", c("Beta(3, 9) posterior", "equal-tailed 95%", "HPD 95%"),
         col = PALETTE[c("navy", "red", "green")], lwd = c(2.5, 4, 4), bty = "n")
})


# ---- Exercise 22: NBA clutch free throws -------------------------------------

section("Ex 22  clutch free throws")

nba <- data.frame(
  player  = c("Russell Westbrook", "James Harden", "Kawhi Leonard", "LeBron James",
              "Isaiah Thomas", "Stephen Curry", "Giannis Antetokounmpo", "John Wall",
              "Anthony Davis", "Kevin Durant"),
  overall = c(0.845, 0.847, 0.880, 0.674, 0.909, 0.898, 0.770, 0.801, 0.802, 0.875),
  makes   = c(64, 72, 55, 27, 75, 24, 28, 66, 40, 13),
  tries   = c(75, 95, 63, 39, 83, 26, 41, 82, 54, 16)
)

clutch_summary <- function(a0 = 1, b0 = 1) {
  A <- a0 + nba$makes
  B <- b0 + nba$tries - nba$makes
  data.frame(player  = nba$player,
             overall = nba$overall,
             mean    = round(A / (A + B), 3),
             lower   = round(qbeta(0.025, A, B), 3),
             upper   = round(qbeta(0.975, A, B), 3),
             p_below = round(pbeta(nba$overall, A, B), 3))
}

print(clutch_summary(), row.names = FALSE)
cat("\nprior sensitivity of P(theta < overall):\n")
print(data.frame(player   = nba$player,
                 uniform  = clutch_summary(1, 1)$p_below,
                 jeffreys = clutch_summary(0.5, 0.5)$p_below,
                 beta22   = clutch_summary(2, 2)$p_below), row.names = FALSE)

save_fig("ex22_clutch", width = 7.5, height = 6, function() {
  par(mfrow = c(5, 2), mar = c(2.2, 2.5, 1.6, 0.5), mgp = c(1.3, 0.4, 0), cex.main = 0.9)
  theta <- seq(0.4, 1, length.out = 601)
  for (i in seq_len(nrow(nba))) {
    plot(theta, dbeta(theta, 1 + nba$makes[i], 1 + nba$tries[i] - nba$makes[i]),
         type = "l", lwd = 2, col = PALETTE["navy"], main = nba$player[i],
         xlab = "", ylab = "", yaxt = "n")
    abline(v = nba$overall[i], lty = 2, lwd = 2, col = PALETTE["red"])
  }
})


# ---- Exercise 23: Federalist papers ------------------------------------------

section("Ex 23  Federalist papers")

# Posterior probability of Hamilton for `y` uses of "upon" in `kwords` x 1,000 words.
p_hamilton <- function(y, kwords = 1, prior = 0.5) {
  lik_h <- dpois(y, 3.24 * kwords)
  lik_m <- dpois(y, 0.23 * kwords)
  prior * lik_h / (prior * lik_h + (1 - prior) * lik_m)
}

say("Y =  3 in  1,000 words : %.5f", p_hamilton(3))
say("Y =  1 in  1,000 words : %.5f", p_hamilton(1))
say("Y = 30 in 10,000 words : 1 - %.1e", 1 - p_hamilton(30, kwords = 10))

upon <- 0:20
print(data.frame(Y = upon, p_hamilton = signif(p_hamilton(upon), 5)), row.names = FALSE)

save_fig("ex23_federalist", function() {
  plot(upon, p_hamilton(upon), type = "b", pch = 19, col = PALETTE["navy"], ylim = c(0, 1),
       xlab = "Y = occurrences of \"upon\" per 1,000 words", ylab = "Prob(Hamilton | Y)")
  abline(h = 0.5, lty = 2, col = PALETTE["red"])
})


# ---- Exercise 24: frequentist vs Bayesian test -------------------------------

section("Ex 24  tutoring program")

z <- (0.725 - 0.7) / sqrt(0.7 * 0.3 / 1000)
say("z = %.4f   normal p = %.4f   exact binomial p = %.4f",
    z, pnorm(z, lower.tail = FALSE), pbinom(724, 1000, 0.7, lower.tail = FALSE))
say("Beta(726, 276):  P(H0 | Y) = %.4f   mean = %.4f   95%% CI = (%.4f, %.4f)",
    pbeta(0.7, 726, 276), 726 / 1002, qbeta(0.025, 726, 276), qbeta(0.975, 726, 276))

save_fig("ex24_posterior", function() {
  theta <- seq(0.66, 0.79, length.out = 600)
  dens  <- dbeta(theta, 726, 276)
  plot(theta, dens, type = "l", lwd = 2.5, col = PALETTE["navy"],
       xlab = expression(theta), ylab = "Posterior density")
  null <- theta <= 0.7
  polygon(c(theta[null], 0.7), c(dens[null], 0), col = adjustcolor(PALETTE["red"], 0.35), border = NA)
  abline(v = 0.7, lty = 2, col = PALETTE["red"])
})


# ---- Exercise 25: Moderna vaccine efficacy -----------------------------------

section("Ex 25  vaccine efficacy (NEJM 2021 primary analysis)")

trial <- list(placebo = c(infected = 185, n = 14073),
              vaccine = c(infected = 11,  n = 14134))
draws <- 1e5

posterior_draws <- function(arm) rbeta(draws, 1 + arm["infected"], 1 + arm["n"] - arm["infected"])
theta_p  <- posterior_draws(trial$placebo)
theta_v  <- posterior_draws(trial$vaccine)
efficacy <- 1 - theta_v / theta_p

say("VE mean = %.4f   median = %.4f   95%% CI = (%.4f, %.4f)",
    mean(efficacy), median(efficacy), quantile(efficacy, 0.025), quantile(efficacy, 0.975))
say("P(VE > 0.3) = %.5f   P(VE > 0.9) = %.4f", mean(efficacy > 0.3), mean(efficacy > 0.9))

save_fig("ex25_efficacy", function() {
  hist(efficacy, breaks = 80, freq = FALSE, col = PALETTE["navy"], border = "white", main = "",
       xlab = expression("Vaccine efficacy" ~ 1 - theta[V] / theta[P]), ylab = "Posterior density")
  abline(v = quantile(efficacy, c(0.025, 0.975)), lty = 2, lwd = 2, col = PALETTE["red"])
})
