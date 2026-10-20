# =============================================================================
#  Chapter 3 -- JAGS versions of the models fitted by hand in ch03.R
#  Run by GitHub Actions (see .github/workflows/build.yml), not on the laptop.
#
#  Run from the repository root:   Rscript ch03/code/ch03_jags.R
#  Needs JAGS + the R packages rjags, coda and MASS.
# =============================================================================

suppressPackageStartupMessages({ library(rjags); library(coda) })
set.seed(2026)

# Compile, burn in and sample a JAGS model; returns a coda mcmc.list.
run_jags <- function(model, data, params, inits = NULL, chains = 3, burn = 5000, iters = 20000) {
  m <- jags.model(textConnection(model), data = data, inits = inits,
                  n.chains = chains, quiet = TRUE)
  update(m, burn, progress.bar = "none")
  coda.samples(m, params, n.iter = iters, progress.bar = "none")
}

report <- function(title, samples) {
  cat("\n==== ", title, " ====\n", sep = "")
  draws <- as.matrix(samples)
  tab <- t(apply(draws, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975)))))
  tab <- cbind(tab, ESS = effectiveSize(samples))
  if (nchain(samples) > 1) tab <- cbind(tab, Rhat = gelman.diag(samples, multivariate = FALSE)$psrf[, 1])
  print(signif(tab, 4))
  invisible(draws)
}


# ---- Exercise 4(e) -------------------------------------------------------------

for (a in c(10, 1)) report(sprintf("Ex 4(e)  a = %d", a), run_jags("
model{
  for(i in 1:n){
    Y[i] ~ dnorm(0, prec[i])         # JAGS uses precision = 1 / variance
    prec[i] ~ dgamma(a, b)           # sigma2[i] ~ InvGamma(a, b)
    sigma2[i] <- 1 / prec[i]
  }
  b ~ dgamma(1, 1)
}", list(Y = 1:10, n = 10, a = a), c("b", "sigma2[1]", "sigma2[10]")))


# ---- Exercise 5(d) -------------------------------------------------------------

ex5 <- read.csv(file.path("ch03", "code", "ex05_data.csv"))
report("Ex 5(d)  two-group normal", run_jags("
model{
  for(i in 1:N){
    Y[i] ~ dnorm(mu + delta * group[i], tau)
  }
  mu    ~ dnorm(0, 1 / 100^2)
  delta ~ dnorm(0, 1 / 100^2)
  tau   ~ dgamma(0.01, 0.01)
  sigma <- 1 / sqrt(tau)
}", list(Y = ex5$y, group = ex5$group, N = nrow(ex5)), c("mu", "delta", "sigma")))


# ---- Exercise 7(e) -------------------------------------------------------------

report("Ex 7(e)  NBA clutch shooting", run_jags("
model{
  for(i in 1:10){
    Y[i] ~ dbin(theta[i], n[i])
    theta[i] ~ dbeta(exp(m) * q[i], exp(m) * (1 - q[i]))
  }
  m ~ dnorm(0, 0.1)                  # variance 10
}", list(Y = c(64, 72, 55, 27, 75, 24, 28, 66, 40, 13), n = c(75, 95, 63, 39, 83, 26, 41, 82, 54, 16),
         q = c(0.845, 0.847, 0.880, 0.674, 0.909, 0.898, 0.770, 0.801, 0.802, 0.875)),
  c("theta", "m"), iters = 30000))


# ---- Exercise 9(b) -------------------------------------------------------------

galaxies <- MASS::galaxies
report("Ex 9(b)  t model for galaxies", run_jags("
model{
  for(i in 1:n){
    Y[i] ~ dt(mu, tau, k)
  }
  mu  ~ dnorm(0, 1 / 10000^2)
  tau ~ dgamma(0.01, 0.01)
  k   ~ dunif(1, 30)
  sigma <- 1 / sqrt(tau)
}", list(Y = galaxies, n = length(galaxies)), c("mu", "sigma", "k"),
  inits = list(list(mu = 15000, tau = 1 / 5000^2, k = 25), list(mu = 25000, tau = 1 / 1000^2, k = 3),
               list(mu = 20000, tau = 1 / 3000^2, k = 10))))


# ---- Exercise 11 ---------------------------------------------------------------

report("Ex 11  Reggie Jackson", run_jags("
model{
  Y_reg ~ dpois(N_reg * lambda_reg)
  Y_ws  ~ dpois(N_ws * lambda_ws)
  lambda_reg ~ dunif(0, 10)
  lambda_ws  ~ dunif(0, 10)
  ratio <- lambda_ws / lambda_reg
}", list(Y_reg = 563, N_reg = 2820, Y_ws = 10, N_ws = 27), c("lambda_reg", "lambda_ws", "ratio")))


# ---- Exercise 12 ---------------------------------------------------------------

draws12 <- report("Ex 12  Poisson trend", run_jags("
model{
  for(t in 1:6){
    Y[t] ~ dpois(lambda[t])
    log(lambda[t]) <- alpha + beta * t
  }
  alpha ~ dnorm(0, 0.01)             # variance 10^2
  beta  ~ dnorm(0, 0.01)
}", list(Y = c(64, 13, 33, 18, 30, 20)), c("alpha", "beta")))
cat("P(beta < 0) =", mean(draws12[, "beta"] < 0), "\n")


# ---- Exercise 14 ---------------------------------------------------------------

draws14 <- report("Ex 14  three-arm trial", run_jags("
model{
  for(j in 1:3){
    Y[j] ~ dbin(p[j], 100)
    p[j] ~ dunif(0, 1)
  }
  low_best <- step(p[2] - p[1]) * step(p[2] - p[3])
}", list(Y = c(52, 60, 54)), c("p", "low_best")))


# ---- Exercise 15(e) ------------------------------------------------------------

set.seed(27695)
B <- rbinom(30, 1, 0.5)
Y15 <- rnorm(30, B * 4, 1)
report("Ex 15(e)  normal mixture via latent B", run_jags("
model{
  for(i in 1:n){
    Y[i] ~ dnorm(B[i] * theta, 1)
    B[i] ~ dbern(0.5)
  }
  theta ~ dnorm(0, 0.01)
}", list(Y = Y15, n = 30), "theta", inits = function() list(theta = mean(Y15))))


# ---- Exercise 17 ---------------------------------------------------------------

report("Ex 17  T. rex growth curve", run_jags("
model{
  for(i in 1:n){
    mass[i] ~ dnorm(mu[i], tau)
    mu[i] <- theta1 + theta2 * pow(age[i], theta3)
  }
  theta1 ~ dnorm(0, 1 / 100^2)
  theta2 ~ dunif(0, 20000)
  theta3 ~ dnorm(0, 1)
  tau    ~ dgamma(0.01, 0.01)
  sigma <- 1 / sqrt(tau)
}", list(mass = c(29.9, 1761, 1807, 2984, 3230, 5040, 5654), age = c(2, 15, 14, 16, 18, 22, 28), n = 7),
  c("theta1", "theta2", "theta3", "sigma"),
  inits = list(list(theta1 = 0, theta2 = 200, theta3 = 1, tau = 1e-6),
               list(theta1 = 0, theta2 = 1000, theta3 = 0.5, tau = 1e-6),
               list(theta1 = 0, theta2 = 20, theta3 = 2, tau = 1e-6)), iters = 50000))


# ---- Exercise 18 ---------------------------------------------------------------

for (ab in c(1, 10)) report(sprintf("Ex 18  unknown n, a = b = %d", ab), run_jags("
model{
  Y ~ dbin(p, n)
  n ~ dpois(10)
  p ~ dbeta(a, b)
  theta <- n * p
}", list(Y = 10, a = ab, b = ab), c("n", "p", "theta"),
  inits = list(list(n = 10, p = 0.95), list(n = 30, p = 0.3), list(n = 15, p = 0.6))))


# ---- Exercise 19 ---------------------------------------------------------------

cat("\n==== Ex 19  effective sample sizes ====\n")
for (n in c(5, 25)) for (eps in c(0.1, 10)) {
  s <- run_jags("
model{
  for(i in 1:n){
    Y[i] ~ dnorm(mu[i], prec[i])
    mu[i] ~ dnorm(0, theta1)
    prec[i] ~ dgamma(theta2, theta3)
  }
  theta1 ~ dgamma(eps, eps)
  theta2 ~ dgamma(eps, eps)
  theta3 ~ dgamma(eps, eps)
}", list(Y = 1:n, n = n, eps = eps), c("theta1", "theta2", "theta3"), chains = 1)
  cat(sprintf("n = %2d  eps = %4.1f   ESS: %s\n", n, eps,
              paste(sprintf("%s %.0f", c("theta1", "theta2", "theta3"), effectiveSize(s)), collapse = "   ")))
}
