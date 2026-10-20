# =============================================================================
#  Chapter 7 -- JAGS versions of Exercises 3(b), 4 and 6
#  Run by GitHub Actions, not on the laptop. Data files are written by ch07.R.
#  Run from the repository root:   Rscript ch07/code/ch07_jags.R
# =============================================================================

suppressPackageStartupMessages({ library(rjags); library(coda) })
set.seed(2026)
DIR <- file.path("ch07", "code")

run_jags <- function(model, data, params, inits = NULL, chains = 3, burn = 5000, iters = 25000) {
  m <- jags.model(textConnection(model), data = data, inits = inits, n.chains = chains, quiet = TRUE)
  update(m, burn, progress.bar = "none")
  coda.samples(m, params, n.iter = iters, progress.bar = "none")
}
report <- function(title, samples, keep = NULL) {
  cat("\n==== ", title, " ====\n", sep = "")
  draws <- as.matrix(samples)
  if (!is.null(keep)) { draws <- draws[, keep, drop = FALSE]; samples <- samples[, keep, drop = FALSE] }
  tab <- t(apply(draws, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975)))))
  tab <- cbind(tab, ESS = effectiveSize(samples), Rhat = gelman.diag(samples, multivariate = FALSE)$psrf[, 1])
  print(signif(tab, 4))
  invisible(draws)
}


# ---- Exercise 3(b): bivariate random-effects meta-analysis ---------------------------

coch <- read.csv(file.path(DIR, "cochrane.csv"))
d3 <- report("Ex 3(b)  Cochrane meta-analysis", run_jags("
model{
  for(i in 1:K){
    ev_ctrl[i] ~ dbin(ilogit(alpha[i, 1]), n_ctrl[i])
    ev_trt[i]  ~ dbin(ilogit(alpha[i, 2]), n_trt[i])
    alpha[i, 1:2] ~ dmnorm(mu[1:2], Omega[1:2, 1:2])
  }
  for(k in 1:2){ mu[k] ~ dnorm(0, 1 / 100) }
  Omega[1:2, 1:2] ~ dwish(I2[1:2, 1:2], 3)        # Sigma = Omega^-1 ~ InvWishart(3, I)
  Sigma[1:2, 1:2] <- inverse(Omega[1:2, 1:2])
  lor <- mu[2] - mu[1]
}", list(ev_ctrl = coch$ev_ctrl, n_ctrl = coch$n_ctrl, ev_trt = coch$ev_trt, n_trt = coch$n_trt,
         K = nrow(coch), I2 = diag(2)), c("mu", "lor", "Sigma"), iters = 40000))
cat("P(lor < 0) =", mean(d3[, "lor"] < 0), "\n")


# ---- Exercise 4(c): hierarchical AR(1) pacing model, 3 chains x 25,000 ---------------

Y4 <- as.matrix(read.csv(file.path(DIR, "marathon_female_complete.csv")))
report("Ex 4(c)  marathon pacing model", run_jags("
model{
  for(i in 1:n){
    Y[i, 1] ~ dnorm(mu[i], 1 / sigma0^2)
    for(j in 2:J){
      Y[i, j] ~ dnorm(mu[i] + rho[i] * (Y[i, j-1] - mu[i]), prec[i])
    }
    mu[i]   ~ dnorm(theta1, 1 / theta2)
    rho[i]  ~ dnorm(theta3, 1 / theta4)
    prec[i] ~ dgamma(theta5, theta6)               # sigma_i^2 ~ InvGamma(theta5, theta6)
  }
  theta1 ~ dnorm(0, 1 / 100^2)
  theta2 <- 1 / itheta2;  itheta2 ~ dgamma(0.1, 0.1)
  theta3 ~ dnorm(0, 1 / 10^2)
  theta4 <- 1 / itheta4;  itheta4 ~ dgamma(0.1, 0.1)
  theta5 ~ dgamma(0.1, 0.1)
  theta6 ~ dgamma(0.1, 0.1)
  sigma0 <- 1 / sqrt(isigma0);  isigma0 ~ dgamma(0.1, 0.1)
}", list(Y = Y4, n = nrow(Y4), J = ncol(Y4)), c(paste0("theta", 1:6), "sigma0"),
  inits = list(list(theta5 = 2, theta6 = 0.1), list(theta5 = 0.5, theta6 = 0.01), list(theta5 = 5, theta6 = 1))))


# ---- Exercise 6(b): spatial model for village malaria prevalence --------------------

vill <- read.csv(file.path(DIR, "gambia_villages.csv"))
D6 <- as.matrix(dist(vill[, c("x", "y")]))
d6 <- report("Ex 6(b)  spatial malaria model (short run)", run_jags("
model{
  for(v in 1:V){
    Y[v] ~ dbin(ilogit(alpha[v]), n[v])
    m[v] <- mu
    for(u in 1:V){ Cov[v, u] <- s2 * exp(-D[v, u] / rho) }
  }
  alpha[1:V] ~ dmnorm.vcov(m[1:V], Cov[1:V, 1:V])
  mu  ~ dnorm(0, 1 / 100)
  s2  <- 1 / is2;  is2 ~ dgamma(0.1, 0.1)
  sigma <- sqrt(s2)
  rho ~ dunif(0, dstar)
}", list(Y = vill$Y, n = vill$n, V = nrow(vill), D = D6, dstar = max(D6)),
  c("mu", "sigma", "rho"), burn = 1000, iters = 4000))     # JAGS is slow here: 65 x 65 inverse per update
