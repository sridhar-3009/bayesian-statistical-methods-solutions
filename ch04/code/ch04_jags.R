# =============================================================================
#  Chapter 4 -- JAGS version of Exercise 9 (AR(1) model for "Sophia")
#  Run by GitHub Actions, not on the laptop.
#  Run from the repository root:   Rscript ch04/code/ch04_jags.R
# =============================================================================

suppressPackageStartupMessages({ library(rjags); library(coda) })
set.seed(2026)

dat <- read.csv(file.path("ch04", "code", "ex09_sophia_logodds.csv"))   # written by ch04.R

model <- "
model{
  for(t in 1:T){ mu[t] <- alpha + beta * t }
  for(t in 2:T){
    Y[t] ~ dnorm(mu[t] + rho * (Y[t-1] - mu[t-1]), prec)
  }
  alpha ~ dnorm(0, 1 / 100^2)
  beta  ~ dnorm(0, 1 / 100^2)
  rho   ~ dunif(-1, 1)
  prec  ~ dgamma(0.1, 0.1)          # sigma^2 ~ InvGamma(0.1, 0.1)
  sigma <- 1 / sqrt(prec)
}"
m <- jags.model(textConnection(model), data = list(Y = dat$Y, T = nrow(dat)),
                inits = list(list(rho = -0.5), list(rho = 0), list(rho = 0.9)),
                n.chains = 3, quiet = TRUE)
update(m, 10000, progress.bar = "none")
s <- coda.samples(m, c("alpha", "beta", "rho", "sigma"), n.iter = 30000, progress.bar = "none")

cat("\n==== Ex 9(b)  AR(1) model for Sophia ====\n")
draws <- as.matrix(s)
tab <- t(apply(draws, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975)))))
print(signif(cbind(tab, ESS = effectiveSize(s), Rhat = gelman.diag(s, multivariate = FALSE)$psrf[, 1]), 4))
