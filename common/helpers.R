# =============================================================================
#  Shared helpers for every chapter script (base R only).
#  Each chapter script calls  source("common/helpers.R")  and sets FIG_DIR.
# =============================================================================

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

# Equal-tailed 95% interval of a distribution given its quantile function.
ci95 <- function(qfun, ...) qfun(c(0.025, 0.975), ...)
fmt_ci <- function(ci, digits = 3) sprintf(paste0("(%.", digits, "f, %.", digits, "f)"), ci[1], ci[2])


# ---- MCMC diagnostics (no packages needed) -----------------------------------

# Effective sample size via Geyer's initial positive sequence of autocorrelations.
ess <- function(x) {
  x <- as.numeric(x)
  if (var(x) == 0) return(NA_real_)
  rho   <- acf(x, lag.max = min(length(x) - 1, 5000), plot = FALSE)$acf[, 1, 1]
  pairs <- rho[seq(1, length(rho) - 1, by = 2)] + rho[seq(2, length(rho), by = 2)]
  keep  <- if (any(pairs <= 0)) seq_len(which(pairs <= 0)[1] - 1) else seq_along(pairs)
  length(x) / max(1, -1 + 2 * sum(pairs[keep]))
}

# Gelman-Rubin potential scale reduction factor; `chains` is iterations x chains.
rhat <- function(chains) {
  n <- nrow(chains)
  W <- mean(apply(chains, 2, var))
  B <- n * var(colMeans(chains))
  sqrt(((n - 1) / n * W + B / n) / W)
}

# Overlaid trace plot of several chains (iterations x chains matrix).
trace_plot <- function(chains, ylab, main = "") {
  matplot(chains, type = "l", lty = 1, lwd = 0.6, col = adjustcolor(PALETTE, 0.75),
          xlab = "iteration", ylab = ylab, main = main)
}

# Posterior summary row: mean, sd and 95% interval of a vector of draws.
post_summary <- function(x, digits = 4) {
  round(c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))), digits)
}
