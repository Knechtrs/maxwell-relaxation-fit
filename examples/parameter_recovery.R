# Parameter recovery and model-selection checks on synthetic data.
# Run from the repository root:  Rscript examples/parameter_recovery.R
# Produces man/figures/recovery.png and prints summary tables.

source("R/maxwell.R")

t <- seq(0.5, 3000, by = 0.5)                # 2 Hz for 3000 s, as in the experiments
A_true <- c(0.3, 0.5)
tau_true <- c(20, 800)
truth <- maxwell_model(t, A_true, tau_true)
n_rep <- 30

# 1. Recovery with independent noise ------------------------------------------
est <- t(sapply(seq_len(n_rep), function(s) {
  set.seed(s)
  y <- truth + rnorm(length(t), 0, 0.003)
  f <- fit_maxwell(t, y, n_elements = 2, n_starts = 10, seed = s)
  c(f$A, f$tau)
}))
colnames(est) <- c("A1", "A2", "tau1", "tau2")
cat("\nRecovered parameters, independent noise (true: A = 0.3, 0.5; tau = 20, 800 s)\n")
print(round(rbind(mean = colMeans(est), sd = apply(est, 2, sd)), 3))

# 2. Model selection with independent vs autocorrelated noise -----------------
# Filtered instrument signals have strongly autocorrelated noise. Information
# criteria assume independent errors, so they over-reward extra elements.
select_n <- function(noise_fun, thin) {
  sapply(seq_len(n_rep), function(s) {
    set.seed(s)
    y <- truth + noise_fun()
    idx <- if (thin) log_thin(t, 200) else seq_along(t)
    cmp <- compare_maxwell(t[idx], y[idx], max_elements = 3, n_starts = 10, seed = s)
    cmp$n_elements[which.min(cmp$bic)]
  })
}
iid <- function() rnorm(length(t), 0, 0.003)
ar1 <- function() as.numeric(arima.sim(list(ar = 0.995), n = length(t), sd = 0.0003))

choices <- list(
  "Independent\nall points"   = select_n(iid, FALSE),
  "Correlated\nall points"    = select_n(ar1, FALSE),
  "Correlated\nlog-thinned"   = select_n(ar1, TRUE)
)
tab <- sapply(choices, function(x) table(factor(x, levels = 1:3)) / length(x))
cat("\nFraction of replicates where BIC selects 1, 2 or 3 elements (truth = 2)\n")
print(round(tab, 2))

# 3. Figure --------------------------------------------------------------------
set.seed(1)
y_ex <- truth + rnorm(length(t), 0, 0.003)
f1 <- fit_maxwell(t, y_ex, 1, seed = 1)
f2 <- fit_maxwell(t, y_ex, 2, seed = 1)
show <- log_thin(t, 400)

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
png("man/figures/recovery.png", width = 2100, height = 680, res = 200)
par(mfrow = c(1, 3), mar = c(4.4, 4.2, 2.6, 1), cex = 0.85, font.main = 1)

plot(t[show], y_ex[show], log = "x", pch = 16, cex = 0.35, col = "grey55",
     xlab = "Time (s, log scale)", ylab = "Normalized stress",
     main = "A  Synthetic relaxation curve")
lines(t, f1$fitted, col = "#d95f02", lwd = 1.5)
lines(t, f2$fitted, col = "#1b9e77", lwd = 1.5, lty = 2)
points(t[show], y_ex[show], pch = 16, cex = 0.3, col = "grey40")
legend("bottomleft", bty = "n", lwd = 1.5, lty = c(1, 2, NA), pch = c(NA, NA, 16),
       col = c("#d95f02", "#1b9e77", "grey40"),
       legend = c("1 element", "2 elements", "data"))

hist(est[, "tau2"], breaks = 12, col = "#1b9e77", border = "white",
     xlab = expression(hat(tau)[2] ~ "(s)"), main = expression("B  Recovery of" ~ tau[2]))
abline(v = tau_true[2], lwd = 2, lty = 2)
legend("topright", bty = "n", lty = 2, lwd = 2, legend = "true value")

barplot(tab, col = c("#d95f02", "#1b9e77", "#7570b3"), border = NA,
        ylab = "Fraction of replicates", main = "C  Elements chosen by BIC",
        names.arg = rep("", 3), ylim = c(0, 1.18), axes = FALSE)
axis(2, at = seq(0, 1, 0.2))
mtext(names(choices), side = 1, at = c(0.7, 1.9, 3.1), line = 1.6, cex = 0.62)
legend("top", horiz = TRUE, bty = "n", fill = c("#d95f02", "#1b9e77", "#7570b3"),
       border = NA, legend = c("1 element", "2 elements", "3 elements"), cex = 0.8)
dev.off()
cat("\nSaved man/figures/recovery.png\n")
