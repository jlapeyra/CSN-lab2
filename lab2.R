# ============================================================
# Fit-and-plot degree distributions from *_degree_sequence.txt
# Models: Geometric (support 1..), Zero-truncated Poisson, Zeta (power law),
#         Altmann-Zeta (shifted zeta), Truncated Zeta (power law w/ cutoff)
# Output: ./plots/<Language>_pmf.png and ./plots/<Language>_ccdf.png
# ============================================================

# --- Packages ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gsl)        # for Riemann/Hurwitz zeta and polylog
})

# --- IO ---
DATA_DIR  <- "./data"
PLOT_DIR  <- "./plots"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Helpers ---
empirical_pmf <- function(x) {
  x <- x[x >= 1L]
  tab <- as.integer(table(factor(x, levels = sort(unique(x)))))
  k   <- sort(unique(x))
  data.table(k = k, pmf = tab / sum(tab))
}

empirical_ccdf <- function(x) {
  x <- x[x >= 1L]
  k  <- sort(unique(x))
  cc <- vapply(k, function(t) mean(x >= t), numeric(1))
  data.table(k = k, ccdf = cc)
}

# log-factorial via lgamma
lfact <- function(k) lgamma(k + 1)

# --- Model PMFs (all on support k = 1,2,...) ---
pmf_geometric <- function(k, p) { pmax(p * (1 - p)^(k - 1), .Machine$double.eps) }             # 0 < p < 1
pmf_ztpois    <- function(k, lambda) {                                                           # lambda > 0
  denom <- 1 - exp(-lambda)
  pmax(exp(-lambda) * lambda^k / exp(lfact(k)) / denom, .Machine$double.eps)
}
pmf_zeta <- function(k, alpha) {                                                                 # alpha > 1
  Z <- zeta(alpha)
  pmax(k^(-alpha) / Z, .Machine$double.eps)
}
pmf_altmann <- function(k, alpha, delta) {                                                       # alpha > 1, delta >= 0
  S <- hzeta(alpha, 1 + delta)   # Hurwitz zeta: sum_{k=0..∞} (k+1+delta)^(-alpha)
  pmax((k + delta)^(-alpha) / S, .Machine$double.eps)
}
pmf_trunczeta <- function(k, alpha, tau) {                                                       # alpha > 0, tau > 0
  # Power-law with exponential cutoff: ∝ k^{-alpha} * exp(-tau k)
  # Normalizer is polylog(alpha, e^{-tau}) = Li_{alpha}(exp(-tau))
  Z <- polylog(alpha, exp(-tau))
  pmax((k^(-alpha)) * exp(-tau * k) / Z, .Machine$double.eps)
}

# --- Negative log-likelihoods (for MLE) ---
nll_geom <- function(theta, x) {
  p <- plogis(theta[1])                    # map R -> (0,1)
  -sum(log(pmf_geometric(x, p)))
}
nll_ztpois <- function(theta, x) {
  lambda <- exp(theta[1])                  # map R -> (0,∞)
  -sum(log(pmf_ztpois(x, lambda)))
}
nll_zeta <- function(theta, x) {
  alpha <- 1 + exp(theta[1])               # map R -> (1,∞)
  -sum(log(pmf_zeta(x, alpha)))
}
nll_altmann <- function(theta, x) {
  alpha <- 1 + exp(theta[1])               # (1,∞)
  delta <- exp(theta[2]) - 1               # [0,∞)
  -sum(log(pmf_altmann(x, alpha, delta)))
}
nll_trunczeta <- function(theta, x) {
  alpha <- exp(theta[1]) + 1e-6            # (0,∞) avoid zero
  tau   <- exp(theta[2]) + 1e-6            # (0,∞)
  -sum(log(pmf_trunczeta(x, alpha, tau)))
}

# --- Fit wrappers (robust starts + optim) ---
fit_geometric <- function(x) {
  p0 <- 1 / mean(x)                        # MLE for geometric on {1..∞}
  o  <- optim(qlogis(p0), nll_geom, x = x, method = "BFGS")
  list(p = plogis(o$par[1]), nll = o$value, conv = o$convergence)
}
fit_ztpois <- function(x) {
  # Solve E[X | X>=1] = m = lambda / (1 - e^{-lambda}) for start
  m <- mean(x)
  fn <- function(l) l / (1 - exp(-l)) - m
  l0 <- try(uniroot(fn, c(1e-6, max(10, 2*m)))$root, silent = TRUE)
  l0 <- if (inherits(l0, "try-error")) m else l0
  o  <- optim(log(l0), nll_ztpois, x = x, method = "BFGS")
  list(lambda = exp(o$par[1]), nll = o$value, conv = o$convergence)
}
fit_zeta <- function(x) {
  # crude start via Hill estimator for alpha on tail
  x_tail <- x[x >= quantile(x, 0.9)]
  alpha0 <- if (length(unique(x_tail)) > 1) 1 + length(x_tail) / sum(log(x_tail / min(x_tail))) else 2
  o <- optim(log(alpha0 - 1), nll_zeta, x = x, method = "BFGS")
  list(alpha = 1 + exp(o$par[1]), nll = o$value, conv = o$convergence)
}
fit_altmann <- function(x) {
  # starts: alpha from zeta fit; delta near 0
  a0 <- fit_zeta(x)$alpha
  o  <- optim(c(log(a0 - 1), log(1 + 1e-3)), nll_altmann, x = x, method = "BFGS")
  list(alpha = 1 + exp(o$par[1]), delta = exp(o$par[2]) - 1, nll = o$value, conv = o$convergence)
}
fit_trunczeta <- function(x) {
  # starts: alpha from zeta, tau small
  a0 <- fit_zeta(x)$alpha
  o  <- optim(c(log(max(1e-3, a0)), log(1e-2)), nll_trunczeta, x = x, method = "BFGS")
  list(alpha = exp(o$par[1]) + 1e-6, tau = exp(o$par[2]) + 1e-6, nll = o$value, conv = o$convergence)
}

# --- Plotting ---
make_pmf_plot <- function(lang, x, fits) {
  emp <- empirical_pmf(x)
  kmax <- max(quantile(x, 0.99), 20)
  kseq <- 1:as.integer(kmax)
  
  curves <- rbind(
    data.table(k = kseq, y = pmf_geometric(kseq, fits$geom$p),        model = "Geometric"),
    data.table(k = kseq, y = pmf_ztpois(kseq,    fits$ztp$lambda),     model = "Poisson (ZTP)"),
    data.table(k = kseq, y = pmf_zeta(kseq,      fits$zeta$alpha),     model = "Zeta"),
    data.table(k = kseq, y = pmf_altmann(kseq,   fits$alt$alpha, fits$alt$delta), model = "Altmann-Zeta"),
    data.table(k = kseq, y = pmf_trunczeta(kseq, fits$tz$alpha, fits$tz$tau),     model = "Truncated Zeta")
  )
  
  ggplot() +
    geom_col(data = emp, aes(k, pmf), width = 0.9, alpha = 0.35) +
    geom_line(data = curves, aes(k, y, linetype = model)) +
    scale_y_continuous(trans = "log10") +
    coord_cartesian(xlim = c(1, kmax)) +
    labs(title = paste0(lang, " — PMF (log scale)"),
         x = "Degree k", y = "P(K = k)", linetype = "Model") +
    theme_minimal(base_size = 12)
}

make_ccdf_plot <- function(lang, x, fits) {
  emp <- empirical_ccdf(x)
  kmax <- max(quantile(x, 0.99), 20)
  kseq <- 1:as.integer(kmax)
  
  # compute CCDFs from PMFs
  ccdf_from_pmf <- function(p) rev(cumsum(rev(p)))
  cc_curves <- rbind(
    data.table(k = kseq, y = ccdf_from_pmf(pmf_geometric(kseq, fits$geom$p)),        model = "Geometric"),
    data.table(k = kseq, y = ccdf_from_pmf(pmf_ztpois(kseq,    fits$ztp$lambda)),     model = "Poisson (ZTP)"),
    data.table(k = kseq, y = ccdf_from_pmf(pmf_zeta(kseq,      fits$zeta$alpha)),     model = "Zeta"),
    data.table(k = kseq, y = ccdf_from_pmf(pmf_altmann(kseq,   fits$alt$alpha, fits$alt$delta)), model = "Altmann-Zeta"),
    data.table(k = kseq, y = ccdf_from_pmf(pmf_trunczeta(kseq, fits$tz$alpha, fits$tz$tau)),     model = "Truncated Zeta")
  )
  
  ggplot() +
    geom_point(data = emp, aes(k, ccdf), size = 1.2, alpha = 0.7) +
    geom_line(data = cc_curves, aes(k, y, linetype = model)) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    coord_cartesian(xlim = c(1, kmax)) +
    labs(title = paste0(lang, " — CCDF (log–log)"),
         x = "Degree k", y = "P(K ≥ k)", linetype = "Model") +
    theme_minimal(base_size = 12)
}

# --- Main loop: read, fit, plot ---
files <- list.files(DATA_DIR, pattern = "_degree_sequence\\.txt$", full.names = TRUE)
stopifnot(length(files) > 0)

for (f in files) {
  lang <- sub("_degree_sequence\\.txt$", "", basename(f))
  x    <- scan(f, what = integer(), quiet = TRUE)
  x    <- x[x >= 1L]
  
  if (!length(x)) {
    message("[skip] Empty sequence in: ", f)
    next
  }
  
  # Fit all models
  fits <- list(
    geom = fit_geometric(x),
    ztp  = fit_ztpois(x),
    zeta = fit_zeta(x),
    alt  = fit_altmann(x),
    tz   = fit_trunczeta(x)
  )
  
  # Save PMF plot
  g_pmf  <- make_pmf_plot(lang, x, fits)
  ggsave(file.path(PLOT_DIR, paste0(lang, "_pmf.png")), g_pmf, width = 7, height = 5, dpi = 150)
  
  # Save CCDF plot
  g_ccdf <- make_ccdf_plot(lang, x, fits)
  ggsave(file.path(PLOT_DIR, paste0(lang, "_ccdf.png")), g_ccdf, width = 7, height = 5, dpi = 150)
  
  message("[ok] Plots written for ", lang)
}

# --- (Optional) print a small fit summary table to console ---
summ <- rbindlist(lapply(seq_along(files), function(i) {
  f   <- files[i]
  lang <- sub("_degree_sequence\\.txt$", "", basename(f))
  x    <- scan(f, what = integer(), quiet = TRUE); x <- x[x >= 1L]
  if (!length(x)) return(NULL)
  g <- fit_geometric(x); p <- fit_ztpois(x); z <- fit_zeta(x); a <- fit_altmann(x); t <- fit_trunczeta(x)
  data.table(
    Language = lang,
    Geometric_NLL = g$nll,
    PoissonZTP_NLL = p$nll,
    Zeta_NLL = z$nll,
    AltmannZeta_NLL = a$nll,
    TruncZeta_NLL = t$nll
  )
}))
if (nrow(summ)) print(summ)

