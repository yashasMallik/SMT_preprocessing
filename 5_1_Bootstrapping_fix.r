# =============================================================================
# Bootstrap Confidence Intervals for Two-Component Residence Time CDF
# and Statistical Comparison Between Two Conditions
#
# Model: P(t) = alpha * exp(-ks * t) + (1 - alpha) * exp(-kns * t)
#
# Usage:
#   1. Set the file paths and column names in the CONFIGURATION section
#   2. Adjust starting parameter guesses if your data differs substantially
#      from the defaults (ks~0.18, kns~0.73, alpha~0.26)
#   3. Run the full script — outputs are saved to your working directory
#
# Required packages: minpack.lm, ggplot2, dplyr, tidyr
#   Install with: install.packages(c("minpack.lm", "ggplot2", "dplyr", "tidyr"))
# =============================================================================

library(minpack.lm)
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# CONFIGURATION — edit this section
# =============================================================================

# Paths to your DiaTrack CSV output files (one per condition)
# Expected columns: trajectory number, frame, x, y, z
# Column names are detected automatically — see COLUMN NAMES below
FILE_A <- "2-12-26- SRCAP HALO - 25PM.csv"   # e.g. DMSO / untreated
FILE_B <- "4-20-26- NLS SNAP - 25PM.csv"   # e.g. dBET-1 / MZ-1 treated

# --- Column names (from DiaTrack CSV output) ---
# DiaTrack exports: unnamed row index, Trajectory, Frame, x, y, z
# These match the exact column headers in the DiaTrack CSV.
COL_TRAJECTORY <- "Trajectory"
COL_FRAME      <- "Frame"

# Frame interval in seconds (i.e. time between frames)
# e.g. 10 ms exposure → FRAME_INTERVAL <- 0.01
FRAME_INTERVAL <- 0.2   # seconds per frame — CHANGE THIS to match your acquisition

# Minimum trajectory length to include (frames).
# Trajectories of 1 frame have undefined dwell time — exclude them.
MIN_TRAJ_FRAMES <- 2

# Labels for the two conditions (used in plots and output tables)
LABEL_A <- "Condition A"
LABEL_B <- "Condition B"

# Number of bootstrap replicates (1000 is fast; 5000 is more precise)
B <- 2000

# Maximum time (s) to plot on the x-axis
T_MAX <- 100

# Output directory (set to "." for current working directory)
OUT_DIR <- "."

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load a DiaTrack CSV and compute per-trajectory dwell times (in seconds).
#' DiaTrack format: unnamed row index, Trajectory, Frame, x, y, z
load_times <- function(path,
                       col_traj       = COL_TRAJECTORY,
                       col_frame      = COL_FRAME,
                       frame_interval = FRAME_INTERVAL,
                       min_frames     = MIN_TRAJ_FRAMES) {

  df <- read.csv(path, stringsAsFactors = FALSE, row.names = 1)
  message(sprintf("Read %d rows from %s", nrow(df), path))

  # Validate columns exist
  for (col in c(col_traj, col_frame)) {
    if (!col %in% names(df)) {
      stop(sprintf("Column '%s' not found. Available: %s",
                   col, paste(names(df), collapse = ", ")))
    }
  }

  traj_ids <- df[[col_traj]]
  frames   <- as.numeric(df[[col_frame]])

  # Dwell time = (max_frame - min_frame + 1) * frame_interval per trajectory
  traj_lengths <- tapply(frames, traj_ids, function(f) max(f) - min(f) + 1)
  traj_lengths <- traj_lengths[traj_lengths >= min_frames]

  dwell_times <- as.numeric(traj_lengths) * frame_interval
  dwell_times <- dwell_times[is.finite(dwell_times) & dwell_times > 0]

  message(sprintf("  Trajectories total:  %d", length(unique(traj_ids))))
  message(sprintf("  Trajectories kept:   %d  (>= %d frames)", length(dwell_times), min_frames))
  message(sprintf("  Frame interval:      %.4f s", frame_interval))
  message(sprintf("  Dwell time range:    %.3f – %.3f s\n", min(dwell_times), max(dwell_times)))

  dwell_times
}

#' Build a survival (1-CDF) data frame from a vector of dwell times
make_survival <- function(times) {
  t_sorted <- sort(times)
  n        <- length(t_sorted)
  surv     <- 1 - (seq_len(n) / n)
  data.frame(t = t_sorted, surv = surv)
}

#' Find the best starting parameters for a survival curve using a broad grid
#' search. Returns a list of the top N candidate start sets ranked by RSS,
#' so that fit_two_component can try each in turn until one converges.
candidate_starts <- function(surv_df, n_candidates = 10) {
  t    <- surv_df$t
  surv <- surv_df$surv

  # Broad grid — covers fast/slow populations across a wide range of rates
  alpha_grid <- seq(0.05, 0.95, by = 0.10)
  ks_grid    <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
  kns_grid   <- c(0.05, 0.1,  0.3, 0.7, 1.5, 3.0, 5.0, 10.0)

  results <- list()
  for (a in alpha_grid) {
    for (ks in ks_grid) {
      for (kns in kns_grid) {
        # Enforce the two components to be meaningfully separated
        if (kns <= ks * 1.5) next
        pred <- a * exp(-ks * t) + (1 - a) * exp(-kns * t)
        rss  <- sum((surv - pred)^2)
        if (is.finite(rss)) {
          results[[length(results) + 1]] <- list(
            alpha = a, ks = ks, kns = kns, rss = rss
          )
        }
      }
    }
  }

  # Sort by RSS and return top N
  rss_vals <- sapply(results, `[[`, "rss")
  top      <- results[order(rss_vals)[seq_len(min(n_candidates, length(results)))]]
  top
}

#' Fit the two-component model to a survival data frame.
#' Tries multiple candidate starting points and returns the best converging fit.
#' Returns a named numeric vector c(alpha, ks, kns) or NAs if all fail.
fit_two_component <- function(surv_df, starts = NULL) {
  candidates <- if (!is.null(starts)) list(starts) else candidate_starts(surv_df)

  best_rss <- Inf
  best_fit <- c(alpha = NA_real_, ks = NA_real_, kns = NA_real_)

  for (s in candidates) {
    result <- tryCatch({
      fit <- nlsLM(
        surv ~ alpha * exp(-ks * t) + (1 - alpha) * exp(-kns * t),
        data    = surv_df,
        start   = list(alpha = s$alpha, ks = s$ks, kns = s$kns),
        lower   = c(alpha = 0,   ks = 1e-6, kns = 1e-6),
        upper   = c(alpha = 1,   ks = Inf,  kns = Inf),
        control = nls.lm.control(maxiter = 500)
      )
      list(params = coef(fit), rss = sum(residuals(fit)^2))
    }, error = function(e) NULL)

    if (!is.null(result) && result$rss < best_rss) {
      best_rss <- result$rss
      best_fit <- result$params
    }
  }
  best_fit
}

#' Evaluate the two-component model at a vector of time points
predict_two_comp <- function(t, alpha, ks, kns) {
  alpha * exp(-ks * t) + (1 - alpha) * exp(-kns * t)
}

#' Run B bootstrap replicates for a vector of dwell times.
#' Candidate starts are found once from the full observed survival curve,
#' then reused across replicates for speed. Each replicate still tries all
#' candidates so convergence is robust even for unusual resamples.
#' Returns a B x 3 matrix of fitted parameters.
bootstrap_fit <- function(times, B = 2000) {
  n <- length(times)
  result <- matrix(NA_real_, nrow = B, ncol = 3,
                   dimnames = list(NULL, c("alpha", "ks", "kns")))

  # Find candidate starts from the full observed data once
  obs_surv   <- make_survival(times)
  candidates <- candidate_starts(obs_surv)
  best_start <- candidates[[1]]
  message(sprintf("  Best grid start: alpha=%.3f  ks=%.4f  kns=%.4f",
                  best_start$alpha, best_start$ks, best_start$kns))

  message(sprintf("Running %d bootstrap replicates...", B))
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  for (i in seq_len(B)) {
    resampled   <- sample(times, n, replace = TRUE)
    surv_df     <- make_survival(resampled)
    params      <- fit_two_component(surv_df, starts = NULL)
    result[i, ] <- params
    setTxtProgressBar(pb, i)
  }
  close(pb)

  n_failed <- sum(is.na(result[, "alpha"]))
  if (n_failed > 0) {
    message(sprintf("  Note: %d / %d replicates failed to converge (%.1f%%) — excluded from CIs.",
                    n_failed, B, 100 * n_failed / B))
  }
  result
}

#' Compute bootstrap CI envelope for the survival curve over a time grid
bootstrap_curve_ci <- function(boot_params, t_grid, level = 0.95) {
  probs   <- c((1 - level) / 2, 1 - (1 - level) / 2)
  valid   <- complete.cases(boot_params)
  params  <- boot_params[valid, , drop = FALSE]

  curves <- apply(params, 1, function(p) {
    predict_two_comp(t_grid, p["alpha"], p["ks"], p["kns"])
  })
  # curves is length(t_grid) x n_valid
  lo <- apply(curves, 1, quantile, probs = probs[1])
  hi <- apply(curves, 1, quantile, probs = probs[2])
  med <- apply(curves, 1, median)
  data.frame(t = t_grid, lo = lo, hi = hi, median = med)
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

# --- Load data ----------------------------------------------------------------
times_A <- load_times(FILE_A)
times_B <- load_times(FILE_B)

# --- Fit observed data --------------------------------------------------------
surv_A      <- make_survival(times_A)
surv_B      <- make_survival(times_B)
obs_fit_A   <- fit_two_component(surv_A)
obs_fit_B   <- fit_two_component(surv_B)

cat("\n--- Observed parameter estimates ---\n")
cat(sprintf("%-20s  alpha = %.4f  ks = %.4f  kns = %.4f\n",
            LABEL_A, obs_fit_A["alpha"], obs_fit_A["ks"], obs_fit_A["kns"]))
cat(sprintf("%-20s  alpha = %.4f  ks = %.4f  kns = %.4f\n\n",
            LABEL_B, obs_fit_B["alpha"], obs_fit_B["ks"], obs_fit_B["kns"]))

# Abort early with a clear message if either observed fit failed
for (lbl_fit in list(list(LABEL_A, obs_fit_A), list(LABEL_B, obs_fit_B))) {
  if (any(is.na(lbl_fit[[2]]))) {
    stop(sprintf(paste0(
      "\nObserved fit FAILED for condition '%s'.\n",
      "This usually means the two-component model does not fit this CDF well.\n",
      "Suggestions:\n",
      "  1. Plot the raw survival curve to inspect its shape.\n",
      "  2. The data may be better described by a single exponential.\n",
      "  3. Try adjusting MIN_TRAJ_FRAMES or FRAME_INTERVAL.\n",
      "  4. Check that FILE_A / FILE_B point to the correct condition."
    ), lbl_fit[[1]]))
  }
}

# --- Bootstrap ----------------------------------------------------------------
set.seed(42)
boot_A <- bootstrap_fit(times_A, B)
boot_B <- bootstrap_fit(times_B, B)

# --- Parameter CIs ------------------------------------------------------------
ci_A <- apply(boot_A, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_B <- apply(boot_B, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

cat("--- 95% Bootstrap CIs ---\n")
for (param in c("alpha", "ks", "kns")) {
  cat(sprintf("%s  [%s]: obs = %.4f  95%% CI [%.4f, %.4f]\n",
              param, LABEL_A, obs_fit_A[param], ci_A["2.5%", param], ci_A["97.5%", param]))
  cat(sprintf("%s  [%s]: obs = %.4f  95%% CI [%.4f, %.4f]\n\n",
              param, LABEL_B, obs_fit_B[param], ci_B["2.5%", param], ci_B["97.5%", param]))
}

# --- Between-condition comparison ---------------------------------------------
cat("--- Between-condition comparison (A - B) ---\n")
param_summary <- data.frame(
  parameter  = c("alpha", "ks", "kns"),
  obs_A      = obs_fit_A[c("alpha","ks","kns")],
  obs_B      = obs_fit_B[c("alpha","ks","kns")],
  obs_diff   = obs_fit_A[c("alpha","ks","kns")] - obs_fit_B[c("alpha","ks","kns")]
)

for (param in c("alpha", "ks", "kns")) {
  delta      <- boot_A[, param] - boot_B[, param]
  ci_delta   <- quantile(delta, c(0.025, 0.975), na.rm = TRUE)
  obs_d      <- obs_fit_A[param] - obs_fit_B[param]
  # Two-sided bootstrap p-value
  p_val      <- 2 * min(mean(delta < 0, na.rm = TRUE),
                        mean(delta > 0, na.rm = TRUE))
  p_val      <- max(p_val, 1 / sum(complete.cases(boot_A)))  # floor at 1/B

  cat(sprintf("%s: obs diff = %+.4f  95%% CI [%+.4f, %+.4f]  p = %.4f  %s\n",
              param, obs_d, ci_delta[1], ci_delta[2], p_val,
              ifelse(ci_delta[1] > 0 | ci_delta[2] < 0, "*", "ns")))
}

# =============================================================================
# PLOTS
# =============================================================================

t_grid <- seq(0, T_MAX, length.out = 500)

# --- Survival curve CI envelopes ----------------------------------------------
env_A <- bootstrap_curve_ci(boot_A, t_grid)
env_B <- bootstrap_curve_ci(boot_B, t_grid)
env_A$condition <- LABEL_A
env_B$condition <- LABEL_B
env_all <- rbind(env_A, env_B)

# Observed raw survival for plotting
obs_plot_A <- make_survival(times_A) %>% mutate(condition = LABEL_A)
obs_plot_B <- make_survival(times_B) %>% mutate(condition = LABEL_B)
obs_plot   <- rbind(obs_plot_A, obs_plot_B)

# Observed fitted curves
obs_curve_A <- data.frame(
  t = t_grid,
  surv = predict_two_comp(t_grid, obs_fit_A["alpha"], obs_fit_A["ks"], obs_fit_A["kns"]),
  condition = LABEL_A
)
obs_curve_B <- data.frame(
  t = t_grid,
  surv = predict_two_comp(t_grid, obs_fit_B["alpha"], obs_fit_B["ks"], obs_fit_B["kns"]),
  condition = LABEL_B
)
obs_curves <- rbind(obs_curve_A, obs_curve_B)

cols <- c("#378ADD", "#D85A30")
names(cols) <- c(LABEL_A, LABEL_B)

p1 <- ggplot() +
  # CI ribbon
  geom_ribbon(data = env_all,
              aes(x = t, ymin = lo, ymax = hi, fill = condition),
              alpha = 0.20) +
  # Raw survival (subsampled for speed)
  geom_point(data = obs_plot %>% group_by(condition) %>%
               slice(round(seq(1, n(), length.out = 300))),
             aes(x = t, y = surv, colour = condition),
             size = 0.6, alpha = 0.5) +
  # Fitted curve
  geom_line(data = obs_curves,
            aes(x = t, y = surv, colour = condition),
            linewidth = 0.9) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, T_MAX), expand = c(0, 0)) +
  labs(
    title    = "Residence time survival curves with 95% bootstrap CI",
    x        = "Time (s)",
    y        = "Uncorrected survival probability (1-CDF)",
    colour   = NULL, fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.80, 0.80))

ggsave(file.path(OUT_DIR, "survival_curves_bootstrap_CI.pdf"),
       p1, width = 7, height = 5)
message("Saved: survival_curves_bootstrap_CI.pdf")

# --- Bootstrap parameter distributions with CIs ------------------------------
boot_long <- rbind(
  as.data.frame(boot_A) %>% mutate(condition = LABEL_A),
  as.data.frame(boot_B) %>% mutate(condition = LABEL_B)
) %>%
  pivot_longer(cols = c("alpha", "ks", "kns"),
               names_to = "parameter", values_to = "value") %>%
  filter(!is.na(value))

# Observed values for annotation
obs_vals <- rbind(
  data.frame(condition = LABEL_A, parameter = names(obs_fit_A), value = obs_fit_A),
  data.frame(condition = LABEL_B, parameter = names(obs_fit_B), value = obs_fit_B)
) %>% filter(parameter %in% c("alpha", "ks", "kns"))

p2 <- ggplot(boot_long, aes(x = value, fill = condition, colour = condition)) +
  geom_density(alpha = 0.30, linewidth = 0.7) +
  geom_vline(data = obs_vals,
             aes(xintercept = value, colour = condition),
             linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ parameter, scales = "free",
             labeller = as_labeller(c(alpha = "alpha (transient fraction)",
                                      ks    = "k[s] (s⁻¹)",
                                      kns   = "k[ns] (s⁻¹)"))) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  labs(
    title  = "Bootstrap distributions of fitted parameters",
    x      = "Parameter value",
    y      = "Density",
    fill   = NULL, colour = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave(file.path(OUT_DIR, "bootstrap_parameter_distributions.pdf"),
       p2, width = 9, height = 4)
message("Saved: bootstrap_parameter_distributions.pdf")

# --- Difference distributions (A - B) ----------------------------------------
diff_df <- data.frame(
  alpha = boot_A[, "alpha"] - boot_B[, "alpha"],
  ks    = boot_A[, "ks"]    - boot_B[, "ks"],
  kns   = boot_A[, "kns"]   - boot_B[, "kns"]
) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "diff") %>%
  filter(!is.na(diff))

obs_diff_df <- data.frame(
  parameter = c("alpha", "ks", "kns"),
  diff      = c(obs_fit_A["alpha"] - obs_fit_B["alpha"],
                obs_fit_A["ks"]    - obs_fit_B["ks"],
                obs_fit_A["kns"]   - obs_fit_B["kns"])
)

p3 <- ggplot(diff_df, aes(x = diff)) +
  geom_density(fill = "#7A5AF8", colour = "#7A5AF8", alpha = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid",  colour = "black",  linewidth = 0.6) +
  geom_vline(data = obs_diff_df,
             aes(xintercept = diff),
             linetype = "dashed", colour = "#D85A30", linewidth = 0.9) +
  facet_wrap(~ parameter, scales = "free",
             labeller = as_labeller(c(alpha = "Δ alpha",
                                      ks    = "Δ k[s] (s⁻¹)",
                                      kns   = "Δ k[ns] (s⁻¹)"))) +
  labs(
    title    = paste("Bootstrap difference distributions:", LABEL_A, "−", LABEL_B),
    subtitle = "Dashed red = observed difference | Solid black = zero (no difference)",
    x        = "Difference",
    y        = "Density"
  ) +
  theme_classic(base_size = 13)

ggsave(file.path(OUT_DIR, "bootstrap_difference_distributions.pdf"),
       p3, width = 9, height = 4)
message("Saved: bootstrap_difference_distributions.pdf")

# --- Save numeric summary to CSV ---------------------------------------------
summary_rows <- list()
for (param in c("alpha", "ks", "kns")) {
  for (cond in list(list(label = LABEL_A, obs = obs_fit_A, boot = boot_A, ci = ci_A),
                    list(label = LABEL_B, obs = obs_fit_B, boot = boot_B, ci = ci_B))) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      condition  = cond$label,
      parameter  = param,
      observed   = cond$obs[param],
      ci_lo_95   = cond$ci["2.5%",  param],
      ci_hi_95   = cond$ci["97.5%", param],
      boot_mean  = mean(cond$boot[, param], na.rm = TRUE),
      boot_sd    = sd(cond$boot[,   param], na.rm = TRUE)
    )
  }
}
summary_df <- do.call(rbind, summary_rows)

# Add difference rows
for (param in c("alpha", "ks", "kns")) {
  delta    <- boot_A[, param] - boot_B[, param]
  ci_delta <- quantile(delta, c(0.025, 0.975), na.rm = TRUE)
  p_val    <- 2 * min(mean(delta < 0, na.rm = TRUE), mean(delta > 0, na.rm = TRUE))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    condition  = paste(LABEL_A, "-", LABEL_B),
    parameter  = param,
    observed   = obs_fit_A[param] - obs_fit_B[param],
    ci_lo_95   = ci_delta[1],
    ci_hi_95   = ci_delta[2],
    boot_mean  = mean(delta, na.rm = TRUE),
    boot_sd    = sd(delta,   na.rm = TRUE)
  )
}

final_summary <- do.call(rbind, summary_rows)
write.csv(final_summary, file.path(OUT_DIR, "bootstrap_summary.csv"), row.names = FALSE)
message("Saved: bootstrap_summary.csv")

cat("\n=== Analysis complete ===\n")
cat("Output files:\n")
cat("  survival_curves_bootstrap_CI.pdf\n")
cat("  bootstrap_parameter_distributions.pdf\n")
cat("  bootstrap_difference_distributions.pdf\n")
cat("  bootstrap_summary.csv\n")
