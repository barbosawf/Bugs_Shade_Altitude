# 03 - Sensitivity analysis: Monalonion with all sampling occasions ------------------
# Refits the Monalonion models including the first three occasions
# (days 0, 16 and 31), in which only one Monalonion individual was collected in
# 72 trees, and compares tests, pairwise comparisons and fitted curves with the
# main analysis (days 46-167). Models are the same as in 01_Monalonion.R.

source(here::here("R_Scripts", "functions.R"))

run_bootstrap <- TRUE # set FALSE to skip (slow: 1000 simulations per test)

bugs <- read_bugs()
start <- min(bugs$Date)

sets <- list(`Days 46-167 (main analysis)` = bugs |> dplyr::filter(Time >= 46),
             `All occasions (days 0-167)` = bugs)

days_main <- sort(unique(sets[[1]]$Time)) # sampling days of the main analysis
day_ref <- mean(sets[[1]]$Time) # day used for the overall comparisons

lrt <- list()
blocks <- list()
tukey <- list()
per_day <- list()
curves <- list()

# Each data set is analysed in turn ---------------------------------------------------
for (nm in names(sets)) {
  d <- sets[[nm]]
  time_center <- mean(d$Time)
  time_scale <- sd(d$Time)
  d <- d |> dplyr::mutate(sT = (Time - time_center) / time_scale)

  # Final model: (Time + Time^2) x Elevation + Elevation x Shading
  m_final <- glmmTMB::glmmTMB(Monalonion ~ (sT + I(sT^2)) * Altitude + Altitude * Condition +
                                (1 | Farm) + (1 | Plot) + (1 | Event),
                              family = glmmTMB::nbinom1(), data = d)
  # All two-way interactions (linear and quadratic time)
  m_2way <- glmmTMB::glmmTMB(Monalonion ~ (sT + I(sT^2)) * (Altitude + Condition) + Altitude:Condition +
                      (1 | Farm) + (1 | Plot) + (1 | Event),
                    family = glmmTMB::nbinom1(), data = d)
  # Full model: (Time + Time^2) x Elevation x Shading
  m_full <- glmmTMB::glmmTMB(Monalonion ~ (sT + I(sT^2)) * Altitude * Condition +
                  (1 | Farm) + (1 | Plot) + (1 | Event),
                family = glmmTMB::nbinom1(), data = d)
  # Models of steps 2 and 3 of the backward elimination in 01_Monalonion.R
  m_s2 <- update(m_full, . ~ . - I(sT^2):Altitude:Condition - I(sT^2):Condition)
  m_s3 <- update(m_s2, . ~ . - sT:Altitude:Condition)

  # Type II LRT of the final model and test of the three-way block
  m_noTA <- update(m_final, . ~ . - sT:Altitude)
  m_noQA <- update(m_final, . ~ . - I(sT^2):Altitude)
  m_noAC <- update(m_final, . ~ . - Altitude:Condition)
  m_noAlt <- update(m_final, . ~ . - sT:Altitude - I(sT^2):Altitude - Altitude:Condition)
  lrt[[nm]] <- dplyr::bind_rows(
    lrt_row("sT", update(m_noTA, . ~ . - sT), m_noTA),
    lrt_row("I(sT^2)", update(m_noQA, . ~ . - I(sT^2)), m_noQA),
    lrt_row("Altitude", update(m_noAlt, . ~ . - Altitude), m_noAlt),
    lrt_row("Condition", update(m_noAC, . ~ . - Condition), m_noAC),
    lrt_row("sT:Altitude", m_noTA, m_final),
    lrt_row("I(sT^2):Altitude", m_noQA, m_final),
    lrt_row("Altitude:Condition", m_noAC, m_final)
  ) |>
    dplyr::mutate(Model = "Final") |>
    dplyr::bind_rows(lrt_row("Three-way block (Time + Time^2) x Elevation x Shading", m_2way, m_full) |>
                dplyr::mutate(Model = "Full vs two-way")) |>
    dplyr::mutate(Data = nm, .before = 1)

  # Critical tests, with parametric bootstrap (same tests as 01_Monalonion.R)
  block_pairs <- list(
    `Three-way: Time(linear+quadratic) x Elevation x Shading (vs all two-way)` = list(m_2way, m_full),
    `Three-way: Time(linear) x Elevation x Shading (step 3)` = list(m_s3, m_s2),
    `Quadratic time x Elevation (final model)` = list(m_noQA, m_final)
  )
  blocks[[nm]] <- purrr::imap_dfr(block_pairs, function(p, tn) {
    if (run_bootstrap) {
      set.seed(2026) # simulations of the parametric bootstrap
      pb_lrt(p[[1]], p[[2]], d, "Monalonion") |> dplyr::mutate(Test = tn, .before = 1)
    } else {
      lrt_row(tn, p[[1]], p[[2]]) |> dplyr::rename(Test = Term, P_asymptotic = P)
    }
  }) |>
    dplyr::mutate(Data = nm, .before = 1)

  # Pairwise comparisons (final model), at the same day in both data sets
  sT_ref <- (day_ref - time_center) / time_scale
  tukey[[nm]] <- dplyr::bind_rows(
    cld_tbl(emmeans::emmeans(m_final, ~ Altitude | Condition, at = list(sT = sT_ref), type = "response")) |>
      dplyr::mutate(Comparison = "Elevation within shading"),
    emmeans::contrast(emmeans::emmeans(m_final, ~ Condition | Altitude, at = list(sT = sT_ref)), "pairwise") |>
      summary(type = "response") |>
      tibble::as_tibble() |>
      dplyr::mutate(Comparison = "Shaded vs unshaded within elevation")
  ) |>
    dplyr::mutate(Data = nm, .before = 1)

  # Elevation comparisons per sampling day (final model)
  per_day[[nm]] <- cld_tbl(emmeans::emmeans(m_final, ~ Altitude | sT,
                                   at = list(sT = (days_main - time_center) / time_scale),
                                   type = "response")) |>
    dplyr::mutate(Time = round(sT * time_scale + time_center), Data = nm)

  # Fitted curves (final model)
  grid_days <- seq(min(d$Time), 167, 1)
  curves[[nm]] <- emmeans::emmeans(m_final, ~ Altitude * Condition * sT,
                          at = list(sT = (grid_days - time_center) / time_scale),
                          type = "response") |>
    tibble::as_tibble() |>
    dplyr::mutate(Time = sT * time_scale + time_center, Data = nm) |>
    dplyr::rename(Predicted = response, Lower = asymp.LCL, Upper = asymp.UCL)
}

# Tables ----------------------------------------------------------------------------------

lrt <- dplyr::bind_rows(lrt)
readr::write_csv(lrt, file.path(dir_tables, "03_Mon_sensitivity_LRT.csv"))
print(lrt, n = Inf)

blocks <- dplyr::bind_rows(blocks)
readr::write_csv(blocks, file.path(dir_tables, "03_Mon_sensitivity_interaction_blocks.csv"))
print(blocks, width = Inf)

tukey <- dplyr::bind_rows(tukey)
readr::write_csv(tukey, file.path(dir_tables, "03_Mon_sensitivity_pairwise_final_model.csv"))

letters_cmp <- dplyr::bind_rows(per_day) |>
  dplyr::select(Data, Time, Altitude, .group) |>
  tidyr::pivot_wider(names_from = Data, values_from = .group)
readr::write_csv(letters_cmp, file.path(dir_tables, "03_Mon_sensitivity_letters_per_day_final_model.csv"))
print(letters_cmp, n = Inf)

curves <- dplyr::bind_rows(curves)
peaks <- curves |>
  dplyr::group_by(Data, Altitude, Condition) |>
  dplyr::slice_max(Predicted, n = 1) |>
  dplyr::ungroup() |>
  dplyr::select(Data, Altitude, Condition, Time, Predicted, Lower, Upper)
readr::write_csv(peaks, file.path(dir_tables, "03_Mon_sensitivity_peaks_final_model.csv"))
print(peaks, n = Inf)

# Figure S6 -------------------------------------------------------------------------------

obs_means <- bugs |>
  dplyr::group_by(Time, Altitude, Condition) |>
  dplyr::summarise(Mean = mean(Monalonion), SE = sd(Monalonion) / sqrt(dplyr::n()), .groups = "drop")

set_cols <- c(`Days 46-167 (main analysis)` = "#00AFBB", `All occasions (days 0-167)` = "#C0392B")
fig_s <- ggplot2::ggplot(curves, ggplot2::aes(Time, Predicted, colour = Data, fill = Data)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_pointrange(data = obs_means, ggplot2::aes(Time, Mean, ymin = Mean - SE, ymax = Mean + SE),
                  inherit.aes = FALSE, size = 0.2) +
  ggplot2::facet_grid(Condition ~ Altitude, scales = "free_y") +
  ggplot2::scale_colour_manual(values = set_cols) +
  ggplot2::scale_fill_manual(values = set_cols) +
  ggplot2::labs(x = "Time (days)", y = "Number of Monalonion per tree", colour = NULL, fill = NULL,
       caption = "Final model (NB1 GLMM). Points: observed mean (± SE) per day.") +
  theme_bugs()
save_fig(fig_s, "S6_Monalonion_sensitivity_all_occasions", 11, 6.5)
