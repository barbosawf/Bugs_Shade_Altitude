# 01 - Monalonion ---------------------------------------------------------------
#
# Response: number of Monalonion individuals per tree (tree-level counts).
# Model: negative binomial GLMM (glmmTMB) with random intercepts for
#   Farm  - farms (two per elevation)
#   Plot  - farm x shading condition (group of trees observed on every occasion)
#   Event - farm x sampling date (both conditions sampled on the same day)
# Fixed effects are tested by Type II likelihood-ratio tests (each term is tested
# against the model with all terms that do not contain it). The final model was
# selected by backward elimination with LRTs (section 3), and the critical tests
# were checked by parametric bootstrap (section 4).
#
# Final model: (Time + Time^2) x Elevation + Elevation x Shading. The three-way
# interaction (Time + Time^2) x Elevation x Shading received only weak support
# and was not retained.
#
# Analyses use occasions from day 46 onwards (see 03_Sensitivity_Monalonion.R
# for the analysis with all occasions).

source(here::here("R_Scripts", "functions.R"))

run_bootstrap <- TRUE # parametric bootstrap (~2 min with 16 cores)
first_day <- 46

# Data -------------------------------------------------------------------------------

bugs <- read_bugs()
start <- min(bugs$Date)
dM <- bugs |> dplyr::filter(Time >= first_day)

# Time is standardized (sT) to improve the numerical stability of the GLMMs
time_center <- mean(dM$Time)
time_scale <- sd(dM$Time)
dM <- dM |> dplyr::mutate(sT = (Time - time_center) / time_scale)

# 1. Error distribution and random-effect structure ----------------------------------
# All candidates have the full fixed part (Time + Time^2) x Elevation x Shading.

m_pois_olre <- glmmTMB::glmmTMB(
  Monalonion ~ (sT + I(sT^2)) * Altitude * Condition +
    (1 | Farm) + (1 | Plot) + (1 | Event) + (1 | Tree),
  family = poisson(), data = dM)

m_nb1 <- glmmTMB::glmmTMB(
  Monalonion ~ (sT + I(sT^2)) * Altitude * Condition +
    (1 | Farm) + (1 | Plot) + (1 | Event),
  family = glmmTMB::nbinom1(), data = dM)

m_nb2 <- glmmTMB::glmmTMB(
  Monalonion ~ (sT + I(sT^2)) * Altitude * Condition +
    (1 | Farm) + (1 | Plot) + (1 | Event),
  family = glmmTMB::nbinom2(), data = dM)

m_nb1_glm <- glmmTMB::glmmTMB(
  Monalonion ~ (sT + I(sT^2)) * Altitude * Condition,
  family = glmmTMB::nbinom1(), data = dM)

m_nb1_ar1 <- glmmTMB::glmmTMB(
  Monalonion ~ (sT + I(sT^2)) * Altitude * Condition +
    (1 | Farm) + ar1(fTime + 0 | Plot),
  family = glmmTMB::nbinom1(), data = dM)

cands <- list(`Poisson GLMM + observation-level RE` = m_pois_olre,
              `NB1 GLMM` = m_nb1,
              `NB2 GLMM` = m_nb2,
              `NB1 GLM (no random effects)` = m_nb1_glm,
              `NB1 GLMM, AR(1) over dates within plot` = m_nb1_ar1)
tab_cands <- purrr::imap_dfr(cands, ~ tibble::tibble(Model = .y, AIC = AIC(.x),
                                      logLik = as.numeric(logLik(.x)),
                                      df = attr(logLik(.x), "df"))) |>
  dplyr::mutate(dAIC = AIC - min(AIC)) |>
  dplyr::arrange(AIC)
readr::write_csv(tab_cands, file.path(dir_tables, "01_Mon_candidate_models_AIC.csv"))
print(tab_cands)

# 2. Fixed-effect selection: Type II LRT of the full model ---------------------------

m_full <- m_nb1 # (Time + Time^2) x Elevation x Shading, NB1 GLMM

# Models without the three-way terms (needed as reference for the two-way terms)
m_no3lin <- update(m_full, . ~ . - sT:Altitude:Condition)
m_no3quad <- update(m_full, . ~ . - I(sT^2):Altitude:Condition)
m_no3 <- update(m_full, . ~ . - sT:Altitude:Condition - I(sT^2):Altitude:Condition)

# Models used to test the main effects
# (reduced model first, then the term is added at the end of the formula)
m_sT_0 <- update(m_full, . ~ . - sT - sT:Altitude - sT:Condition - sT:Altitude:Condition)
m_sT_1 <- update(m_sT_0, . ~ . + sT)
m_sT2_0 <- update(m_full, . ~ . - I(sT^2) - I(sT^2):Altitude - I(sT^2):Condition -
                    I(sT^2):Altitude:Condition)
m_sT2_1 <- update(m_sT2_0, . ~ . + I(sT^2))
m_alt_0 <- glmmTMB::glmmTMB(Monalonion ~ (sT + I(sT^2)) * Condition +
                     (1 | Farm) + (1 | Plot) + (1 | Event),
                   family = glmmTMB::nbinom1(), data = dM)
m_alt_1 <- update(m_alt_0, . ~ . + Altitude)
m_cond_0 <- glmmTMB::glmmTMB(Monalonion ~ (sT + I(sT^2)) * Altitude +
                      (1 | Farm) + (1 | Plot) + (1 | Event),
                    family = glmmTMB::nbinom1(), data = dM)
m_cond_1 <- update(m_cond_0, . ~ . + Condition)

lrt_full <- dplyr::bind_rows(
  lrt_row("sT", m_sT_0, m_sT_1),
  lrt_row("I(sT^2)", m_sT2_0, m_sT2_1),
  lrt_row("Altitude", m_alt_0, m_alt_1),
  lrt_row("Condition", m_cond_0, m_cond_1),
  lrt_row("sT:Altitude", update(m_no3lin, . ~ . - sT:Altitude), m_no3lin),
  lrt_row("I(sT^2):Altitude", update(m_no3quad, . ~ . - I(sT^2):Altitude), m_no3quad),
  lrt_row("sT:Condition", update(m_no3lin, . ~ . - sT:Condition), m_no3lin),
  lrt_row("I(sT^2):Condition", update(m_no3quad, . ~ . - I(sT^2):Condition), m_no3quad),
  lrt_row("Altitude:Condition", update(m_no3, . ~ . - Altitude:Condition), m_no3),
  lrt_row("sT:Altitude:Condition", m_no3lin, m_full),
  lrt_row("I(sT^2):Altitude:Condition", m_no3quad, m_full)
)

# 3. Backward elimination by likelihood-ratio tests -----------------------------------
# At each step, every term that is not contained in a higher-order term is tested
# by an LRT (model with vs without the term); the term with the largest P is
# removed if P > 0.05. Main effects are kept while they are part of an interaction.

# Step 1: full model
step1 <- dplyr::bind_rows(
  lrt_row("sT:Altitude:Condition", m_no3lin, m_full),
  lrt_row("I(sT^2):Altitude:Condition", m_no3quad, m_full)
) |> dplyr::mutate(Step = 1, .before = 1)
stopifnot(step1$Term[which.max(step1$P)] == "I(sT^2):Altitude:Condition", max(step1$P) > 0.05)
m_s1 <- m_no3quad # I(sT^2):Altitude:Condition removed

# Step 2
step2 <- dplyr::bind_rows(
  lrt_row("sT:Altitude:Condition", update(m_s1, . ~ . - sT:Altitude:Condition), m_s1),
  lrt_row("I(sT^2):Altitude", update(m_s1, . ~ . - I(sT^2):Altitude), m_s1),
  lrt_row("I(sT^2):Condition", update(m_s1, . ~ . - I(sT^2):Condition), m_s1)
) |> dplyr::mutate(Step = 2, .before = 1)
stopifnot(step2$Term[which.max(step2$P)] == "I(sT^2):Condition", max(step2$P) > 0.05)
m_s2 <- update(m_s1, . ~ . - I(sT^2):Condition)

# Step 3
step3 <- dplyr::bind_rows(
  lrt_row("sT:Altitude:Condition", update(m_s2, . ~ . - sT:Altitude:Condition), m_s2),
  lrt_row("I(sT^2):Altitude", update(m_s2, . ~ . - I(sT^2):Altitude), m_s2)
) |> dplyr::mutate(Step = 3, .before = 1)
stopifnot(step3$Term[which.max(step3$P)] == "sT:Altitude:Condition", max(step3$P) > 0.05)
m_s3 <- update(m_s2, . ~ . - sT:Altitude:Condition)

# Step 4
step4 <- dplyr::bind_rows(
  lrt_row("sT:Altitude", update(m_s3, . ~ . - sT:Altitude), m_s3),
  lrt_row("I(sT^2):Altitude", update(m_s3, . ~ . - I(sT^2):Altitude), m_s3),
  lrt_row("sT:Condition", update(m_s3, . ~ . - sT:Condition), m_s3),
  lrt_row("Altitude:Condition", update(m_s3, . ~ . - Altitude:Condition), m_s3)
) |> dplyr::mutate(Step = 4, .before = 1)
stopifnot(step4$Term[which.max(step4$P)] == "sT:Condition", max(step4$P) > 0.05)
m_s4 <- update(m_s3, . ~ . - sT:Condition)

# Step 5: all remaining interactions have P < 0.05, so elimination stops
step5 <- dplyr::bind_rows(
  lrt_row("sT:Altitude", update(m_s4, . ~ . - sT:Altitude), m_s4),
  lrt_row("I(sT^2):Altitude", update(m_s4, . ~ . - I(sT^2):Altitude), m_s4),
  lrt_row("Altitude:Condition", update(m_s4, . ~ . - Altitude:Condition), m_s4)
) |> dplyr::mutate(Step = 5, .before = 1)
stopifnot(max(step5$P) < 0.05)

tab_steps <- dplyr::bind_rows(step1, step2, step3, step4, step5) |>
  dplyr::group_by(Step) |>
  dplyr::mutate(Removed = P == max(P) & P > 0.05) |>
  dplyr::ungroup()
readr::write_csv(tab_steps, file.path(dir_tables, "01_Mon_backward_elimination.csv"))
print(tab_steps, n = Inf)

# Final model: (Time + Time^2) x Elevation + Elevation x Shading
m_final <- glmmTMB::glmmTMB(Monalonion ~ (sT + I(sT^2)) * Altitude + Altitude * Condition +
                              (1 | Farm) + (1 | Plot) + (1 | Event),
                            family = glmmTMB::nbinom1(), data = dM)
stopifnot(abs(logLik(m_final) - logLik(m_s4)) < 1e-6)

# 4. Parametric bootstrap of the critical tests ------------------------------------------

# All two-way interactions (linear and quadratic time), no three-way
m_2way <- glmmTMB::glmmTMB(Monalonion ~ (sT + I(sT^2)) * (Altitude + Condition) + Altitude:Condition +
                    (1 | Farm) + (1 | Plot) + (1 | Event),
                  family = glmmTMB::nbinom1(), data = dM)

block_pairs <- list(
  `Three-way: Time(linear+quadratic) x Elevation x Shading (vs all two-way)` = list(m_2way, m_full),
  `Three-way: Time(linear) x Elevation x Shading (step 3)` = list(m_s3, m_s2),
  `Quadratic time x Elevation (final model)` = list(update(m_final, . ~ . - I(sT^2):Altitude), m_final)
)
tab_blocks <- purrr::imap_dfr(block_pairs, function(p, nm) {
  if (run_bootstrap) {
    set.seed(2026) # simulations of the parametric bootstrap
    pb_lrt(p[[1]], p[[2]], dM, "Monalonion") |> dplyr::mutate(Test = nm, .before = 1)
  } else {
    lrt_row(nm, p[[1]], p[[2]]) |> dplyr::rename(Test = Term, P_asymptotic = P)
  }
}) |>
  dplyr::mutate(AIC_final = AIC(m_final), AIC_all_two_way = AIC(m_2way), AIC_full = AIC(m_full))
readr::write_csv(tab_blocks, file.path(dir_tables, "01_Mon_interaction_block_tests.csv"))
print(tab_blocks, width = Inf)

# =====================================================================================
# FINAL MODEL: (Time + Time^2) x Elevation + Elevation x Shading
# =====================================================================================

saveRDS(m_final, file.path(dir_models, "Mon_model_final.rds"))

# Type II LRT
m_noTA <- update(m_final, . ~ . - sT:Altitude)
m_noQA <- update(m_final, . ~ . - I(sT^2):Altitude)
m_noAC <- update(m_final, . ~ . - Altitude:Condition)
m_noAlt <- update(m_final, . ~ . - sT:Altitude - I(sT^2):Altitude - Altitude:Condition)
lrt_final <- dplyr::bind_rows(
  lrt_row("sT", update(m_noTA, . ~ . - sT), m_noTA),
  lrt_row("I(sT^2)", update(m_noQA, . ~ . - I(sT^2)), m_noQA),
  lrt_row("Altitude", update(m_noAlt, . ~ . - Altitude), m_noAlt),
  lrt_row("Condition", update(m_noAC, . ~ . - Condition), m_noAC),
  lrt_row("sT:Altitude", m_noTA, m_final),
  lrt_row("I(sT^2):Altitude", m_noQA, m_final),
  lrt_row("Altitude:Condition", m_noAC, m_final)
)
readr::write_csv(lrt_final, file.path(dir_tables, "01_Mon_final_LRT.csv"))
print(lrt_final)

# Variance components of the random effects
vc_final <- glmmTMB::VarCorr(m_final)$cond
tibble::tibble(Group = names(vc_final), Variance = purrr::map_dbl(vc_final, ~ .x[1, 1])) |>
  dplyr::mutate(SD = sqrt(Variance)) |>
  readr::write_csv(file.path(dir_tables, "01_Mon_final_variance_components.csv"))

# Residual diagnostics (DHARMa)
set.seed(2026)
res_final <- DHARMa::simulateResiduals(m_final, n = 1000)
save_dharma_plot(res_final, "S4_Mon_DHARMa", "Monalonion - final model (NB1 GLMM)")
tests_final <- list(Dispersion = DHARMa::testDispersion(res_final, plot = FALSE),
               `Uniformity (KS)` = DHARMa::testUniformity(res_final, plot = FALSE),
               `Zero inflation` = DHARMa::testZeroInflation(res_final, plot = FALSE),
               Outliers = DHARMa::testOutliers(res_final, type = "bootstrap", plot = FALSE))
tibble::tibble(Test = names(tests_final),
       Statistic = purrr::map_dbl(tests_final, ~ .x$statistic[[1]]),
       P = purrr::map_dbl(tests_final, "p.value")) |>
  readr::write_csv(file.path(dir_tables, "01_Mon_final_DHARMa_tests.csv"))

# Elevation within shading condition (at the mean sampling day, sT = 0) ------------
em_ec <- emmeans::emmeans(m_final, ~ Altitude | Condition, type = "response")
cld_ec <- cld_tbl(em_ec)
# Shading within elevation
ct_ce <- emmeans::contrast(emmeans::emmeans(m_final, ~ Condition | Altitude), "pairwise") |>
  summary(infer = TRUE, type = "response") |>
  tibble::as_tibble()
readr::write_csv(cld_ec, file.path(dir_tables, "01_Mon_final_Tukey_elevation_within_condition.csv"))
readr::write_csv(ct_ce, file.path(dir_tables, "01_Mon_final_shaded_vs_unshaded_within_elevation.csv"))

# Elevation within sampling day (averaged over shading) ----------------------------
days_obs <- sort(unique(dM$Time))
sT_obs <- (days_obs - time_center) / time_scale
em_et <- emmeans::emmeans(m_final, ~ Altitude | sT, at = list(sT = sT_obs), type = "response")
cld_et <- cld_tbl(em_et) |> dplyr::mutate(Time = round(sT * time_scale + time_center))
readr::write_csv(cld_et, file.path(dir_tables, "01_Mon_final_Tukey_elevation_within_day.csv"))

# Curves and peaks ---------------------------------------------------------------------
# Estimated mean number of insects per tree for every day (fixed effects)
grid_days <- seq(first_day, max(dM$Time), by = 1)
sT_grid <- (grid_days - time_center) / time_scale

curv_ec <- emmeans::emmeans(m_final, ~ Altitude * Condition * sT, at = list(sT = sT_grid),
                    type = "response") |>
  tibble::as_tibble() |>
  dplyr::mutate(Time = sT * time_scale + time_center) |>
  dplyr::rename(Predicted = response, Lower = asymp.LCL, Upper = asymp.UCL)

curv_e <- emmeans::emmeans(m_final, ~ Altitude * sT, at = list(sT = sT_grid), type = "response") |>
  tibble::as_tibble() |>
  dplyr::mutate(Time = sT * time_scale + time_center) |>
  dplyr::rename(Predicted = response, Lower = asymp.LCL, Upper = asymp.UCL)

peaks <- dplyr::bind_rows(
  curv_e |> dplyr::group_by(Altitude) |> dplyr::slice_max(Predicted, n = 1) |> dplyr::mutate(Condition = "Average"),
  curv_ec |> dplyr::group_by(Altitude, Condition) |> dplyr::slice_max(Predicted, n = 1) |>
    dplyr::mutate(Condition = as.character(Condition))) |>
  dplyr::mutate(Date = start + Time) |>
  dplyr::select(Altitude, Condition, Time, Date, Predicted, Lower, Upper)
readr::write_csv(peaks, file.path(dir_tables, "01_Mon_final_peaks.csv"))
readr::write_csv(curv_ec, file.path(dir_tables, "01_Mon_final_curves_elevation_condition.csv"))
readr::write_csv(curv_e, file.path(dir_tables, "01_Mon_final_curves_elevation.csv"))
print(peaks)

# Figures --------------------------------------------------------------------------------

# Upper limit of the y axes: the multiple of 15 above the largest confidence limit
# (ggplot drops values outside the limits, which would cut ribbons and error bars)
ymax2 <- 15 * ceiling(max(curv_ec$Upper, cld_ec$asymp.UCL) / 15) # Figure 2
ymax3 <- 15 * ceiling(max(75, curv_e$Upper, cld_et$asymp.UCL) / 15) # Figure 3

curve_panel <- function(cv, facet = TRUE, show_month = TRUE, ymax) {
  p <- ggplot2::ggplot(cv, ggplot2::aes(Time, Predicted, ymin = Lower, ymax = Upper,
                      fill = Altitude, linetype = Altitude)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = first_day, ymin = -Inf, ymax = Inf, alpha = .15) +
    ggplot2::geom_vline(xintercept = first_day, linetype = "dashed", colour = "darkred", linewidth = 0.4) +
    ggplot2::geom_ribbon(alpha = 0.5) + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_fill_manual(values = elev_cols) + ggplot2::scale_linetype_manual(values = elev_lty) +
    ggplot2::scale_x_continuous(limits = c(0, 170), breaks = seq(0, 175, 25),
                       sec.axis = if (show_month) month_axis(start) else ggplot2::waiver()) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, 15)) +
    ggplot2::labs(x = "Time (days)", y = "Number of insects", fill = "Elevation", linetype = "Elevation") +
    theme_bugs()
  if (facet) p <- p + ggplot2::facet_wrap(~Condition)
  p
}

# Figure 2 --------------------------------------------------------------------
fig2a <- curve_panel(curv_ec, ymax = ymax2)

ast <- ct_ce |> dplyr::transmute(Altitude, star = ifelse(p.value < 0.05, "*", ""))
bars_ec <- cld_ec |>
  dplyr::left_join(ast, by = "Altitude") |>
  dplyr::mutate(Altitude = factor(Altitude, levels = rev(elev_levels)))

fig2b_side <- function(cond, reverse) {
  dd <- dplyr::filter(bars_ec, Condition == cond)
  p <- ggplot2::ggplot(dd, ggplot2::aes(Altitude, response, fill = Altitude)) +
    ggplot2::geom_col(alpha = 0.75, colour = "black", linewidth = 0.2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = .group, y = asymp.UCL), hjust = if (reverse) 1.6 else -0.6, size = 4) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = elev_cols, breaks = elev_levels) +
    ggplot2::labs(y = "Number of insects", x = NULL, fill = "Elevation", subtitle = cond) +
    theme_bugs() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
          plot.subtitle = ggplot2::element_text(hjust = 0.5))
  if (reverse) {
    p + ggplot2::scale_y_reverse(limits = c(ymax2, 0), breaks = seq(0, ymax2, 15)) +
      # asterisk (shaded vs unshaded) to the left of the Tukey letter, clear of the error bar
      ggplot2::geom_text(ggplot2::aes(label = star, y = asymp.UCL + 0.07 * ymax2), hjust = 1, size = 6, vjust = 0.75)
  } else {
    p + ggplot2::scale_y_continuous(limits = c(0, ymax2), breaks = seq(0, ymax2, 15))
  }
}
centre <- ggplot2::ggplot(tibble::tibble(x = 1, y = factor(rev(elev_levels), levels = rev(elev_levels))),
                 ggplot2::aes(x, y, label = y)) +
  ggplot2::geom_text(size = 4) +
  ggplot2::theme_void()
fig2b <- (fig2b_side("Shaded", TRUE) | centre | fig2b_side("Unshaded", FALSE)) +
  patchwork::plot_layout(widths = c(1, 0.22, 1), guides = "collect") & ggplot2::theme(legend.position = "top")

fig2 <- patchwork::wrap_elements(fig2a) / patchwork::wrap_elements(fig2b) + patchwork::plot_layout(heights = c(1.1, 1)) +
  patchwork::plot_annotation(tag_levels = list(c("(a)", "(b)")))
save_fig(fig2, "Figure_2_Monalonion", 9, 10)

# Figure 3 --------------------------------------------------------------------
fig3a <- curve_panel(curv_e, facet = FALSE, ymax = ymax3)
fig3b <- cld_et |>
  dplyr::mutate(Altitude = factor(Altitude, levels = elev_levels)) |>
  ggplot2::ggplot(ggplot2::aes(factor(Time), response, fill = Altitude)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(0.9), alpha = 0.75, colour = "black", linewidth = 0.2) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL), position = ggplot2::position_dodge(0.9), width = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = .group, y = asymp.UCL), position = ggplot2::position_dodge(0.9), vjust = -0.5, size = 3) +
  ggplot2::scale_fill_manual(values = elev_cols) +
  ggplot2::scale_y_continuous(limits = c(0, ymax3), breaks = seq(0, ymax3, 15)) +
  ggplot2::labs(x = "Time (days)", y = "Number of insects", fill = "Elevation") +
  theme_bugs()
fig3 <- (fig3a | fig3b) + patchwork::plot_annotation(tag_levels = list(c("(a)", "(b)")))
save_fig(fig3, "Figure_3_Monalonion", 12, 5.5)

sessioninfo <- capture.output(sessionInfo())
writeLines(sessioninfo, file.path(dir_tables, "sessionInfo.txt"))
