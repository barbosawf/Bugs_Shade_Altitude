# 02 - Antiteuchus ---------------------------------------------------------------
#
# Response: number of Antiteuchus individuals per tree (all 12 occasions).
# Model: negative binomial GLMM (glmmTMB) with random intercepts for
#   Farm  - farms (two per elevation)
#   Plot  - farm x shading condition (group of trees observed on every occasion)
#   Event - farm x sampling date
# Fixed effects are tested by Type II likelihood-ratio tests (each term is tested
# against the model with all terms that do not contain it). Non-significant
# interactions are removed (final model: main effects of time, time^2, elevation
# and shading).

source(here::here("R_Scripts", "functions.R"))

# Data -------------------------------------------------------------------------------

bugs <- read_bugs()
start <- min(bugs$Date)
dA <- bugs

# Time is standardized (sT) to improve the numerical stability of the GLMMs
time_center <- mean(dA$Time)
time_scale <- sd(dA$Time)
dA <- dA |> dplyr::mutate(sT = (Time - time_center) / time_scale)

# 1. Error distribution and random-effect structure ----------------------------------
# All candidates have the full fixed part (Time + Time^2) x Elevation x Shading.

m_pois_olre <- glmmTMB::glmmTMB(
  Antiteuchus ~ (sT + I(sT^2)) * Altitude * Condition +
    (1 | Farm) + (1 | Plot) + (1 | Event) + (1 | Tree),
  family = poisson(), data = dA)

m_nb1 <- glmmTMB::glmmTMB(
  Antiteuchus ~ (sT + I(sT^2)) * Altitude * Condition +
    (1 | Farm) + (1 | Plot) + (1 | Event),
  family = glmmTMB::nbinom1(), data = dA)

m_nb2 <- glmmTMB::glmmTMB(
  Antiteuchus ~ (sT + I(sT^2)) * Altitude * Condition +
    (1 | Farm) + (1 | Plot) + (1 | Event),
  family = glmmTMB::nbinom2(), data = dA)

m_nb1_glm <- glmmTMB::glmmTMB(
  Antiteuchus ~ (sT + I(sT^2)) * Altitude * Condition,
  family = glmmTMB::nbinom1(), data = dA)

m_nb1_ar1 <- glmmTMB::glmmTMB(
  Antiteuchus ~ (sT + I(sT^2)) * Altitude * Condition +
    (1 | Farm) + ar1(fTime + 0 | Plot),
  family = glmmTMB::nbinom1(), data = dA)

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
readr::write_csv(tab_cands, file.path(dir_tables, "02_Ant_candidate_models_AIC.csv"))
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
m_alt_0 <- glmmTMB::glmmTMB(Antiteuchus ~ (sT + I(sT^2)) * Condition +
                     (1 | Farm) + (1 | Plot) + (1 | Event),
                   family = glmmTMB::nbinom1(), data = dA)
m_alt_1 <- update(m_alt_0, . ~ . + Altitude)
m_cond_0 <- glmmTMB::glmmTMB(Antiteuchus ~ (sT + I(sT^2)) * Altitude +
                      (1 | Farm) + (1 | Plot) + (1 | Event),
                    family = glmmTMB::nbinom1(), data = dA)
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
readr::write_csv(lrt_full, file.path(dir_tables, "02_Ant_LRT_full_model.csv"))
print(lrt_full)

# All interactions jointly: main-effects model vs full model
m_main <- glmmTMB::glmmTMB(Antiteuchus ~ sT + I(sT^2) + Altitude + Condition +
                    (1 | Farm) + (1 | Plot) + (1 | Event),
                  family = glmmTMB::nbinom1(), data = dA)
tab_int <- lrt_row("All interactions jointly (vs main-effects model)", m_main, m_full) |>
  dplyr::rename(Test = Term) |>
  dplyr::mutate(AIC_main = AIC(m_main), AIC_full = AIC(m_full))
readr::write_csv(tab_int, file.path(dir_tables, "02_Ant_all_interactions_test.csv"))
print(tab_int)

# 3. Final model: sT + sT^2 + Altitude + Condition -----------------------------------

saveRDS(m_main, file.path(dir_models, "Ant_model_final.rds"))

# Type II LRT (no interactions, so each term is dropped from the final model)
lrt_main <- dplyr::bind_rows(
  lrt_row("sT", update(m_main, . ~ . - sT), m_main),
  lrt_row("I(sT^2)", update(m_main, . ~ . - I(sT^2)), m_main),
  lrt_row("Altitude", update(m_main, . ~ . - Altitude), m_main),
  lrt_row("Condition", update(m_main, . ~ . - Condition), m_main)
)
readr::write_csv(lrt_main, file.path(dir_tables, "02_Ant_final_LRT.csv"))
print(lrt_main)

# Variance components of the random effects
vc <- glmmTMB::VarCorr(m_main)$cond
tibble::tibble(Group = names(vc), Variance = purrr::map_dbl(vc, ~ .x[1, 1])) |>
  dplyr::mutate(SD = sqrt(Variance)) |>
  readr::write_csv(file.path(dir_tables, "02_Ant_variance_components.csv"))

# Residual diagnostics (DHARMa)
set.seed(2026)
res <- DHARMa::simulateResiduals(m_main, n = 1000)
save_dharma_plot(res, "S5_Ant_DHARMa", "Antiteuchus - final model (NB1 GLMM)")
tests <- list(Dispersion = DHARMa::testDispersion(res, plot = FALSE),
              `Uniformity (KS)` = DHARMa::testUniformity(res, plot = FALSE),
              `Zero inflation` = DHARMa::testZeroInflation(res, plot = FALSE),
              Outliers = DHARMa::testOutliers(res, type = "bootstrap", plot = FALSE))
tibble::tibble(Test = names(tests),
       Statistic = purrr::map_dbl(tests, ~ .x$statistic[[1]]),
       P = purrr::map_dbl(tests, "p.value")) |>
  readr::write_csv(file.path(dir_tables, "02_Ant_DHARMa_tests.csv"))

# Tukey comparisons (at the mean sampling day, sT = 0)
em_c <- emmeans::emmeans(m_main, ~Condition, type = "response")
cld_c <- cld_tbl(em_c)
em_e <- emmeans::emmeans(m_main, ~Altitude, type = "response")
cld_e <- cld_tbl(em_e)
readr::write_csv(cld_c, file.path(dir_tables, "02_Ant_Tukey_condition.csv"))
readr::write_csv(cld_e, file.path(dir_tables, "02_Ant_Tukey_elevation.csv"))
readr::write_csv(tibble::as_tibble(summary(pairs(em_e), infer = TRUE)),
          file.path(dir_tables, "02_Ant_pairs_elevation.csv"))

# Curves and peaks: estimated mean number of insects per tree for every day
grid_days <- seq(0, max(dA$Time), by = 1)
sT_grid <- (grid_days - time_center) / time_scale

curv_c <- emmeans::emmeans(m_main, ~ Condition * sT, at = list(sT = sT_grid), type = "response") |>
  tibble::as_tibble() |>
  dplyr::mutate(Time = sT * time_scale + time_center) |>
  dplyr::rename(Predicted = response, Lower = asymp.LCL, Upper = asymp.UCL)

curv_e <- emmeans::emmeans(m_main, ~ Altitude * sT, at = list(sT = sT_grid), type = "response") |>
  tibble::as_tibble() |>
  dplyr::mutate(Time = sT * time_scale + time_center) |>
  dplyr::rename(Predicted = response, Lower = asymp.LCL, Upper = asymp.UCL)

peaks <- dplyr::bind_rows(curv_c |> dplyr::group_by(Condition) |> dplyr::slice_max(Predicted, n = 1) |>
                     dplyr::mutate(Group = as.character(Condition)),
                   curv_e |> dplyr::group_by(Altitude) |> dplyr::slice_max(Predicted, n = 1) |>
                     dplyr::mutate(Group = as.character(Altitude))) |>
  dplyr::ungroup() |>
  dplyr::mutate(Date = start + Time) |>
  dplyr::select(Group, Time, Date, Predicted, Lower, Upper)
readr::write_csv(peaks, file.path(dir_tables, "02_Ant_peaks.csv"))
readr::write_csv(curv_c, file.path(dir_tables, "02_Ant_curves_condition.csv"))
readr::write_csv(curv_e, file.path(dir_tables, "02_Ant_curves_elevation.csv"))
print(peaks)

# Figure 4 -------------------------------------------------------------------------------

ymax <- 75
curve_plot <- function(cv, grp, cols, ltys, lab) {
  ggplot2::ggplot(cv, ggplot2::aes(Time, Predicted, ymin = Lower, ymax = Upper,
                 fill = .data[[grp]], linetype = .data[[grp]])) +
    ggplot2::geom_ribbon(alpha = 0.5) + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_linetype_manual(values = ltys) +
    ggplot2::scale_x_continuous(limits = c(0, 170), breaks = seq(0, 175, 25),
                       sec.axis = month_axis(start)) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, 15)) +
    ggplot2::labs(x = "Time (days)", y = "Number of insects", fill = lab, linetype = lab) +
    theme_bugs()
}
bar_plot <- function(cl, grp, cols, lab) {
  cl |> dplyr::mutate(dplyr::across(dplyr::all_of(grp), ~ factor(.x, levels = names(cols)))) |>
    ggplot2::ggplot(ggplot2::aes(.data[[grp]], response, fill = .data[[grp]])) +
    ggplot2::geom_col(alpha = 0.75, colour = "black", linewidth = 0.2, show.legend = FALSE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = .group, y = asymp.UCL), vjust = -0.8, size = 4) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, 15)) +
    ggplot2::labs(x = lab, y = "Number of insects") +
    theme_bugs()
}

fig4 <- (curve_plot(curv_c, "Condition", cond_cols, cond_lty, "Condition") |
           bar_plot(cld_c, "Condition", cond_cols, "Condition")) /
  (curve_plot(curv_e, "Altitude", elev_cols, elev_lty, "Elevation") |
     bar_plot(cld_e, "Altitude", elev_cols, "Elevation")) +
  patchwork::plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)")))
save_fig(fig4, "Figure_4_Antiteuchus", 10, 9)
