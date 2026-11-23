# Shared settings and small helpers -----------------------------------------------
# Sourced by the analysis scripts. Model fitting, tests and estimates are written
# out in the scripts themselves; this file only holds packages, folders, graphic
# settings and helpers that format or save results.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(glmmTMB)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(DHARMa)
  library(patchwork)
  library(svglite)
})

emmeans::emm_options(lmer.df = "asymptotic")

# Folders ------------------------------------------------------------------------

dir_tidy <- here::here("Data", "tidy_data")
dir_tables <- here::here("Outputs", "tables")
dir_models <- here::here("Outputs", "models")
dir_fig <- here::here("Outputs", "figures")
dir_diag <- here::here("Outputs", "diagnostics")

for (d in c(dir_tidy, dir_tables, dir_models,
            file.path(dir_fig, c("svg", "pdf", "jpg")),
            file.path(dir_diag, c("svg", "pdf", "jpg")))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Graphic settings -------------------------------------------------------------------

elev_levels <- c("High", "Medium", "Low")
elev_cols <- c(High = "#00AFBB", Medium = "#E7B800", Low = "#999999")
elev_lty <- c(High = "solid", Medium = "dotted", Low = "dashed")
cond_cols <- c(Shaded = "#00AFBB", Unshaded = "#E7B800")
cond_lty <- c(Shaded = "solid", Unshaded = "dotted")

# Month labels at the 15th of each month (day 0 = first sampling date)
month_axis <- function(start_date) {
  m15 <- seq(as.Date(start_date), by = "month", length.out = 6)
  ggplot2::sec_axis(~., breaks = as.numeric(m15 - as.Date(start_date)),
           labels = format(m15, "%b"), name = NULL)
}

theme_bugs <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top",
          axis.title = ggplot2::element_text(size = 13),
          axis.text = ggplot2::element_text(size = 11, colour = "black"),
          axis.ticks.length = ggplot2::unit(.15, "cm"),
          legend.title = ggplot2::element_text(size = 11),
          legend.text = ggplot2::element_text(size = 10),
          strip.text = ggplot2::element_text(size = 11),
          plot.tag = ggplot2::element_text(size = 15))
}

# Saves a ggplot as svg, pdf and jpg (600 dpi)
save_fig <- function(plot, name, width, height, folder = dir_fig) {
  ggplot2::ggsave(file.path(folder, "svg", paste0(name, ".svg")), plot,
         width = width, height = height, device = svglite::svglite)
  ggplot2::ggsave(file.path(folder, "pdf", paste0(name, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(file.path(folder, "jpg", paste0(name, ".jpg")), plot,
         width = width, height = height, dpi = 600)
  invisible(plot)
}

# Data -------------------------------------------------------------------------------

read_bugs <- function() {
  readr::read_csv(file.path(dir_tidy, "bugs_tidy.csv"), show_col_types = FALSE) |>
    dplyr::mutate(Altitude = factor(Altitude, levels = elev_levels),
           Condition = factor(Condition, levels = c("Shaded", "Unshaded")),
           Farm = factor(Farm),
           Plot = factor(Plot), # farm x shading condition (repeated unit)
           Event = factor(Event), # farm x sampling date
           fTime = factor(Time))
}

# Result formatting ----------------------------------------------------------------

# One row with the likelihood-ratio test of anova(m0, m1)
lrt_row <- function(term, m0, m1) {
  a <- anova(m0, m1)
  tibble::tibble(Term = term, LR_chisq = a$Chisq[2], df = a$`Chi Df`[2], P = a$`Pr(>Chisq)`[2])
}

# Compact letter display of an emmGrid as a tibble
# (Tukey-adjusted pairwise tests; unadjusted 95% confidence intervals)
cld_tbl <- function(em) {
  multcomp::cld(em, Letters = letters, reversed = TRUE) |>
    tibble::as_tibble() |>
    dplyr::mutate(.group = stringr::str_trim(.group))
}

# Parametric bootstrap of the likelihood-ratio test between two nested fitted
# models: responses are simulated from m0, both models are refitted to each
# simulated data set, and the observed LR is compared with the simulated LRs.
# Call set.seed() before it for reproducible results.
pb_lrt <- function(m0, m1, data, response, nsim = 1000,
                   cores = max(1, parallel::detectCores() - 2)) {
  obs <- as.numeric(2 * (logLik(m1) - logLik(m0)))
  sims <- simulate(m0, nsim = nsim)
  lr <- unlist(parallel::mclapply(seq_len(nsim), function(i) {
    dd <- data
    dd[[response]] <- sims[[i]]
    a <- try(suppressWarnings(update(m0, data = dd)), silent = TRUE)
    b <- try(suppressWarnings(update(m1, data = dd)), silent = TRUE)
    if (inherits(a, "try-error") || inherits(b, "try-error")) return(NA_real_)
    max(0, as.numeric(2 * (logLik(b) - logLik(a))))
  }, mc.cores = cores))
  tibble::tibble(LR_chisq = obs,
         df = attr(logLik(m1), "df") - attr(logLik(m0), "df"),
         P_asymptotic = pchisq(obs, df, lower.tail = FALSE),
         P_bootstrap = (sum(lr >= obs, na.rm = TRUE) + 1) / (sum(!is.na(lr)) + 1),
         n_sim_ok = sum(!is.na(lr)))
}

# Saves the DHARMa residual plots (svg, pdf, jpg)
save_dharma_plot <- function(res, name, title) {
  draw <- function() {
    plot(res)
    title(title, outer = TRUE, line = -1)
  }
  svglite::svglite(file.path(dir_diag, "svg", paste0(name, ".svg")), 10, 5); draw(); dev.off()
  cairo_pdf(file.path(dir_diag, "pdf", paste0(name, ".pdf")), 10, 5); draw(); dev.off()
  jpeg(file.path(dir_diag, "jpg", paste0(name, ".jpg")), 10, 5, units = "in", res = 300); draw(); dev.off()
}
