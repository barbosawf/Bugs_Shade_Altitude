# 00 - Data preparation ------------------------------------------------------
# Reads the raw field spreadsheet and writes the tidy data set used by all
# analyses (Data/tidy_data/bugs_tidy.csv), plus design/descriptive tables.

source(here::here("R_Scripts", "functions.R"))

# Raw data -----------------------------------------------------------------------

raw <- readxl::read_excel(here::here("Data", "raw_data", "Dados_Percevejos_Antiteuchus_Monalonion.xlsx"),
                  sheet = "Planilha1")

bugs <- raw |>
  dplyr::transmute(
    Date = as.Date(Date),
    Month,
    Time = as.integer(Time), # days since the first sampling
    Altitude = factor(Altitude, levels = elev_levels),
    Farm,
    Condition = factor(Condition, levels = c("Shaded", "Unshaded")),
    Tree = ID, # each tree sampled only once
    Plot = paste(Farm, Condition, sep = "_"),
    Event = paste(Farm, sprintf("%03d", Time), sep = "_"),
    Monalonion = as.integer(Monalonion),
    Antiteuchus = as.integer(Antiteuchus),
    Temperature, Humidity
  )

stopifnot(!any(duplicated(bugs$Tree))) # each tree sampled only once

readr::write_csv(bugs, file.path(dir_tidy, "bugs_tidy.csv"))

# Design summary ---------------------------------------------------------------

design <- tibble::tibble(
  Item = c("Farms", "Farms per elevation", "Sampling occasions",
           "Days between occasions", "First / last sampling day",
           "Trees per farm x condition x occasion",
           "Trees per occasion", "Total trees sampled",
           "Repeated unit (farm x condition)"),
  Value = c(dplyr::n_distinct(bugs$Farm),
            paste(unique(dplyr::count(dplyr::distinct(bugs, Farm, Altitude), Altitude)$n), collapse = "/"),
            dplyr::n_distinct(bugs$Time),
            paste(range(diff(sort(unique(bugs$Time)))), collapse = "-"),
            paste(range(bugs$Time), collapse = " / "),
            paste(unique(dplyr::count(bugs, Farm, Condition, Time)$n), collapse = "/"),
            paste(unique(dplyr::count(bugs, Time)$n), collapse = "/"),
            nrow(bugs),
            dplyr::n_distinct(bugs$Plot))
)
readr::write_csv(design, file.path(dir_tables, "00_design_summary.csv"))
print(design)

# Sampling calendar
calendar <- bugs |> dplyr::distinct(Time, Date, Month) |> dplyr::arrange(Time)
readr::write_csv(calendar, file.path(dir_tables, "00_sampling_calendar.csv"))

# Descriptive statistics per occasion, elevation and condition (per tree)
descr <- bugs |>
  tidyr::pivot_longer(c(Monalonion, Antiteuchus), names_to = "Genus", values_to = "Count") |>
  dplyr::group_by(Genus, Time, Altitude, Condition) |>
  dplyr::summarise(n_trees = dplyr::n(), Total = sum(Count), Mean = mean(Count),
            SD = sd(Count), Zeros = sum(Count == 0), .groups = "drop")
readr::write_csv(descr, file.path(dir_tables, "00_descriptive_counts.csv"))

# Zero counts per occasion
zeros <- bugs |>
  dplyr::group_by(Time, Date) |>
  dplyr::summarise(n_trees = dplyr::n(),
            Monalonion_total = sum(Monalonion),
            Monalonion_zero_trees = sum(Monalonion == 0),
            Antiteuchus_total = sum(Antiteuchus),
            Antiteuchus_zero_trees = sum(Antiteuchus == 0), .groups = "drop")
readr::write_csv(zeros, file.path(dir_tables, "00_zero_counts_by_occasion.csv"))
print(zeros, n = Inf)

# Descriptive plots ---------------------------------------------------------------

p_descr <- bugs |>
  tidyr::pivot_longer(c(Monalonion, Antiteuchus), names_to = "Genus", values_to = "Count") |>
  ggplot2::ggplot(ggplot2::aes(factor(Time), Count, colour = Altitude)) +
  ggplot2::geom_boxplot(outlier.size = 0.8) +
  ggplot2::facet_grid(Genus ~ Condition, scales = "free_y") +
  ggplot2::scale_colour_manual(values = elev_cols) +
  ggplot2::labs(x = "Time (days)", y = "Number of insects per tree", colour = "Elevation") +
  theme_bugs()
save_fig(p_descr, "00_descriptive_boxplots", 11, 7)
