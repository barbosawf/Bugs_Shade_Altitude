# Runs the complete analysis in order ----------------------------------------------
scripts <- c("00_Data_Preparation.R", "01_Monalonion.R", "02_Antiteuchus.R",
             "03_Sensitivity_Monalonion.R")
for (s in scripts) {
  message("Running ", s)
  source(here::here("R_Scripts", s), local = new.env())
}
