# Bugs_Shade_Altitude

R code for the analysis of the effects of elevation, shading condition and time on
the abundance of the cacao pest bugs *Monalonion* (Miridae) and *Antiteuchus*
(Pentatomidae) in Amazonas, Peru.

## Data

The field data are not included in this repository; they are available from the
corresponding author on reasonable request. To run the analysis, place the
spreadsheet in

```
Data/raw_data/Dados_Percevejos_Antiteuchus_Monalonion.xlsx
```

The expected sheet and columns are described in
[`Data/raw_data/README.md`](Data/raw_data/README.md).

## Project structure

```
Data/
  raw_data/        field spreadsheet (not included, see above)
  tidy_data/       bugs_tidy.csv (one row per tree), created by 00_Data_Preparation.R
R_Scripts/
  functions.R                    packages, folders, colors and helpers that format/save results
  00_Data_Preparation.R          tidy data, design summary, descriptive tables
  01_Monalonion.R                Monalonion: NB GLMMs, backward elimination and final model
  02_Antiteuchus.R               Antiteuchus: NB GLMMs, backward elimination and final model
  03_Sensitivity_Monalonion.R    Monalonion with all 12 sampling occasions
  run_all.R                      runs everything in order
Outputs/
  tables/          csv results (tests, Tukey comparisons, peaks, curves, AIC, diagnostics)
  figures/         svg/ pdf/ jpg/  figures
  diagnostics/     svg/ pdf/ jpg/  DHARMa residual plots
  models/          fitted models (rds; created when the scripts are run)
```

## Statistical analysis

Tree-level counts were modeled with negative binomial (NB1) GLMMs (`glmmTMB`) with
random intercepts for farm, shading condition within farm and sampling date within
farm. The main effects of time (linear and quadratic), elevation and
shading were kept in all models, and non-significant interactions were removed by backward
elimination with likelihood-ratio tests (the same procedure for both species); the
critical tests for *Monalonion* were also checked by parametric bootstrap. Residuals were checked with `DHARMa`, and Tukey
comparisons were made with `emmeans`.

For *Monalonion* (days 46–167) the final model is (Time + Time²) × Elevation +
Elevation × Shading. The three-way interaction (Time + Time²) × Elevation × Shading
received only weak support (parametric bootstrap P = 0.050; P = 0.090 with all
sampling occasions) and was not retained. For *Antiteuchus* all interactions were
removed, leaving the main effects of time, time², elevation and shading.

## How to run

1. Place the data file as described above.
2. Open `Manuel.Rproj` in RStudio (or start R in the project folder).
3. Run `source("R_Scripts/run_all.R")`.

Required packages: tidyverse, here, readxl, glmmTMB, emmeans, multcomp,
multcompView, DHARMa, patchwork and svglite. The parametric bootstrap runs in
parallel on a PSOCK cluster (`parallel::makeCluster`), which works on any operating
system and IDE; set `run_bootstrap <- FALSE` in `01_Monalonion.R` and
`03_Sensitivity_Monalonion.R` to skip it.
