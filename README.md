# E-Governance Data Analysis in R: UN E-Government Development Index (EGDI) 2024

Four-week internship project: acquisition, transformation, visualization and predictive modeling of public e-governance data.

**Main dataset:** UN DESA E-Government Development Index 2024 (193 UN Member States)
Source: https://publicadministration.un.org/egovkb/Data-Center
**Week 4 additions:** World Bank GDP and population (via github.com/datasets/gdp and github.com/datasets/population)

## Structure
| Path | Description |
|---|---|
| `egov_analysis.R` | Week 1: import, cleaning, EDA, charts, Wilcoxon test |
| `week2_transform.R` | Week 2: dplyr/tidyr pipeline (joins, derived variables, aggregation, reshaping) |
| `week3_visualize.R` | Week 3: seven ggplot2 charts with design iteration |
| `week4_modeling.R` | Week 4: confidence intervals, ANOVA, t-test, correlation, linear and logistic models, cross-validation, diagnostics |
| `data/` | Raw, lookup, cleaned, transformed and model-output CSVs |
| `figures/` | All charts (fig* Week 1, w2_* pipeline, w3_* Week 3, w4_* Week 4) |
| `reports/` | Word reports for Weeks 1 to 4 |

## How to run (from the project folder, in this order)
```r
install.packages(c("readr", "dplyr", "tidyr", "ggplot2", "scales", "stringr"))
source("egov_analysis.R")
source("week2_transform.R")
source("week3_visualize.R")
source("week4_modeling.R")
```

## Notes
- `data/country_regions.csv`, `data/south_asia_peers.csv` and `data/wb_name_map.csv` are lookup tables compiled for this project.
- Week 4 uses the latest World Bank GDP value from 2019-2023 for each country; five countries without recent data are excluded from the models.
- `set.seed(2024)` makes cross-validation and the train/test split reproducible.
