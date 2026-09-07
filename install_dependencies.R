required_packages <- c(
  "conflicted",
  "countrycode",
  "data.table",
  "dplyr",
  "future",
  "future.apply",
  "GWmodel",
  "here",
  "mice",
  "purrr",
  "readr",
  "readxl",
  "rnaturalearth",
  "rnaturalearthdata",
  "scales",
  "sf",
  "sn",
  "sp",
  "spdep",
  "stringr",
  "tibble",
  "tidyr",
  "WDI"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  install.packages(missing_packages, dependencies = TRUE)
} else {
  message("All required R packages are already installed.")
}
