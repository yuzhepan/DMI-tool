# DMI-tool

DMI-tool provides an interactive browser calculator and the reproducible R analysis code accompanying the manuscript **“A probabilistic framework reveals life course heterogeneity in global diet and healthy longevity associations.”**

## About the Dietary Matching Index

The Dietary Matching Index (DMI) is a probabilistic framework for evaluating concordance with established dietary patterns. For each dietary component, observed intake is mapped to a guideline-concordance probability using pattern-specific reference values. These component probabilities are then combined with a weighted generalized mean. The aggregation gives greater influence to poorly matched components and reduces the extent to which strong performance in one component can offset a marked shortfall in another.

The analysis implements DMI versions of eight dietary patterns: AHEI, MSDPS, hPDI, DASH, DRRD, WCRF, PHDI and DBI. It evaluates age-, sex-, time- and geography-specific variation and examines country-level associations with population-weighted healthy life expectancy (HALE).

## Repository contents

- `index.html`: browser-based DMI calculator with component probabilities and improvement summaries.
- `DMI_reproducible_analysis.R`: complete paper analysis, from country-level dietary inputs through DMI construction, uncertainty propagation, spatial analyses and geographically and temporally weighted regression (GTWR).
- `install_dependencies.R`: installs the R packages required by the analysis.
- `DATA.md`: input-data structure, filenames and source notes.
- `CITATION.cff`: citation metadata for this repository.
- `data/cn_foods/`: food-composition resources used by the browser calculator.
- `data/vendor/`: local JavaScript dependencies used by the browser calculator.

## Analysis workflow

The R script performs the following steps:

1. Reads and harmonizes country-level dietary estimates.
2. Constructs conventional scores and DMIs for eight dietary patterns.
3. Propagates dietary-input uncertainty with 500 Monte Carlo draws.
4. Reconstructs population-weighted HALE and life expectancy (LE) for the study age groups.
5. Handles missing covariate data using multiple imputation by chained equations with predictive mean matching (`m = 20`).
6. Calculates uncertainty-aware global Moran's I and local indicators of spatial association.
7. Fits corrected great-circle GTWR models in each completed dataset and combines estimates using Rubin's rules.

The primary GTWR outcome is HALE. The GTWR distance matrix uses great-circle spatial distance with `longlat = TRUE`, `lambda = 0.5` and `ksi = 0.5`.

## Requirements

The analysis was prepared with R 4.5.1. Install the required packages with:

```r
source("install_dependencies.R")
```

See [DATA.md](DATA.md) for the required input files and directory structure.

## Running the analysis

Open the repository as the working directory and run:

```powershell
& "D:\R\R-4.5.1\bin\Rscript.exe" .\DMI_reproducible_analysis.R
```

On other systems:

```bash
Rscript DMI_reproducible_analysis.R
```

By default, source data are read from `data/raw`. Alternative locations can be supplied without editing the script:

```powershell
$env:GDD_RAW_ROOT = "D:\path\to\raw-data"
$env:GDD_AUX_ROOT = "D:\path\to\auxiliary-data"
```

The formal GTWR run uses 500 uncertainty draws. The number of draws and parallel workers can be configured through `GDD_GTWR_BOOT_B` and `GDD_GTWR_WORKERS`.

## Main outputs

- `multiple_imputation/`: the `mice` object and 20 completed analysis datasets.
- `GTWR_MI20/`: per-imputation GTWR objects and Rubin-pooled local and stratum-level estimates.
- `outputs/health_data/`: reconstructed HALE and LE datasets.
- `outputs/logs/DMI_reproducible_analysis_log.txt`: console output and warnings.
- `global_moran_*.rds` and `local_moran_bootstrap_*.rds`: year-specific spatial results.

## Citation

Please use the citation information in [CITATION.cff](CITATION.cff). The manuscript citation and DOI can be added after publication.

## License

The source code is released under the [MIT License](LICENSE). Input datasets remain subject to the terms of their respective providers.

## Validation statement

The DMI is a research framework whose effectiveness and practical utility require further validation in independent populations and prospective applications.
