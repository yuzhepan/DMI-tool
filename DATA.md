# Data requirements

The paper analysis reads country-level dietary estimates, health estimates, population data and contextual covariates. Place the files under `data/raw` as shown below, or set `GDD_RAW_ROOT` and `GDD_AUX_ROOT` to directories containing the same structure.

```text
data/raw/
├── Country-level estimates/
│   ├── v01_cnty.csv
│   ├── v02_cnty.csv
│   ├── v03_cnty.csv
│   ├── v04_cnty.csv
│   ├── v05_cnty.csv
│   ├── v06_cnty.csv
│   ├── v07_cnty.csv
│   ├── v08_cnty.csv
│   ├── v09_cnty.csv
│   ├── v10_cnty.csv
│   ├── v11_cnty.csv
│   ├── v12_cnty.csv
│   ├── v13_cnty.csv
│   ├── v14_cnty.csv
│   ├── v15_cnty.csv
│   ├── v16_cnty.csv
│   ├── v17_cnty.csv
│   ├── v18_cnty.csv
│   ├── v27_cnty.csv
│   ├── v29_cnty.csv
│   ├── v30_cnty.csv
│   ├── v31_cnty.csv
│   ├── v34_cnty.csv
│   ├── v35_cnty.csv
│   ├── v37_cnty.csv
│   └── v57_cnty.csv
├── GBD2023HALE/
│   ├── HALE2018.csv
│   └── LE2018.csv
├── WPP2022_Population1JanuaryByAge5GroupSex_Medium.csv
├── SDI.xlsx
└── data_qihou.rds
```

## Dietary estimates

The `v*_cnty.csv` files are country-level Global Dietary Database estimates. The analysis uses the country code, age, sex, year, median estimate and 95% uncertainty limits. These source files are distributed separately from this repository.

## Health estimates

`HALE2018.csv` and `LE2018.csv` contain age- and sex-specific healthy life expectancy and life expectancy estimates. The script aggregates them into ages 2–19, 20–64 and 65+ using the WPP medium-variant population data.

## Contextual covariates

- `SDI.xlsx` supplies the Socio-demographic Index.
- `data_qihou.rds` supplies annual temperature, precipitation and elevation variables.
- Urbanization, population density, health expenditure and PM2.5 are retrieved through the World Bank `WDI` package.

## Path overrides

`GDD_RAW_ROOT` is used for the dietary and GBD directories. `GDD_AUX_ROOT` is used for the population, SDI and climate files. When neither variable is set, both default to `data/raw` in the repository.

Example for Windows PowerShell:

```powershell
$env:GDD_RAW_ROOT = "D:\path\to\raw-data"
$env:GDD_AUX_ROOT = "D:\path\to\auxiliary-data"
```

The source datasets remain subject to the access and reuse terms of their respective providers.
