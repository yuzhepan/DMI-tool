# DMI Reproducible Analysis Pipeline
### This script reconstructs the Dietary Matching Index and runs the paper analyses.

# 1. Load packages and resolve function conflicts ================================================
## 1.1 Load packages ================
library(conflicted)
library(countrycode)
library(data.table)
library(dplyr)
library(future)
library(future.apply)
library(GWmodel)
library(here)
library(mice)
library(purrr)
library(readr)
library(readxl)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(sn)
library(sp)
library(spdep)
library(stringr)
library(tibble)
library(tidyr)
library(WDI)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(dplyr::lag)
conflicted::conflicts_prefer(dplyr::select)

## 1.2 Record console output ================
PIPELINE_LOG_DIR <- here::here("outputs", "logs")
dir.create(PIPELINE_LOG_DIR, recursive = TRUE, showWarnings = FALSE)
PIPELINE_LOG_FILE <- file.path(PIPELINE_LOG_DIR, "DMI_reproducible_analysis_log.txt")
pipeline_log_connection <- file(PIPELINE_LOG_FILE, open = "wt")
prior_output_sinks <- sink.number(type = "output")
prior_message_sink <- sink.number(type = "message")
sink(pipeline_log_connection, split = TRUE)
sink(pipeline_log_connection, type = "message")
cat("Pipeline started:", format(Sys.time()), "\n")

tryCatch({

# 2. Construct the Dietary Matching Indices =======================================================
## 2.1 Preprocess country-level dietary data ================
### Define energy references for the three age groups.
GDD_RAW_ROOT <- Sys.getenv("GDD_RAW_ROOT", unset = here::here("data", "raw"))
GDD_AUX_ROOT <- Sys.getenv("GDD_AUX_ROOT", unset = GDD_RAW_ROOT)
COUNTRY_ESTIMATE_DIR <- file.path(GDD_RAW_ROOT, "Country-level estimates")

energy_ref <- tibble::tibble(
  age_group = c("2-19","20-64","65+"),
  energy_kcal = c(1700, 2000, 1700)
)

data_1  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v08_cnty.csv"))
data_2  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v02_cnty.csv"))
data_3  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v01_cnty.csv"))
data_4  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v15_cnty.csv"))
data_5  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v16_cnty.csv"))
data_6  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v06_cnty.csv"))
data_7  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v05_cnty.csv"))
data_8  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v10_cnty.csv"))
data_9  <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v09_cnty.csv"))
data_10 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v11_cnty.csv"))
data_11 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v29_cnty.csv"))
data_12 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v30_cnty.csv"))
data_13 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v31_cnty.csv"))
data_14 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v37_cnty.csv"))

### Inspect the extra 2020 records in the two meat files after the files have been loaded.
table(data_8$year)
table(data_9$year)

clean_data <- function(df){
  df %>%
    filter(female != 999,
           age != 999, age != 0.5, age != 1.5,
           urban == 999, edu == 999) %>%
    select(-urban, -edu)
}

data_list <- list(
  whole_grains    = clean_data(data_1),
  nonstarchy_veg  = clean_data(data_2),
  fruit           = clean_data(data_3),
  ssb             = clean_data(data_4),
  fruit_juice     = clean_data(data_5),
  nuts            = clean_data(data_6),
  legumes         = clean_data(data_7),
  red_meat        = clean_data(data_8),
  processed_meat  = clean_data(data_9),
  seafood         = clean_data(data_10),
  omega6_pctE     = clean_data(data_11),
  seafood_omega3  = clean_data(data_12),   # mg/d -> %E
  plant_omega3    = clean_data(data_13),   # mg/d -> %E
  sodium_mg       = clean_data(data_14)    # mg/d
)

grp_0_19  <- c(3.5, 7.5, 12.5, 17.5)
grp_20_64 <- c(22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5)
grp_65p   <- c(67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5)

fast_group <- function(df){
  df %>%
    filter(female %in% c(0,1), age != 999) %>%
    mutate(age_group = case_when(
      age %in% grp_0_19  ~ "2-19",
      age %in% grp_20_64 ~ "20-64",
      age %in% grp_65p   ~ "65+",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(age_group)) %>%
    mutate(sd_est = (upperci_95 - lowerci_95) / (2*1.96)) %>%
    group_by(iso3, year, female, age_group) %>%
    summarise(
      n_rows = n(),
      k_sd   = sum(is.finite(sd_est) & is.finite(median)),
      mean   = mean(median, na.rm = TRUE),
      se_cl  = ifelse(k_sd > 0,
                      sqrt(sum(sd_est[is.finite(sd_est) & is.finite(median)]^2)) / k_sd,
                      NA_real_),
      se_fallback = {
        m_ok <- median[is.finite(median)]
        if(length(m_ok) > 1) sd(m_ok)/sqrt(length(m_ok)) else NA_real_
      },
      se_use = ifelse(is.finite(se_cl), se_cl, se_fallback),
      ci_lo  = mean - 1.96 * se_use,
      ci_hi  = mean + 1.96 * se_use,
      .groups = "drop"
    ) %>%
    left_join(energy_ref, by = "age_group")
}

processed <- purrr::imap(data_list, ~ fast_group(.x) %>% mutate(var = .y))

keys <- c("iso3","year","female","age_group","energy_kcal")

std_intake <- function(df, stem){
  df %>%
    select(all_of(keys), mean, ci_lo, ci_hi) %>%
    mutate(sd = (ci_hi - ci_lo)/(2*1.96)) %>%
    rename(
      "{stem}"    := mean,
      "{stem}_lo" := ci_lo,
      "{stem}_hi" := ci_hi,
      "{stem}_sd" := sd
    )
}

ahei_whole_grains   <- std_intake(processed$whole_grains,   "whole_grains_g")
ahei_vegetables     <- std_intake(processed$nonstarchy_veg, "vegetables_g")
ahei_fruits         <- std_intake(processed$fruit,          "fruits_g")
ahei_nuts           <- std_intake(processed$nuts,           "nuts_g")
ahei_legumes        <- std_intake(processed$legumes,        "legumes_g")
ahei_red_meat       <- std_intake(processed$red_meat,       "red_meat_g")
ahei_processed_meat <- std_intake(processed$processed_meat, "processed_meat_g")
ahei_sodium_mg      <- std_intake(processed$sodium_mg,      "sodium_mg")

ahei_ssb_juice <- std_intake(processed$ssb, "ssb_g") %>%
  full_join(std_intake(processed$fruit_juice, "juice_g"),
            by = keys) %>%
  mutate(
    ssb_juice_g    = coalesce(ssb_g,0) + coalesce(juice_g,0),
    ssb_juice_sd   = sqrt(coalesce(ssb_g_sd,0)^2 + coalesce(juice_g_sd,0)^2),
    ssb_juice_lo   = ssb_juice_g - 1.96*ssb_juice_sd,
    ssb_juice_hi   = ssb_juice_g + 1.96*ssb_juice_sd
  ) %>%
  select(all_of(keys), ssb_juice_g, ssb_juice_lo, ssb_juice_hi, ssb_juice_sd)

mg_to_pctE <- function(mg, kcal_total){
  ### Convert mg to g, then to kcal (x9), and finally to percentage of energy.
  g    <- mg/1000
  kcal <- g * 9
  100 * kcal / kcal_total
}

ahei_lc_n3_pctE <- std_intake(processed$seafood_omega3, "lc_n3_mg") %>%
  mutate(
    lc_n3_pctE    = mg_to_pctE(lc_n3_mg,    energy_kcal),
    lc_n3_pctE_sd = mg_to_pctE(lc_n3_mg_sd, energy_kcal),
    lc_n3_pctE_lo = lc_n3_pctE - 1.96*lc_n3_pctE_sd,
    lc_n3_pctE_hi = lc_n3_pctE + 1.96*lc_n3_pctE_sd
  ) %>%
  select(all_of(keys), lc_n3_pctE, lc_n3_pctE_lo, lc_n3_pctE_hi, lc_n3_pctE_sd)

ahei_plant_n3_pctE <- std_intake(processed$plant_omega3, "plant_n3_mg") %>%
  mutate(
    plant_n3_pctE    = mg_to_pctE(plant_n3_mg,    energy_kcal),
    plant_n3_pctE_sd = mg_to_pctE(plant_n3_mg_sd, energy_kcal),
    plant_n3_pctE_lo = plant_n3_pctE - 1.96*plant_n3_pctE_sd,
    plant_n3_pctE_hi = plant_n3_pctE + 1.96*plant_n3_pctE_sd
  ) %>%
  select(all_of(keys), plant_n3_pctE, plant_n3_pctE_lo, plant_n3_pctE_hi, plant_n3_pctE_sd)

### PUFA (% energy) is the sum of omega-6 and both omega-3 components.
ahei_pufa_pctE <- std_intake(processed$omega6_pctE, "o6_pctE") %>%
  full_join(ahei_lc_n3_pctE,    by = keys) %>%
  full_join(ahei_plant_n3_pctE, by = keys) %>%
  mutate(
    lc_n3_pctE    = coalesce(lc_n3_pctE, 0),
    plant_n3_pctE = coalesce(plant_n3_pctE, 0),
    lc_n3_pctE_sd    = coalesce(lc_n3_pctE_sd, 0),
    plant_n3_pctE_sd = coalesce(plant_n3_pctE_sd, 0),
    pufa_pctE     = o6_pctE + lc_n3_pctE + plant_n3_pctE,
    pufa_pctE_sd  = sqrt(o6_pctE_sd^2 + lc_n3_pctE_sd^2 + plant_n3_pctE_sd^2),
    pufa_pctE_lo  = pufa_pctE - 1.96*pufa_pctE_sd,
    pufa_pctE_hi  = pufa_pctE + 1.96*pufa_pctE_sd
  ) %>%
  select(all_of(keys), pufa_pctE, pufa_pctE_lo, pufa_pctE_hi, pufa_pctE_sd)

ahei_list <- list(
  ahei_whole_grains, ahei_vegetables, ahei_fruits,
  ahei_nuts, ahei_legumes,
  ahei_red_meat, ahei_processed_meat,
  ahei_sodium_mg,
  ahei_ssb_juice,
  ahei_lc_n3_pctE, ahei_plant_n3_pctE, ahei_pufa_pctE
)

AHEI_wide <- reduce(ahei_list, ~ full_join(.x, .y, by = keys))

AHEI_wide %>% select(iso3, year, female, age_group, energy_kcal) %>% head()
AHEI_wide %>% select(starts_with("pufa_pctE")) %>% summary()
AHEI_wide %>% select(starts_with("ssb_juice_g")) %>% summary()

sapply(AHEI_wide, function(x) sum(is.na(x)))
AHEI_wide <- AHEI_wide[complete.cases(AHEI_wide), ]
length(unique(AHEI_wide$iso3)) # 185

### Additional MSDPS components include eggs, processed meat and unprocessed red meat.
## 2.2 Construct MSDPS component data ================
data_15 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v03_cnty.csv"))
data_16 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v12_cnty.csv"))
data_17 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v57_cnty.csv"))
data_18 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v14_cnty.csv"))
data_19 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v13_cnty.csv"))
data_20 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v04_cnty.csv"))

processed_new <- list(
  potatoes      = clean_data(data_15)      |> fast_group() |> dplyr::mutate(var = "potatoes"),
  eggs          = clean_data(data_16)      |> fast_group() |> dplyr::mutate(var = "eggs"),
  milk          = clean_data(data_17)      |> fast_group() |> dplyr::mutate(var = "milk"),
  yoghurt       = clean_data(data_18)      |> fast_group() |> dplyr::mutate(var = "yoghurt"),
  cheese        = clean_data(data_19)      |> fast_group() |> dplyr::mutate(var = "cheese"),
  other_starchy = clean_data(data_20)      |> fast_group() |> dplyr::mutate(var = "other_starchy")
)


msdps_whole_grains   <- std_intake(processed$whole_grains,   "whole_grains_g")
msdps_veg_nonstarchy <- std_intake(processed$nonstarchy_veg, "vegetables_g")
msdps_fruits         <- std_intake(processed$fruit,          "fruits_g")
msdps_beans_legumes  <- std_intake(processed$legumes,        "legumes_g")
msdps_nuts_seeds     <- std_intake(processed$nuts,           "nuts_g")
msdps_seafood        <- std_intake(processed$seafood,        "seafood_g")
msdps_eggs           <- std_intake(processed_new$eggs,       "eggs_g")

### Dairy equals milk plus yoghurt plus cheese.
msdps_dairy <- std_intake(processed_new$milk, "milk_g") %>%
  dplyr::full_join(std_intake(processed_new$yoghurt, "yoghurt_g"), by = keys) %>%
  dplyr::full_join(std_intake(processed_new$cheese,  "cheese_g"),  by = keys) %>%
  dplyr::mutate(
    dairy_g  = dplyr::coalesce(milk_g,0) + dplyr::coalesce(yoghurt_g,0) + dplyr::coalesce(cheese_g,0),
    dairy_sd = sqrt(dplyr::coalesce(milk_g_sd,0)^2 + dplyr::coalesce(yoghurt_g_sd,0)^2 + dplyr::coalesce(cheese_g_sd,0)^2),
    dairy_lo = dairy_g - 1.96*dairy_sd,
    dairy_hi = dairy_g + 1.96*dairy_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), dairy_g, dairy_lo, dairy_hi, dairy_sd)

### Starchy vegetables equal potatoes plus other starchy vegetables.
msdps_starchy <- std_intake(processed_new$potatoes,      "potatoes_g") %>%
  dplyr::full_join(std_intake(processed_new$other_starchy, "other_starchy_g"), by = keys) %>%
  dplyr::mutate(
    starchy_g  = dplyr::coalesce(potatoes_g,0) + dplyr::coalesce(other_starchy_g,0),
    starchy_sd = sqrt(dplyr::coalesce(potatoes_g_sd,0)^2 + dplyr::coalesce(other_starchy_g_sd,0)^2),
    starchy_lo = starchy_g - 1.96*starchy_sd,
    starchy_hi = starchy_g + 1.96*starchy_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), starchy_g, starchy_lo, starchy_hi, starchy_sd)

### Total meat equals red meat plus processed meat.
msdps_meats <- std_intake(processed$red_meat,       "red_meat_g") %>%
  dplyr::full_join(std_intake(processed$processed_meat, "processed_meat_g"), by = keys) %>%
  dplyr::mutate(
    meats_g  = dplyr::coalesce(red_meat_g,0) + dplyr::coalesce(processed_meat_g,0),
    meats_sd = sqrt(dplyr::coalesce(red_meat_g_sd,0)^2 + dplyr::coalesce(processed_meat_g_sd,0)^2),
    meats_lo = meats_g - 1.96*meats_sd,
    meats_hi = meats_g + 1.96*meats_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), meats_g, meats_lo, meats_hi, meats_sd)

msdps_list <- list(
  msdps_whole_grains,
  msdps_veg_nonstarchy,
  msdps_fruits,
  msdps_dairy,
  msdps_seafood,
  msdps_beans_legumes,
  msdps_nuts_seeds,
  msdps_starchy,
  msdps_eggs,
  msdps_meats
)

MSDPS_wide <- purrr::reduce(msdps_list, ~ dplyr::full_join(.x, .y, by = keys))

sapply(MSDPS_wide, function(x) sum(is.na(x)))
MSDPS_wide <- MSDPS_wide[complete.cases(MSDPS_wide), ]
length(unique(MSDPS_wide$iso3)) # 185

### The remaining hPDI components cover beverages, dairy, eggs, seafood and meat.

## 2.3 Construct hPDI component data ================
data_21 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v07_cnty.csv"))
data_22 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v18_cnty.csv"))
data_23 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v17_cnty.csv"))


processed_hpdi <- list(
  refined_grains = clean_data(data_21) |> fast_group() |> dplyr::mutate(var = "refined_grains"),
  tea_cups       = clean_data(data_22) |> fast_group() |> dplyr::mutate(var = "tea_cups"),
  coffee_cups    = clean_data(data_23) |> fast_group() |> dplyr::mutate(var = "coffee_cups")
)

hpdi_whole_grains   <- std_intake(processed$whole_grains,   "whole_grains_g")
hpdi_veg_nonstarchy <- std_intake(processed$nonstarchy_veg, "vegetables_g")
hpdi_fruits         <- std_intake(processed$fruit,          "fruits_g")
hpdi_nuts_seeds     <- std_intake(processed$nuts,           "nuts_g")
hpdi_beans_legumes  <- std_intake(processed$legumes,        "legumes_g")
hpdi_ssb            <- std_intake(processed$ssb,            "ssb_g")
hpdi_fruit_juice    <- std_intake(processed$fruit_juice,    "fruit_juice_g")
hpdi_seafood        <- std_intake(processed$seafood,        "seafood_g")
hpdi_red_meat       <- std_intake(processed$red_meat,       "red_meat_g")
hpdi_proc_meat      <- std_intake(processed$processed_meat, "processed_meat_g")
hpdi_potatoes       <- std_intake(processed_new$potatoes, "potatoes_g")
hpdi_refined_grains <- std_intake(processed_hpdi$refined_grains, "refined_grains_g")
hpdi_eggs           <- std_intake(processed_new$eggs, "eggs_g")


### Convert tea and coffee from cups/day to g/day (248 g/cup for coffee; 237 g/cup for tea).
cups_to_g <- function(x, grams_per_cup) x * grams_per_cup

hpdi_tea_g <- std_intake(processed_hpdi$tea_cups, "tea_cups") %>%
  dplyr::mutate(
    tea_g    = cups_to_g(tea_cups,    237),
    tea_g_sd = cups_to_g(tea_cups_sd, 237),
    tea_g_lo = tea_g - 1.96*tea_g_sd,
    tea_g_hi = tea_g + 1.96*tea_g_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), tea_g, tea_g_lo, tea_g_hi, tea_g_sd)

hpdi_coffee_g <- std_intake(processed_hpdi$coffee_cups, "coffee_cups") %>%
  dplyr::mutate(
    coffee_g    = cups_to_g(coffee_cups,    248),
    coffee_g_sd = cups_to_g(coffee_cups_sd, 248),
    coffee_g_lo = coffee_g - 1.96*coffee_g_sd,
    coffee_g_hi = coffee_g + 1.96*coffee_g_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), coffee_g, coffee_g_lo, coffee_g_hi, coffee_g_sd)

hpdi_tea_coffee <- hpdi_tea_g %>%
  dplyr::full_join(hpdi_coffee_g, by = keys) %>%
  dplyr::mutate(
    tea_coffee_g  = dplyr::coalesce(tea_g,0) + dplyr::coalesce(coffee_g,0),
    tea_coffee_sd = sqrt(dplyr::coalesce(tea_g_sd,0)^2 + dplyr::coalesce(coffee_g_sd,0)^2),
    tea_coffee_lo = tea_coffee_g - 1.96*tea_coffee_sd,
    tea_coffee_hi = tea_coffee_g + 1.96*tea_coffee_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), tea_coffee_g, tea_coffee_lo, tea_coffee_hi, tea_coffee_sd)

### Dairy equals milk plus yoghurt plus cheese (g/day).
hpdi_dairy <- std_intake(processed_new$milk, "milk_g") %>%
  dplyr::full_join(std_intake(processed_new$yoghurt, "yoghurt_g"), by = keys) %>%
  dplyr::full_join(std_intake(processed_new$cheese,  "cheese_g"),  by = keys) %>%
  dplyr::mutate(
    dairy_g  = dplyr::coalesce(milk_g,0) + dplyr::coalesce(yoghurt_g,0) + dplyr::coalesce(cheese_g,0),
    dairy_sd = sqrt(dplyr::coalesce(milk_g_sd,0)^2 + dplyr::coalesce(yoghurt_g_sd,0)^2 + dplyr::coalesce(cheese_g_sd,0)^2),
    dairy_lo = dairy_g - 1.96*dairy_sd,
    dairy_hi = dairy_g + 1.96*dairy_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), dairy_g, dairy_lo, dairy_hi, dairy_sd)

hpdi_list <- list(
  hpdi_whole_grains,
  hpdi_veg_nonstarchy,
  hpdi_fruits,
  hpdi_nuts_seeds,
  hpdi_beans_legumes,
  hpdi_tea_coffee,
  hpdi_fruit_juice,
  hpdi_refined_grains,
  hpdi_potatoes,
  hpdi_ssb,
  hpdi_dairy,
  hpdi_eggs,
  hpdi_seafood,
  hpdi_red_meat,
  hpdi_proc_meat
)

hPDI_wide <- purrr::reduce(hpdi_list, ~ dplyr::full_join(.x, .y, by = keys))

sapply(hPDI_wide, function(x) sum(is.na(x)))
hPDI_wide <- hPDI_wide[complete.cases(hPDI_wide), ]
length(unique(hPDI_wide$iso3)) # 185


dash_whole_grains <- std_intake(processed$whole_grains,   "whole_grains_g")
dash_vegetables   <- std_intake(processed$nonstarchy_veg, "vegetables_g")
dash_fruits       <- std_intake(processed$fruit,          "fruits_g")
dash_sodium       <- std_intake(processed$sodium_mg,      "sodium_mg")
dash_ssb          <- std_intake(processed$ssb,            "ssb_g")

dash_nuts_legumes <- std_intake(processed$nuts, "nuts_g") %>%
  dplyr::full_join(std_intake(processed$legumes, "legumes_g"), by = keys) %>%
  dplyr::mutate(
    nuts_legumes_g  = dplyr::coalesce(nuts_g, 0) + dplyr::coalesce(legumes_g, 0),
    nuts_legumes_sd = sqrt(dplyr::coalesce(nuts_g_sd, 0)^2 + dplyr::coalesce(legumes_g_sd, 0)^2),
    nuts_legumes_lo = nuts_legumes_g - 1.96 * nuts_legumes_sd,
    nuts_legumes_hi = nuts_legumes_g + 1.96 * nuts_legumes_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), nuts_legumes_g, nuts_legumes_lo, nuts_legumes_hi, nuts_legumes_sd)

### Dairy equals milk plus yoghurt plus cheese (g/day).
dash_dairy <- std_intake(processed_new$milk, "milk_g") %>%
  dplyr::full_join(std_intake(processed_new$yoghurt, "yoghurt_g"), by = keys) %>%
  dplyr::full_join(std_intake(processed_new$cheese,  "cheese_g"),  by = keys) %>%
  dplyr::mutate(
    dairy_g  = dplyr::coalesce(milk_g, 0) + dplyr::coalesce(yoghurt_g, 0) + dplyr::coalesce(cheese_g, 0),
    dairy_sd = sqrt(dplyr::coalesce(milk_g_sd, 0)^2 +
                      dplyr::coalesce(yoghurt_g_sd, 0)^2 +
                      dplyr::coalesce(cheese_g_sd, 0)^2),
    dairy_lo = dairy_g - 1.96 * dairy_sd,
    dairy_hi = dairy_g + 1.96 * dairy_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), dairy_g, dairy_lo, dairy_hi, dairy_sd)

### Total meat equals processed plus unprocessed red meat (g/day).
dash_meats <- std_intake(processed$processed_meat, "processed_meat_g") %>%
  dplyr::full_join(std_intake(processed$red_meat, "red_meat_g"), by = keys) %>%
  dplyr::mutate(
    meats_g  = dplyr::coalesce(processed_meat_g, 0) + dplyr::coalesce(red_meat_g, 0),
    meats_sd = sqrt(dplyr::coalesce(processed_meat_g_sd, 0)^2 +
                      dplyr::coalesce(red_meat_g_sd, 0)^2),
    meats_lo = meats_g - 1.96 * meats_sd,
    meats_hi = meats_g + 1.96 * meats_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), meats_g, meats_lo, meats_hi, meats_sd)

dash_list <- list(
  dash_whole_grains, dash_vegetables, dash_fruits,
  dash_nuts_legumes, dash_dairy, dash_sodium, dash_ssb, dash_meats
)

DASH_wide <- purrr::reduce(dash_list, ~ dplyr::full_join(.x, .y, by = keys))

sapply(DASH_wide , function(x) sum(is.na(x)))
DASH_wide  <- DASH_wide [complete.cases(DASH_wide ), ]
length(unique(DASH_wide$iso3)) # 185

## 2.4 Construct DASH and DRRD component data ================
data_24 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v27_cnty.csv"))

processed_drrd <- list(
  sfa_pctE = clean_data(data_24) |> fast_group() |> dplyr::mutate(var = "sfa_pctE")
)

### Omega-6 is already expressed as a percentage of energy.
o6_pctE <- std_intake(processed$omega6_pctE, "o6_pctE")

lc_n3_pctE <- std_intake(processed$seafood_omega3, "lc_n3_mg") %>%
  dplyr::mutate(
    lc_n3_pctE    = mg_to_pctE(lc_n3_mg,    energy_kcal),
    lc_n3_pctE_sd = mg_to_pctE(lc_n3_mg_sd, energy_kcal),
    lc_n3_pctE_lo = lc_n3_pctE - 1.96*lc_n3_pctE_sd,
    lc_n3_pctE_hi = lc_n3_pctE + 1.96*lc_n3_pctE_sd
  ) %>% dplyr::select(dplyr::all_of(keys), starts_with("lc_n3_pctE"))

plant_n3_pctE <- std_intake(processed$plant_omega3, "plant_n3_mg") %>%
  dplyr::mutate(
    plant_n3_pctE    = mg_to_pctE(plant_n3_mg,    energy_kcal),
    plant_n3_pctE_sd = mg_to_pctE(plant_n3_mg_sd, energy_kcal),
    plant_n3_pctE_lo = plant_n3_pctE - 1.96*plant_n3_pctE_sd,
    plant_n3_pctE_hi = plant_n3_pctE + 1.96*plant_n3_pctE_sd
  ) %>% dplyr::select(dplyr::all_of(keys), starts_with("plant_n3_pctE"))

drrd_pufa_pctE <- o6_pctE %>%
  dplyr::full_join(lc_n3_pctE,    by = keys) %>%
  dplyr::full_join(plant_n3_pctE, by = keys) %>%
  dplyr::mutate(
    across(c(lc_n3_pctE, plant_n3_pctE, lc_n3_pctE_sd, plant_n3_pctE_sd), ~ dplyr::coalesce(.x, 0)),
    pufa_pctE    = o6_pctE + lc_n3_pctE + plant_n3_pctE,
    pufa_pctE_sd = sqrt(o6_pctE_sd^2 + lc_n3_pctE_sd^2 + plant_n3_pctE_sd^2),
    pufa_pctE_lo = pufa_pctE - 1.96*pufa_pctE_sd,
    pufa_pctE_hi = pufa_pctE + 1.96*pufa_pctE_sd
  ) %>% dplyr::select(dplyr::all_of(keys), pufa_pctE, pufa_pctE_lo, pufa_pctE_hi, pufa_pctE_sd)

### Saturated fat is expressed as a percentage of energy.
drrd_sfa_pctE <- std_intake(processed_drrd$sfa_pctE, "sfa_pctE")

### Approximate the ratio SD with the delta method.
drrd_ps_ratio <- drrd_pufa_pctE %>%
  dplyr::full_join(drrd_sfa_pctE, by = keys) %>%
  dplyr::mutate(
    ps_ratio   = pufa_pctE / sfa_pctE,
    sd_ratio   = ps_ratio * sqrt( (pmax(pufa_pctE_sd, 1e-12)/pmax(pufa_pctE, 1e-12))^2 +
                                    (pmax(sfa_pctE_sd,  1e-12)/pmax(sfa_pctE,  1e-12))^2 ),
    ps_ratio_lo = ps_ratio - 1.96*sd_ratio,
    ps_ratio_hi = ps_ratio + 1.96*sd_ratio
  ) %>%
  dplyr::select(dplyr::all_of(keys), ps_ratio, ps_ratio_lo, ps_ratio_hi, sd_ratio)

### Convert coffee from cups/day to g/day using 248 g per cup.
drrd_coffee_g <- std_intake(processed_hpdi$coffee_cups, "coffee_cups") %>%
  dplyr::mutate(
    coffee_g    = cups_to_g(coffee_cups,    248),
    coffee_g_sd = cups_to_g(coffee_cups_sd, 248),
    coffee_g_lo = coffee_g - 1.96*coffee_g_sd,
    coffee_g_hi = coffee_g + 1.96*coffee_g_sd
  ) %>% dplyr::select(dplyr::all_of(keys), coffee_g, coffee_g_lo, coffee_g_hi, coffee_g_sd)

### Retain nuts, seeds and fruits in g/day.
drrd_nuts_g   <- std_intake(processed$nuts,  "nuts_g")
drrd_fruits_g <- std_intake(processed$fruit, "fruits_g")

drrd_ssb_juice <- std_intake(processed$ssb, "ssb_g") %>%
  dplyr::full_join(std_intake(processed$fruit_juice, "juice_g"), by = keys) %>%
  dplyr::mutate(
    ssb_juice_g  = dplyr::coalesce(ssb_g,0) + dplyr::coalesce(juice_g,0),
    ssb_juice_sd = sqrt(dplyr::coalesce(ssb_g_sd,0)^2 + dplyr::coalesce(juice_g_sd,0)^2),
    ssb_juice_lo = ssb_juice_g - 1.96*ssb_juice_sd,
    ssb_juice_hi = ssb_juice_g + 1.96*ssb_juice_sd
  ) %>% dplyr::select(dplyr::all_of(keys), ssb_juice_g, ssb_juice_lo, ssb_juice_hi, ssb_juice_sd)

### Combine processed and unprocessed red meat in g/day.
drrd_meats_g <- std_intake(processed$red_meat, "red_meat_g") %>%
  dplyr::full_join(std_intake(processed$processed_meat, "processed_meat_g"), by = keys) %>%
  dplyr::mutate(
    meats_g  = dplyr::coalesce(red_meat_g,0) + dplyr::coalesce(processed_meat_g,0),
    meats_sd = sqrt(dplyr::coalesce(red_meat_g_sd,0)^2 + dplyr::coalesce(processed_meat_g_sd,0)^2),
    meats_lo = meats_g - 1.96*meats_sd,
    meats_hi = meats_g + 1.96*meats_sd
  ) %>% dplyr::select(dplyr::all_of(keys), meats_g, meats_lo, meats_hi, meats_sd)

drrd_list <- list(
  drrd_ps_ratio,
  drrd_coffee_g,
  drrd_nuts_g,
  drrd_fruits_g,
  drrd_ssb_juice,
  drrd_meats_g
)
DRRD_wide <- purrr::reduce(drrd_list, ~ dplyr::full_join(.x, .y, by = keys))

sapply(DRRD_wide, function(x) sum(is.na(x)))
DRRD_wide <- DRRD_wide[complete.cases(DRRD_wide), ]
length(unique(DRRD_wide$iso3)) # 185


## 2.5 Construct WCRF component data ================
data_25 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v34_cnty.csv"))

processed_wcrf <- list(
  fiber = clean_data(data_25) |> fast_group() |> dplyr::mutate(var = "fiber")
)
wcrf_fiber <- std_intake(processed_wcrf$fiber, "fiber_g")

wcrf_fruit_veg <- std_intake(processed$fruit,          "fruits_g") %>%
  dplyr::full_join(std_intake(processed$nonstarchy_veg, "vegetables_g"), by = keys) %>%
  dplyr::mutate(
    fruits_veg_g  = dplyr::coalesce(fruits_g, 0) + dplyr::coalesce(vegetables_g, 0),
    fruits_veg_sd = sqrt(dplyr::coalesce(fruits_g_sd, 0)^2 + dplyr::coalesce(vegetables_g_sd, 0)^2),
    fruits_veg_lo = fruits_veg_g - 1.96 * fruits_veg_sd,
    fruits_veg_hi = fruits_veg_g + 1.96 * fruits_veg_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), fruits_veg_g, fruits_veg_lo, fruits_veg_hi, fruits_veg_sd)

wcrf_ssb          <- std_intake(processed$ssb,            "ssb_g")
wcrf_proc_meat    <- std_intake(processed$processed_meat, "processed_meat_g")
wcrf_red_meat     <- std_intake(processed$red_meat,       "red_meat_g")

wcrf_list <- list(
  wcrf_fruit_veg,
  wcrf_fiber,
  wcrf_ssb,
  wcrf_proc_meat,
  wcrf_red_meat
)

WCRF_wide <- purrr::reduce(wcrf_list, ~ dplyr::full_join(.x, .y, by = keys))

sapply(WCRF_wide, function(x) sum(is.na(x)))
WCRF_wide <- WCRF_wide[complete.cases(WCRF_wide), ]
length(unique(WCRF_wide$iso3)) # 185

### PHDI additionally includes eggs, seafood, nuts, legumes and added sugar.
## 2.6 Construct PHDI and DBI component data ================
data_26 <- read.csv(file.path(COUNTRY_ESTIMATE_DIR, "v35_cnty.csv"))

processed_phdi <- list(
  added_sugar = clean_data(data_26) |> fast_group() |> dplyr::mutate(var = "added_sugar")
)

phdi_added_sugar <- std_intake(processed_phdi$added_sugar, "added_sugar_g")

### Dairy equals milk plus yoghurt plus cheese (g/day).
phdi_dairy <- std_intake(processed_new$milk, "milk_g") %>%
  dplyr::full_join(std_intake(processed_new$yoghurt, "yoghurt_g"), by = keys) %>%
  dplyr::full_join(std_intake(processed_new$cheese,  "cheese_g"),  by = keys) %>%
  dplyr::mutate(
    dairy_g  = dplyr::coalesce(milk_g, 0) + dplyr::coalesce(yoghurt_g, 0) + dplyr::coalesce(cheese_g, 0),
    dairy_sd = sqrt(dplyr::coalesce(milk_g_sd, 0)^2 +
                      dplyr::coalesce(yoghurt_g_sd, 0)^2 +
                      dplyr::coalesce(cheese_g_sd, 0)^2),
    dairy_lo = dairy_g - 1.96 * dairy_sd,
    dairy_hi = dairy_g + 1.96 * dairy_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), dairy_g, dairy_lo, dairy_hi, dairy_sd)

phdi_starchy <- std_intake(processed_new$potatoes, "potatoes_g") %>%
  dplyr::full_join(std_intake(processed_new$other_starchy, "other_starchy_g"), by = keys) %>%
  dplyr::mutate(
    starchy_g  = dplyr::coalesce(potatoes_g, 0) + dplyr::coalesce(other_starchy_g, 0),
    starchy_sd = sqrt(dplyr::coalesce(potatoes_g_sd, 0)^2 + dplyr::coalesce(other_starchy_g_sd, 0)^2),
    starchy_lo = starchy_g - 1.96 * starchy_sd,
    starchy_hi = starchy_g + 1.96 * starchy_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), starchy_g, starchy_lo, starchy_hi, starchy_sd)

phdi_whole_grains   <- std_intake(processed$whole_grains,   "whole_grains_g")
phdi_veg_nonstarchy <- std_intake(processed$nonstarchy_veg, "vegetables_g")
phdi_fruits         <- std_intake(processed$fruit,          "fruits_g")
phdi_eggs           <- std_intake(processed_new$eggs,           "eggs_g")
phdi_red_meat       <- std_intake(processed$red_meat,       "red_meat_g")
phdi_proc_meat      <- std_intake(processed$processed_meat, "processed_meat_g")
phdi_seafood        <- std_intake(processed$seafood,        "seafood_g")
phdi_nuts_seeds     <- std_intake(processed$nuts,           "nuts_g")
phdi_beans_legumes  <- std_intake(processed$legumes,        "legumes_g")

phdi_list <- list(
  phdi_whole_grains,
  phdi_starchy,          # potatoes + other starchy
  phdi_veg_nonstarchy,
  phdi_fruits,
  phdi_dairy,            # milk + yoghurt + cheese
  phdi_eggs,
  phdi_red_meat,
  phdi_proc_meat,
  phdi_seafood,
  phdi_nuts_seeds,
  phdi_beans_legumes,
  phdi_added_sugar       # added sugar
)

PHDI_wide <- purrr::reduce(phdi_list, ~ dplyr::full_join(.x, .y, by = keys))

sapply(PHDI_wide , function(x) sum(is.na(x)))
PHDI_wide <- PHDI_wide[complete.cases(PHDI_wide), ]
length(unique(PHDI_wide$iso3)) # 185

dbi_whole_grains <- std_intake(processed$whole_grains,   "whole_grains_g")
dbi_vegetables   <- std_intake(processed$nonstarchy_veg, "vegetables_g")
dbi_fruits       <- std_intake(processed$fruit,          "fruits_g")
dbi_red_meat     <- std_intake(processed$red_meat,       "red_meat_g")
dbi_seafood      <- std_intake(processed$seafood,        "seafood_g")
dbi_eggs         <- std_intake(processed_new$eggs,           "eggs_g")
dbi_sodium       <- std_intake(processed$sodium_mg,         "sodium_mg")

### Dairy equals milk plus yoghurt plus cheese (g/day).
dbi_dairy <- std_intake(processed_new$milk, "milk_g") %>%
  dplyr::full_join(std_intake(processed_new$yoghurt, "yoghurt_g"), by = keys) %>%
  dplyr::full_join(std_intake(processed_new$cheese,  "cheese_g"),  by = keys) %>%
  dplyr::mutate(
    dairy_g  = dplyr::coalesce(milk_g, 0) + dplyr::coalesce(yoghurt_g, 0) + dplyr::coalesce(cheese_g, 0),
    dairy_sd = sqrt(dplyr::coalesce(milk_g_sd, 0)^2 +
                      dplyr::coalesce(yoghurt_g_sd, 0)^2 +
                      dplyr::coalesce(cheese_g_sd, 0)^2),
    dairy_lo = dairy_g - 1.96 * dairy_sd,
    dairy_hi = dairy_g + 1.96 * dairy_sd
  ) %>%
  dplyr::select(dplyr::all_of(keys), dairy_g, dairy_lo, dairy_hi, dairy_sd)

dbi_list <- list(
  dbi_whole_grains,
  dbi_vegetables,
  dbi_fruits,
  dbi_dairy,
  dbi_red_meat,
  dbi_seafood,
  dbi_eggs,
  dbi_sodium
)

DBI_wide <- purrr::reduce(dbi_list, ~ dplyr::full_join(.x, .y, by = keys))
sapply(DBI_wide , function(x) sum(is.na(x)))
DBI_wide <- DBI_wide[complete.cases(DBI_wide), ]
length(unique(DBI_wide$iso3)) # 185

## 2.7 Enforce nonnegative dietary values ================
sanitize_nonneg <- function(df, eps = 1e-6) {
  df %>%
    mutate(across(
      where(is.numeric),
      ~ ifelse(. < 0, eps, .)
    ))
}

AHEI_wide  <- sanitize_nonneg(AHEI_wide)
MSDPS_wide <- sanitize_nonneg(MSDPS_wide)
hPDI_wide  <- sanitize_nonneg(hPDI_wide)
DASH_wide  <- sanitize_nonneg(DASH_wide)
DRRD_wide  <- sanitize_nonneg(DRRD_wide)
WCRF_wide  <- sanitize_nonneg(WCRF_wide)
PHDI_wide  <- sanitize_nonneg(PHDI_wide)
DBI_wide   <- sanitize_nonneg(DBI_wide)

check_neg <- function(df) sum(df[sapply(df, is.numeric)] < 0, na.rm = TRUE)
data_frames <- list(
  AHEI = AHEI_wide,
  MSDPS = MSDPS_wide,
  hPDI = hPDI_wide,
  DASH = DASH_wide,
  DRRD = DRRD_wide,
  WCRF = WCRF_wide,
  PHDI = PHDI_wide,
  DBI = DBI_wide
)

sapply(data_frames, check_neg)

## 2.8 Define the probability-based DMI algorithm ================
.impute_U <- function(W0, T0, direction, kappa = 0.2){
  if (is.na(W0) | is.na(T0)) return(NA_real_)
  d <- tolower(direction)

  if (d == "positive") {
    return(T0 * (1 + kappa))
  } else if (d == "negative") {
    U0 <- T0 * (1 - kappa)
    return(max(U0, 0))
  } else {
    return(NA_real_)
  }
}
.validate_WTU <- function(W, T, U){
  low  <- min(W, U)
  high <- max(W, U)
  is.finite(W) && is.finite(T) && is.finite(U) && (T >= low) && (T <= high)
}

calc_prob <- function(x, W, T, U, direction = "positive", alpha0 = 3, scale_div = 4){
  if (any(is.na(c(x,W,T,U)))) return(NA_real_)
  if (!.validate_WTU(W,T,U))  return(NA_real_)
  xi    <- T
  omega <- max(abs(U - W) / scale_div, 1e-6)
  alpha <- ifelse(tolower(direction)=="positive", -abs(alpha0), abs(alpha0))
  cdf   <- psn(x, xi, omega, alpha)
  p     <- if (tolower(direction)=="positive") cdf else (1 - cdf)
  pmin(pmax(p, 1e-6), 1 - 1e-6)
}

power_mean_w <- function(vals, w = NULL, rho = 0.3, eps = 1e-12){
  keep <- is.finite(vals)
  vals <- vals[keep]
  if (!length(vals)) return(NA_real_)
  if (is.null(w)) w <- rep(1/length(vals), length(vals)) else {
    w <- w[keep]
    sw <- sum(w); if (sw <= 0) return(NA_real_)
    w <- w / sw
  }
  vals <- pmin(pmax(vals, eps), 1 - eps)

  if (abs(rho) < 1e-6) {
    return( exp(sum(w * log(vals))) )
  } else {
    ### Evaluate the generalized mean on the log scale for numerical stability.
    z <- rho * log(vals)
    m <- max(z)
    s <- sum(w * exp(z - m))
    return( exp((log(s) + m) / rho) )
  }
}

build_DMI <- function(data, ref_table,
                      rho = 0.3,
                      rescale_mode = c("pattern","none"),
                      kappa = 0.2,
                      alpha0 = 3, scale_div = 4,
                      base_kcal = 2000,
                      use_max_pts = FALSE,
                      return_probs = TRUE) {

  rescale_mode <- match.arg(rescale_mode)
  ref_use <- ref_table %>% dplyr::filter(col %in% names(data))
  food_cols <- unique(ref_use$col)

  out <- data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      .probs_vec = list({
        probs <- rep(NA_real_, nrow(ref_use))
        wvec  <- rep(NA_real_, nrow(ref_use))
        for (i in seq_len(nrow(ref_use))) {
          r <- ref_use[i, ]
          if (!r$sex %in% c("both","F","M")) next
          if (r$sex == "F" && female != 1) next
          if (r$sex == "M" && female != 0) next
          U0i <- ifelse(is.na(r$U0), .impute_U(r$W0, r$T0, r$direction, kappa), r$U0)
          scale_factor <- if (isTRUE(r$energy_scaled)) energy_kcal / base_kcal else 1
          W <- r$W0 * scale_factor; T <- r$T0 * scale_factor; U <- U0i * scale_factor
          x <- get(r$col)
          probs[i] <- calc_prob(x, W, T, U, r$direction, alpha0 = alpha0, scale_div = scale_div)
          wvec[i]  <- if (isTRUE(use_max_pts)) r$max_pts else 1
        }
        p_by_food <- tapply(probs, ref_use$col, function(v) {
          v_f <- v[is.finite(v)]
          if (length(v_f)) v_f[1] else NA_real_
        })
        p_by_food <- p_by_food[food_cols]
        stats::setNames(as.numeric(p_by_food), food_cols)
      }),
      DMI_raw = {
        probs <- rep(NA_real_, nrow(ref_use))
        wvec  <- rep(NA_real_, nrow(ref_use))
        for (i in seq_len(nrow(ref_use))) {
          r <- ref_use[i, ]
          if (!r$sex %in% c("both","F","M")) next
          if (r$sex == "F" && female != 1) next
          if (r$sex == "M" && female != 0) next
          U0i <- ifelse(is.na(r$U0), .impute_U(r$W0, r$T0, r$direction, kappa), r$U0)
          scale_factor <- if (isTRUE(r$energy_scaled)) energy_kcal / base_kcal else 1
          W <- r$W0 * scale_factor; T <- r$T0 * scale_factor; U <- U0i * scale_factor
          x <- get(r$col)
          probs[i] <- calc_prob(x, W, T, U, r$direction, alpha0 = alpha0, scale_div = scale_div)
          wvec[i]  <- if (isTRUE(use_max_pts)) r$max_pts else 1
        }
        power_mean_w(probs, w = wvec, rho = rho)
      }
    ) %>% dplyr::ungroup()

  out <- if (rescale_mode == "pattern") {
    out %>% dplyr::mutate(DMI = scales::rescale(DMI_raw, to = c(0,100), na.rm = TRUE))
  } else {
    out %>% dplyr::mutate(DMI = DMI_raw * 100)
  }

  if (!return_probs) return(out)

  prob_by_group_wide <- out %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::select(row_id, .probs_vec) %>%
    tidyr::unnest_wider(.probs_vec, names_repair = "check_unique")

  prob_by_group_long <- prob_by_group_wide %>%
    tidyr::pivot_longer(cols = all_of(food_cols),
                        names_to = "food_group", values_to = "prob")

  list(
    out = out,
    prob_by_group_wide = prob_by_group_wide,
    prob_by_group_long = prob_by_group_long,
    food_groups = food_cols
  )
}

## 2.9 Define pattern-specific reference tables and calculate DMI values ================
### 2.9.1 AHEI
ref_AHEI <- tibble::tribble(
  ~col,               ~direction, ~W0,   ~T0,   ~U0,  ~energy_scaled, ~sex,   ~max_pts,
  "whole_grains_g",   "positive",   0,    75,   NA,   TRUE,           "F",    10,
  "whole_grains_g",   "positive",   0,    90,   NA,   TRUE,           "M",    10,
  "vegetables_g",     "positive",   0,   600,   NA,   TRUE,           "both", 10,
  "fruits_g",         "positive",   0,   472,   NA,   TRUE,           "both", 10,
  "nuts_g",           "positive",   0,  28.35,   NA,   TRUE,           "both", 10,
  "legumes_g",        "positive",   0,  28.35,   NA,   TRUE,           "both", 10,
  "red_meat_g",       "negative",   169.5, 0,   NA,   TRUE,           "both", 10,
  "processed_meat_g", "negative",   64.5,  0,   NA,   TRUE,           "both", 10,
  "ssb_juice_g",      "negative",  248,    0,   NA,   TRUE,           "both", 10,
  "sodium_mg",        "negative", 3337, 1112,  NA,   FALSE,           "F",    10,
  "sodium_mg",        "negative", 5271, 1612,  NA,   FALSE,           "M",    10,
  "lc_n3_pctE",       "positive",   0.0, 0.12, NA, FALSE,          "both",  10,
  "pufa_pctE",        "positive",   2.0, 10.0, NA, FALSE,          "both",  10
)
DMI_AHEI_list <- build_DMI(AHEI_wide, ref_AHEI,
                         rho = 0.3, rescale_mode = "none",
                         kappa = 0.2, alpha0 = 3, scale_div = 4,
                         use_max_pts = FALSE)
DMI_AHEI_eq <- DMI_AHEI_list$out
summary(DMI_AHEI_eq$DMI)

DMI_AHEI_list2 <- build_DMI(AHEI_wide, ref_AHEI,
                         rho = 0.3, rescale_mode = "none",
                         kappa = 0.2, alpha0 = 3, scale_div = 4,
                         use_max_pts = TRUE)
DMI_AHEI_wt <- DMI_AHEI_list2$out
summary(DMI_AHEI_wt$DMI)

colnames(MSDPS_wide)
### 2.9.2 MSDPS
ref_MSDPS <- tibble::tribble(
  ~col,               ~direction, ~W0,   ~T0,   ~U0,  ~energy_scaled, ~sex,   ~max_pts,
  "whole_grains_g",   "positive",   0,   120,   NA,   TRUE,           "both", 10,
  "vegetables_g",     "positive",   0,   720,   NA,   TRUE,           "both", 10,
  "fruits_g",         "positive",   0,   354,   NA,   TRUE,           "both", 10,
  "dairy_g",          "positive",   0,   490,   NA,   TRUE,           "both", 10,
  "seafood_g",        "positive",   0,   73,    NA,   TRUE,           "both", 10,
  "legumes_g",        "positive",   0,   16.2,  NA,   TRUE,           "both", 10,
  "nuts_g",           "positive",   0,   16.2,  NA,   TRUE,           "both", 10,
  "starchy_g",        "positive",   0,   68.6,  NA,   TRUE,           "both", 10,
  "eggs_g",           "positive",   0,   23.6,  NA,   TRUE,          "both",  10,
  "meats_g",          "negative",  16.1, 0,     NA,   TRUE,          "both",  10
)
DMI_MSDPS_list <- build_DMI(MSDPS_wide, ref_MSDPS,
                          rho = 0.3, rescale_mode = "none",
                          kappa = 0.2, alpha0 = 3, scale_div = 4,
                          use_max_pts = FALSE)
DMI_MSDPS_eq <- DMI_MSDPS_list$out
summary(DMI_MSDPS_eq$DMI)

DMI_MSDPS_list2 <- build_DMI(MSDPS_wide, ref_MSDPS,
                          rho = 0.3, rescale_mode = "none",
                          kappa = 0.2, alpha0 = 3, scale_div = 4,
                          use_max_pts = TRUE)
DMI_MSDPS_wt <- DMI_MSDPS_list2$out
summary(DMI_MSDPS_wt$DMI)


colnames(hPDI_wide)
### 2.9.3 hPDI
ref_hPDI <- tibble::tribble(
  ~col,               ~direction, ~W0,   ~T0,   ~U0,  ~energy_scaled, ~sex,   ~max_pts,
  "whole_grains_g",   "positive",   4.5,   39,   NA,   TRUE,           "both", 5,
  "vegetables_g",     "positive",   180,   570,   NA,   TRUE,           "both", 5,
  "fruits_g",         "positive",   59,   318.6,   NA,   TRUE,           "both", 5,
  "nuts_g",           "positive",   2.8,  20,   NA,   TRUE,           "both", 5,
  "legumes_g",        "positive",   2.8,  20,   NA,   TRUE,           "both", 5,
  "tea_coffee_g",     "positive",   360,  864,   NA,   TRUE,           "both", 5,
  "fruit_juice_g",    "negative",   198.4,  124,   NA,   TRUE,           "both", 5,
  "refined_grains_g", "negative",   30,   16.5,   NA,   TRUE,          "both",  5,
  "potatoes_g",       "negative",  130.2, 55.8,  NA,   TRUE,            "both", 5,
  "ssb_g",            "negative",  223.2, 24.8,  NA,   TRUE,            "both", 5,
  "dairy_g",          "negative",  563.5, 416.5,   NA,   TRUE,           "both", 5,
  "eggs_g",           "negative",  11,  5.5,   NA,   TRUE,           "both", 5,
  "seafood_g",        "negative",  34,    17,   NA,   TRUE,          "both",  5,
  "red_meat_g",       "negative",  79.1, 33.9,  NA,   TRUE,            "both", 5,
  "processed_meat_g", "negative",  21.5, 8.6,  NA,   TRUE,            "both", 5
)
DMI_hPDI_list <- build_DMI(hPDI_wide, ref_hPDI,
                           rho = 0.3, rescale_mode = "none",
                           kappa = 0.2, alpha0 = 3, scale_div = 4,
                           use_max_pts = FALSE)
DMI_hPDI_eq <- DMI_hPDI_list$out
summary(DMI_hPDI_eq$DMI)

DMI_hPDI_list2 <- build_DMI(hPDI_wide, ref_hPDI,
                            rho = 0.3, rescale_mode = "none",
                            kappa = 0.2, alpha0 = 3, scale_div = 4,
                            use_max_pts = TRUE)
DMI_hPDI_wt <- DMI_hPDI_list2$out
summary(DMI_hPDI_wt$DMI)

colnames(DASH_wide)
### 2.9.4 DASH
ref_DASH <- tibble::tribble(
  ~col,               ~direction, ~W0,   ~T0,   ~U0,  ~energy_scaled, ~sex,   ~max_pts,
  "whole_grains_g",   "positive",   9,    36,   NA,   TRUE,           "both", 5,
  "vegetables_g",     "positive",   264,    552,   NA,   TRUE,        "both", 5,
  "fruits_g",         "positive",   153.4,  448.4,   NA,   TRUE,           "both", 5,
  "nuts_legumes_g",   "positive",   17,     36.9,   NA,   TRUE,           "both", 5,
  "dairy_g",          "positive",   98,     367.5,   NA,   TRUE,           "both", 5,
  "sodium_mg",        "negative",   3500,   2900,   NA,   FALSE,          "both", 5,
  "ssb_g",            "negative",   99.2,   24.8,   NA,   TRUE,           "both", 5,
  "meats_g",          "negative",   117,    39,   NA,   TRUE,           "both", 5
)

DMI_DASH_list <- build_DMI(DASH_wide, ref_DASH,
                           rho = 0.3, rescale_mode = "none",
                           kappa = 0.2, alpha0 = 3, scale_div = 4,
                           use_max_pts = FALSE)
DMI_DASH_eq <- DMI_DASH_list$out
summary(DMI_DASH_eq$DMI)

DMI_DASH_list2 <- build_DMI(DASH_wide, ref_DASH,
                            rho = 0.3, rescale_mode = "none",
                            kappa = 0.2, alpha0 = 3, scale_div = 4,
                            use_max_pts = TRUE)
DMI_DASH_wt <- DMI_DASH_list2$out
summary(DMI_DASH_wt$DMI)

colnames(DRRD_wide)
### 2.9.5 DRRD
ref_DRRD <- tibble::tribble(
  ~col,               ~direction, ~W0,   ~T0,   ~U0,  ~energy_scaled, ~sex,   ~max_pts,
  "ps_ratio",         "positive",   0.5,    0.8,   NA,   FALSE,          "both", 5,
  "coffee_g",         "positive",   297.6,  545.6,   NA,   TRUE,           "both", 5,
  "nuts_g",           "positive",   5.7,    20,   NA,   TRUE,           "both", 5,
  "fruits_g",         "positive",   118,    283.2,   NA,   TRUE,           "both", 5,
  "ssb_juice_g",      "negative",   372,    148.8,   NA,   TRUE,           "both", 5,
  "meats_g",          "negative",   93.6,    31.2,   NA,   TRUE,           "both", 5
)

DMI_DRRD_list <- build_DMI(DRRD_wide, ref_DRRD,
                           rho = 0.3, rescale_mode = "none",
                           kappa = 0.2, alpha0 = 3, scale_div = 4,
                           use_max_pts = FALSE)
DMI_DRRD_eq <- DMI_DRRD_list$out
summary(DMI_DRRD_eq$DMI)

DMI_DRRD_list2 <- build_DMI(DRRD_wide, ref_DRRD,
                            rho = 0.3, rescale_mode = "none",
                            kappa = 0.2, alpha0 = 3, scale_div = 4,
                            use_max_pts = TRUE)
DMI_DRRD_wt <- DMI_DRRD_list2$out
summary(DMI_DRRD_wt$DMI)

colnames(WCRF_wide)
### 2.9.6 WCRF
ref_WCRF <- tibble::tribble(
  ~col,               ~direction, ~W0,   ~T0,   ~U0,  ~energy_scaled, ~sex,   ~max_pts,
  "fruits_veg_g",     "positive",   200,  400,   NA,   TRUE,           "both", 0.5,
  "fiber_g",          "positive",   15,    30,   NA,   TRUE,           "both", 0.5,
  "ssb_g",            "negative",   250,   0,   NA,   TRUE,           "both", 1,
  "processed_meat_g", "negative",   14.2,  0,   NA,   TRUE,           "both", 0.5,
  "red_meat_g",       "negative",   71,    0,   NA,   TRUE,           "both", 0.5
)
DMI_WCRF_list <- build_DMI(WCRF_wide, ref_WCRF,
                           rho = 0.3, rescale_mode = "none",
                           kappa = 0.2, alpha0 = 3, scale_div = 4,
                           use_max_pts = FALSE)
DMI_WCRF_eq <- DMI_WCRF_list$out
summary(DMI_WCRF_eq$DMI)

DMI_WCRF_list2 <- build_DMI(WCRF_wide, ref_WCRF,
                            rho = 0.3, rescale_mode = "none",
                            kappa = 0.2, alpha0 = 3, scale_div = 4,
                            use_max_pts = TRUE)
DMI_WCRF_wt <- DMI_WCRF_list2$out
summary(DMI_WCRF_wt$DMI)

colnames(PHDI_wide)
### 2.9.7 PHDI
ref_PHDI <- tibble::tribble(
  ~col,               ~direction, ~W0,   ~T0,   ~U0,  ~energy_scaled, ~sex,   ~max_pts,
  "whole_grains_g",   "positive",   0,    75,   NA,   TRUE,           "F",    10,
  "whole_grains_g",   "positive",   0,    90,   NA,   TRUE,           "M",    10,
  "starchy_g",        "negative",   200,  50,   NA,   TRUE,           "both", 10,
  "vegetables_g",     "positive",   0,    300,   NA,   TRUE,          "both", 10,
  "fruits_g",         "positive",   0,    200,   NA,   TRUE,          "both", 10,
  "dairy_g",          "negative",   1000, 250,   NA,   TRUE,          "both", 10,
  "eggs_g",           "negative",   120,    13,   NA,   TRUE,         "both", 10,
  "red_meat_g",       "negative",   100,    14,   NA,   TRUE,         "both", 10,
  "processed_meat_g", "negative",   100,    14,   NA,   TRUE,         "both", 10,
  "seafood_g",        "positive",   0,    28,   NA,   TRUE,           "both", 10,
  "nuts_g",           "positive",   0,    50,   NA,   TRUE,           "both", 10,
  "legumes_g",        "positive",   0,    100,   NA,   TRUE,          "both", 5,
  "added_sugar_g",    "negative",   25,    5,   NA,   TRUE,           "both", 10
)

DMI_PHDI_list <- build_DMI(PHDI_wide, ref_PHDI,
                           rho = 0.3, rescale_mode = "none",
                           kappa = 0.2, alpha0 = 3, scale_div = 4,
                           use_max_pts = FALSE)
DMI_PHDI_eq <- DMI_PHDI_list$out
summary(DMI_PHDI_eq$DMI)

DMI_PHDI_list2 <- build_DMI(PHDI_wide, ref_PHDI,
                            rho = 0.3, rescale_mode = "none",
                            kappa = 0.2, alpha0 = 3, scale_div = 4,
                            use_max_pts = TRUE)
DMI_PHDI_wt <- DMI_PHDI_list2$out
summary(DMI_PHDI_wt$DMI)

colnames(DBI_wide)
### 2.9.8 DBI
ref_DBI <- tibble::tribble(
  ~col,               ~direction, ~W0,   ~T0,   ~U0,  ~energy_scaled, ~sex,   ~max_pts,
  "whole_grains_g",   "positive",   5,    225,   NA,   TRUE,           "both", 24,
  "vegetables_g",     "positive",   0,    450,   NA,   TRUE,           "both", 6,
  "fruits_g",         "positive",   0,    300,   NA,   TRUE,           "both", 6,
  "dairy_g",          "positive",   0,    300,   NA,   TRUE,           "both", 6,
  "red_meat_g",       "negative",   150,  90,   NA,   TRUE,           "both", 8,
  "seafood_g",        "positive",   0,    75,   NA,   TRUE,           "both", 4,
  "eggs_g",           "positive",   0,    46,   NA,   TRUE,           "both", 8,
  "sodium_mg",        "negative",   15000, 5000,   NA,   FALSE,          "both", 6
)

DMI_DBI_list <- build_DMI(DBI_wide, ref_DBI,
                          rho = 0.3, rescale_mode = "none",
                          kappa = 0.2, alpha0 = 3, scale_div = 4,
                          use_max_pts = FALSE)
DMI_DBI_eq <- DMI_DBI_list$out
summary(DMI_DBI_eq$DMI)

DMI_DBI_list2 <- build_DMI(DBI_wide, ref_DBI,
                           rho = 0.3, rescale_mode = "none",
                           kappa = 0.2, alpha0 = 3, scale_div = 4,
                           use_max_pts = TRUE)
DMI_DBI_wt <- DMI_DBI_list2$out
summary(DMI_DBI_wt$DMI)

### Continue with the DMI objects created above; no historical RDS cache is required.

merge_keys <- c("iso3", "year", "female", "age_group", "energy_kcal")
DMI_AHEI_combined <- DMI_AHEI_wt %>%
  left_join(DMI_AHEI_eq %>% select(all_of(merge_keys), DMI_2 = DMI),
            by = merge_keys)
DMI_MSDPS_combined <- DMI_MSDPS_wt %>%
  left_join(DMI_MSDPS_eq %>% select(all_of(merge_keys), DMI_2 = DMI),
            by = merge_keys)
DMI_hPDI_combined <- DMI_hPDI_wt %>%
  left_join(DMI_hPDI_eq %>% select(all_of(merge_keys), DMI_2 = DMI),
            by = merge_keys)
DMI_DASH_combined <- DMI_DASH_wt %>%
  left_join(DMI_DASH_eq %>% select(all_of(merge_keys), DMI_2 = DMI),
            by = merge_keys)
DMI_DRRD_combined <- DMI_DRRD_wt %>%
  left_join(DMI_DRRD_eq %>% select(all_of(merge_keys), DMI_2 = DMI),
            by = merge_keys)
DMI_WCRF_combined <- DMI_WCRF_wt %>%
  left_join(DMI_WCRF_eq %>% select(all_of(merge_keys), DMI_2 = DMI),
            by = merge_keys)
DMI_PHDI_combined <- DMI_PHDI_wt %>%
  left_join(DMI_PHDI_eq %>% select(all_of(merge_keys), DMI_2 = DMI),
            by = merge_keys)
DMI_DBI_combined <- DMI_DBI_wt %>%
  left_join(DMI_DBI_eq %>% select(all_of(merge_keys), DMI_2 = DMI),
            by = merge_keys)

keys <- c("iso3","year","female","age_group","energy_kcal")

.scale_WTU <- function(W0, T0, U0, energy_kcal, energy_scaled, base_kcal = 2000){
  r <- if (isTRUE(energy_scaled)) energy_kcal / base_kcal else 1
  W <- W0 * r; T <- T0 * r; U <- ifelse(is.na(U0), NA_real_, U0 * r)
  list(W=W, T=T, U=U)
}
.score_linear <- function(x, W, T, direction = c("positive","negative"), max_pts = 10){
  direction <- match.arg(direction)
  if (any(!is.finite(c(x,W,T)))) return(NA_real_)
  if (direction == "positive"){
    if (x <= W) return(0)
    if (x >= T) return(max_pts)
    return(max_pts * (x - W) / (T - W))
  } else {
    if (x >= W) return(0)
    if (x <= T) return(max_pts)
    return(max_pts * (W - x) / (W - T))
  }
}
.score_bi <- function(x, W, T, direction = c("positive","negative"), max_pts = 10){
  direction <- match.arg(direction)
  if (any(!is.finite(c(x,W,T)))) return(NA_real_)
  gap <- abs(T - W); if (gap == 0) return(NA_real_)
  if (direction == "positive"){
    U <- T + gap
    if (x <= W || x >= U) return(0)
    if (x <= T) return(max_pts * (x - W) / gap)
    return(max_pts * (U - x) / gap)
  } else {
    U <- T - gap
    if (x >= W || x <= U) return(0)
    if (x >= T) return(max_pts * (W - x) / gap)
    return(max_pts * (x - U) / gap)
  }
}
.score_step <- function(x, W, T, direction = c("positive","negative"), max_pts, mid_pts){
  direction <- match.arg(direction)
  if (any(!is.finite(c(x,W,T,max_pts,mid_pts)))) return(NA_real_)
  if (direction == "positive"){
    if (x <  W) return(0)
    if (x >= T) return(max_pts)
    return(mid_pts)
  } else {
    if (x >  W) return(0)
    if (x <= T) return(max_pts)
    return(mid_pts)
  }
}

## 2.10 Calculate the conventional dietary scores ================
score_from_ref <- function(data_wide, ref_table, out_total = "TOTAL_score"){
  ref_use <- ref_table %>% filter(col %in% names(data_wide))
  df <- data_wide
  for (i in seq_len(nrow(ref_use))){
    r <- ref_use[i,]
    cname <- r$col
    sname <- paste0(cname, "_score")

    new_vals <- pmap_dbl(
      list(df[[cname]], df$energy_kcal, df$female),
      function(x, kcal, fem){
        if (!(r$sex %in% c("both","F","M"))) return(NA_real_)
        if (r$sex == "F" && !isTRUE(fem == 1)) return(NA_real_)
        if (r$sex == "M" && !isTRUE(fem == 0)) return(NA_real_)
        th <- .scale_WTU(r$W0, r$T0, r$U0, energy_kcal = kcal, energy_scaled = r$energy_scaled)
        W <- th$W; T <- th$T
        if (r$rule == "linear"){
          .score_linear(x, W, T, direction = r$direction, max_pts = r$max_pts)
        } else if (r$rule == "bi"){
          .score_bi(x, W, T, direction = r$direction, max_pts = r$max_pts)
        } else if (r$rule == "step"){
          .score_step(x, W, T, direction = r$direction, max_pts = r$max_pts, mid_pts = r$mid_pts)
        } else NA_real_
      }
    )

    if (!sname %in% names(df)) {
      df[[sname]] <- new_vals
    } else {
      df[[sname]] <- ifelse(is.na(df[[sname]]), new_vals, df[[sname]])
    }
  }
  df %>% mutate("{out_total}" := rowSums(across(ends_with("_score")), na.rm = TRUE))
}

class(AHEI_wide$female)
ref_AHEI_1 <- tibble::tribble(
  ~col,               ~direction, ~W0,    ~T0,    ~U0,  ~energy_scaled, ~sex,   ~rule,    ~max_pts, ~mid_pts,
  "whole_grains_g",   "positive",   0,      75,    NA,   TRUE,           "F",   "linear",    10,      NA,
  "whole_grains_g",   "positive",   0,      90,    NA,   TRUE,           "M",   "linear",    10,      NA,
  "vegetables_g",     "positive",   0,     600,    NA,   TRUE,         "both", "linear",    10,      NA,
  "fruits_g",         "positive",   0,     472,    NA,   TRUE,         "both", "linear",    10,      NA,
  "nuts_g",           "positive",   0,    28.35,   NA,   TRUE,         "both", "linear",    10,      NA,
  "legumes_g",        "positive",   0,    28.35,   NA,   TRUE,         "both", "linear",    10,      NA,
  "red_meat_g",       "negative", 169.5,    0,     NA,   TRUE,         "both", "linear",    10,      NA,
  "processed_meat_g", "negative",  64.5,    0,     NA,   TRUE,         "both", "linear",    10,      NA,
  "ssb_juice_g",      "negative",  248,     0,     NA,   TRUE,         "both", "linear",    10,      NA,
  "sodium_mg",        "negative", 3337,  1112,     NA,   FALSE,          "F",  "linear",    10,      NA,
  "sodium_mg",        "negative", 5271,  1612,     NA,   FALSE,          "M",  "linear",    10,      NA,
  "lc_n3_pctE",       "positive",   0.0,   0.12,   NA,   FALSE,        "both", "linear",    10,      NA,
  "pufa_pctE",        "positive",   2.0,   10.0,   NA,   FALSE,        "both", "linear",    10,      NA
)

AHEI_scored <- score_from_ref(AHEI_wide, ref_AHEI_1, out_total = "AHEI_score")
summary(AHEI_scored$AHEI_score)

### MSDPS
ref_MSDPS_1 <- tibble::tribble(
  ~col,                ~direction,  ~W0,  ~T0,   ~U0, ~energy_scaled, ~sex,  ~rule, ~max_pts, ~mid_pts,
  "whole_grains_g",    "positive",    0,   120,   NA,  TRUE,          "both","bi", 10,     NA,
  "vegetables_g",      "positive",    0,   720,   NA,  TRUE,          "both","bi", 10,     NA,
  "fruits_g",          "positive",    0,   354,   NA,  TRUE,          "both","bi", 10,     NA,
  "dairy_g",           "positive",    0,   490,   NA,  TRUE,          "both","bi", 10,     NA,
  "seafood_g",         "positive",    0,    73,   NA,  TRUE,          "both","bi", 10,     NA,
  "legumes_g",         "positive",    0,  16.2,   NA,  TRUE,          "both","bi", 10,     NA,
  "nuts_g",            "positive",    0,  16.2,   NA,  TRUE,          "both","bi", 10,     NA,
  "starchy_g",         "positive",    0,  68.6,   NA,  TRUE,          "both","bi", 10,     NA,
  "eggs_g",            "positive",    0,  23.6,   NA,  TRUE,          "both","bi", 10,     NA,
  "meats_g",           "negative", 16.1,     0,   NA,  TRUE,          "both","bi", 10,     NA
)

MSDPS_scored <- score_from_ref(MSDPS_wide, ref_MSDPS_1, out_total = "MSDPS_score")
summary(MSDPS_scored$MSDPS_score)

### WCRF
ref_WCRF_1 <- tibble::tribble(
  ~col,               ~direction,  ~W0,  ~T0,  ~U0, ~energy_scaled, ~sex,  ~rule,   ~max_pts, ~mid_pts,
  "fruits_veg_g",     "positive",   200,  400,  NA,  TRUE,          "both","step", 0.5,       0.25,
  "fiber_g",          "positive",    15,   30,  NA,  TRUE,          "both","step", 0.5,       0.25,
  "ssb_g",            "negative",   250,    0,  NA,  TRUE,          "both","step", 1.0,       0.5,
  "processed_meat_g", "negative",  14.2,    3,  NA,  TRUE,          "both","step", 0.5,       0.25,
  "red_meat_g",       "negative",    71,    35,  NA,  TRUE,          "both","step", 0.5,       0.25
)
WCRF_scored <- score_from_ref(WCRF_wide, ref_WCRF_1, out_total = "WCRF_score")
summary(WCRF_scored$WCRF_score)

### PHDI
ref_PHDI_1 <- tibble::tribble(
  ~col,                   ~direction, ~W0,   ~T0,   ~U0, ~energy_scaled, ~sex,  ~rule,    ~max_pts, ~mid_pts,
  "whole_grains_g",       "positive",   0,     75,   NA,  TRUE,           "F",  "linear",   10,      NA,
  "whole_grains_g",       "positive",   0,     90,   NA,  TRUE,           "M",  "linear",   10,      NA,
  "starchy_g",            "negative", 200,     50,   NA,  TRUE,          "both","linear",   10,      NA,
  "vegetables_g",         "positive",   0,    300,   NA,  TRUE,          "both","linear",   10,      NA,
  "fruits_g",             "positive",   0,    200,   NA,  TRUE,          "both","linear",   10,      NA,
  "dairy_g",              "negative",1000,    250,   NA,  TRUE,          "both","linear",   10,      NA,
  "eggs_g",               "negative", 120,     13,   NA,  TRUE,          "both","linear",   10,      NA,
  "red_meat_g",           "negative", 100,     14,   NA,  TRUE,          "both","linear",   10,      NA,
  "processed_meat_g",     "negative", 100,     14,   NA,  TRUE,          "both","linear",   10,      NA,
  "seafood_g",            "positive",   0,     28,   NA,  TRUE,          "both","linear",   10,      NA,
  "nuts_g",               "positive",   0,     50,   NA,  TRUE,          "both","linear",   10,      NA,
  "legumes_g",            "positive",   0,    100,   NA,  TRUE,          "both","linear",   5,      NA,
  "added_sugar_g",        "negative",  25,      5,   NA,  TRUE,          "both","linear",   10,      NA
)

PHDI_scored <- score_from_ref(PHDI_wide, ref_PHDI_1, out_total = "PHDI_score")
summary(PHDI_scored$PHDI_score)

colnames(hPDI_wide)
colnames(DASH_wide)
colnames(DRRD_wide)

score_quintile <- function(df, value_col, direction = c("positive","negative"),
                           group_vars = c("year","age_group","female"),
                           max_pts = 5L, min_pts = 1L) {
  direction <- match.arg(direction)
  stopifnot(value_col %in% names(df))
  df %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(
      q = ntile(.data[[value_col]], 5L)  # 1..5
    ) %>%
    ungroup() %>%
    mutate(
      !!paste0(value_col,"_score") :=
        if (direction == "positive") {
          scales::rescale(q, to = c(min_pts, max_pts), from = c(1,5))
        } else {
          scales::rescale(6 - q, to = c(min_pts, max_pts), from = c(1,5))
        }
    )
}

ref_hPDI_1 <- tibble::tribble(
  ~col,                 ~direction,
  "whole_grains_g",     "positive",
  "vegetables_g",       "positive",
  "fruits_g",           "positive",
  "nuts_g",             "positive",
  "legumes_g",          "positive",
  "tea_coffee_g",       "positive",
  "fruit_juice_g",      "negative",
  "refined_grains_g",   "negative",
  "potatoes_g",         "negative",
  "ssb_g",              "negative",
  "dairy_g",            "negative",
  "eggs_g",             "negative",
  "seafood_g",          "negative",
  "red_meat_g",         "negative",
  "processed_meat_g",   "negative"
)

score_hPDI <- function(hPDI_wide, ref = ref_hPDI_1,
                       group_vars = c("year","age_group","female")){
  df <- hPDI_wide
  for (i in seq_len(nrow(ref))){
    r <- ref[i,]
    df <- score_quintile(df, r$col, direction = r$direction,
                         group_vars = group_vars, max_pts = 5L, min_pts = 1L)
  }
  df %>%
    mutate(hPDI_score = rowSums(across(ends_with("_score")), na.rm = TRUE))
}
hPDI_scored <- score_hPDI(hPDI_wide, ref_hPDI_1)

### DASH
ref_DASH_1 <- tibble::tribble(
  ~col,                 ~direction,
  "whole_grains_g",     "positive",
  "vegetables_g",       "positive",
  "fruits_g",           "positive",
  "nuts_legumes_g",     "positive",
  "dairy_g",            "positive",
  "sodium_mg",          "negative",
  "ssb_g",              "negative",
  "meats_g",            "negative"
)
score_DASH <- function(DASH_wide, ref = ref_DASH_1,
                       group_vars = c("year","age_group","female")){
  df <- DASH_wide
  for (i in seq_len(nrow(ref))){
    r <- ref[i,]
    df <- score_quintile(df, r$col, direction = r$direction,
                         group_vars = group_vars, max_pts = 5L, min_pts = 1L)
  }
  df %>%
    mutate(DASH_score = rowSums(across(ends_with("_score")), na.rm = TRUE))
}
DASH_scored <- score_DASH(DASH_wide, ref_DASH_1)

### DRRD
ref_DRRD_1 <- tibble::tribble(
  ~col,                 ~direction,
  "ps_ratio",           "positive",
  "coffee_g",           "positive",
  "nuts_g",             "positive",
  "fruits_g",           "positive",
  "ssb_juice_g",        "negative",
  "meats_g",            "negative"
)
score_DRRD <- function(DRRD_wide, ref = ref_DRRD_1,
                       group_vars = c("year","age_group","female")){
  df <- DRRD_wide
  for (i in seq_len(nrow(ref))){
    r <- ref[i,]
    df <- score_quintile(df, r$col, direction = r$direction,
                         group_vars = group_vars, max_pts = 5L, min_pts = 1L)
  }
  df %>%
    mutate(DRRD_score = rowSums(across(ends_with("_score")), na.rm = TRUE))
}
DRRD_scored <- score_DRRD(DRRD_wide, ref_DRRD_1)

score_bins <- function(x, breaks) {
  if (is.na(x)) return(NA_real_)
  for (b in breaks) {
    if (x >= b$lo && x <= b$hi) return(b$pts)
  }
  return(NA_real_)
}
score_dbi_wholegrains <- function(x){
  if (is.na(x)) return(NA_real_)
  if (x < 5)    return(-12)
  if (x <= 225) return(-12 + (x - 5) / (225 - 5) * 12)    # -12 → 0
  if (x <= 275) return(0)
  if (x <= 495) return( (x - 275) / (495 - 275) * 12 )    # 0 → +12
  return(12)
}

score_dbi_veg <- function(x){
  if (is.na(x)) return(NA_real_)
  if (x >= 450) return(0)
  k <- floor((450 - x) / 90) + 1
  pts <- -pmin(k, 6)
  return(pts)
}

score_dbi_fruit <- function(x){
  if (is.na(x)) return(NA_real_)
  if (x >= 300) return(0)
  k <- floor((300 - x) / 60) + 1
  pts <- -pmin(k, 6)
  return(pts)
}

score_dbi_dairy <- function(x){
  if (is.na(x)) return(NA_real_)
  if (x >= 300) return(0)
  k <- floor((300 - x) / 60) + 1
  pts <- -pmin(k, 6)
  return(pts)
}

score_dbi_redmeat <- function(x){
  if (is.na(x)) return(NA_real_)
  brks <- list(
    list(lo=0, hi=0, pts=-4),
    list(lo=  0.001, hi= 20, pts=-3),
    list(lo= 20.001, hi= 40, pts=-2),
    list(lo= 40.001, hi= 60, pts=-1),
    list(lo= 60.001, hi= 90, pts= 0),
    list(lo= 90.001, hi=110, pts= 1),
    list(lo=110.001, hi=130, pts= 2),
    list(lo=130.001, hi=150, pts= 3),
    list(lo=150.001, hi=Inf, pts= 4)
  )
  score_bins(x, brks)
}

score_dbi_seafood <- function(x){
  if (is.na(x)) return(NA_real_)
  brks <- list(
    list(lo=0,  hi= 0,  pts=-4),
    list(lo=0.001,  hi=24,  pts=-3),
    list(lo=24.001, hi=49,  pts=-2),
    list(lo=49.001, hi=74,  pts=-1),
    list(lo=74.001, hi=Inf, pts= 0)
  )
  score_bins(x, brks)
}

score_dbi_eggs <- function(x){
  if (is.na(x)) return(NA_real_)
  brks <- list(
    list(lo=0,  hi= 0,  pts=-4),
    list(lo=0.001,  hi=15,  pts=-3),
    list(lo=15.001, hi=30,  pts=-2),
    list(lo=30.001, hi=45,  pts=-1),
    list(lo=45.001, hi=55,  pts= 0),
    list(lo=55.001, hi=70,  pts= 1),
    list(lo=70.001, hi=85,  pts= 2),
    list(lo=85.001, hi=100, pts= 3),
    list(lo=100.001,hi=Inf, pts= 4)
  )
  score_bins(x, brks)
}

score_dbi_sodium <- function(x_mg){
  if (is.na(x_mg)) return(NA_real_)
  if (x_mg < 5000)  return(0)
  if (x_mg <= 6000) return(1)
  step <- floor((x_mg - 6000) / 1500) + 1
  return( pmin(1 + step, 6) )
}

score_DBI <- function(DBI_wide){
  df <- DBI_wide %>%
    mutate(
      dbi_whole_grains = map_dbl(whole_grains_g, score_dbi_wholegrains),
      dbi_veg          = map_dbl(vegetables_g,   score_dbi_veg),
      dbi_fruits       = map_dbl(fruits_g,       score_dbi_fruit),
      dbi_dairy        = map_dbl(dairy_g,        score_dbi_dairy),
      dbi_red_meat     = map_dbl(red_meat_g,     score_dbi_redmeat),
      dbi_seafood      = map_dbl(seafood_g,      score_dbi_seafood),
      dbi_eggs         = map_dbl(eggs_g,         score_dbi_eggs),
      dbi_sodium       = map_dbl(sodium_mg,      score_dbi_sodium)
    ) %>%
    mutate(
      DBI_score = rowSums(across(starts_with("dbi_")), na.rm = TRUE)
    )
  df
}

DBI_scored <- score_DBI(DBI_wide)

sapply(AHEI_scored , function(x) sum(is.na(x)))
sapply(MSDPS_scored , function(x) sum(is.na(x)))
sapply(WCRF_scored , function(x) sum(is.na(x)))
sapply(PHDI_scored , function(x) sum(is.na(x)))
sapply(hPDI_scored , function(x) sum(is.na(x)))
sapply(DASH_scored  , function(x) sum(is.na(x)))
sapply(DRRD_scored , function(x) sum(is.na(x)))
sapply(DBI_scored , function(x) sum(is.na(x)))

DBI_missing <- DBI_scored[is.na(DBI_scored$dbi_seafood), ]

merge_keys <- c("iso3", "year", "female", "age_group", "energy_kcal")
DMI_AHEI_combined <- DMI_AHEI_combined %>%
  left_join(AHEI_scored %>% select(all_of(merge_keys), AHEI_score),
            by = merge_keys)
DMI_MSDPS_combined <- DMI_MSDPS_combined %>%
  left_join(MSDPS_scored %>% select(all_of(merge_keys), MSDPS_score),
            by = merge_keys)
DMI_hPDI_combined <- DMI_hPDI_combined %>%
  left_join(hPDI_scored %>% select(all_of(merge_keys), hPDI_score),
            by = merge_keys)

DMI_DASH_combined <- DMI_DASH_combined %>%
  left_join(DASH_scored %>% select(all_of(merge_keys), DASH_score),
            by = merge_keys)

DMI_DRRD_combined <- DMI_DRRD_combined %>%
  left_join(DRRD_scored %>% select(all_of(merge_keys), DRRD_score),
            by = merge_keys)

DMI_WCRF_combined <- DMI_WCRF_combined %>%
  left_join(WCRF_scored %>% select(all_of(merge_keys), WCRF_score),
            by = merge_keys)

DMI_PHDI_combined <- DMI_PHDI_combined %>%
  left_join(PHDI_scored %>% select(all_of(merge_keys), PHDI_score),
            by = merge_keys)

DMI_DBI_combined <- DMI_DBI_combined %>%
  left_join(DBI_scored %>% select(all_of(merge_keys), DBI_score),
            by = merge_keys)

all_diet_scores <- DMI_AHEI_combined %>%
  select(all_of(merge_keys),
         AHEI_DMI = DMI, AHEI_DMI_2 = DMI_2, AHEI_score) %>%
  left_join(DMI_MSDPS_combined %>%
              select(all_of(merge_keys),
                     MSDPS_DMI = DMI, MSDPS_DMI_2 = DMI_2, MSDPS_score),
            by = merge_keys) %>%
  left_join(DMI_hPDI_combined %>%
              select(all_of(merge_keys),
                     hPDI_DMI = DMI, hPDI_DMI_2 = DMI_2, hPDI_score),
            by = merge_keys) %>%
  left_join(DMI_DASH_combined %>%
              select(all_of(merge_keys),
                     DASH_DMI = DMI, DASH_DMI_2 = DMI_2, DASH_score),
            by = merge_keys) %>%
  left_join(DMI_DRRD_combined %>%
              select(all_of(merge_keys),
                     DRRD_DMI = DMI, DRRD_DMI_2 = DMI_2, DRRD_score),
            by = merge_keys) %>%
  left_join(DMI_WCRF_combined %>%
              select(all_of(merge_keys),
                     WCRF_DMI = DMI, WCRF_DMI_2 = DMI_2, WCRF_score),
            by = merge_keys) %>%
  left_join(DMI_PHDI_combined %>%
              select(all_of(merge_keys),
                     PHDI_DMI = DMI, PHDI_DMI_2 = DMI_2, PHDI_score),
            by = merge_keys) %>%
  left_join(DMI_DBI_combined %>%
              select(all_of(merge_keys),
                     DBI_DMI = DMI, DBI_DMI_2 = DMI_2, DBI_score),
            by = merge_keys)


dmi_cols <- c("AHEI_DMI", "MSDPS_DMI", "hPDI_DMI", "DASH_DMI",
              "DRRD_DMI", "WCRF_DMI", "PHDI_DMI", "DBI_DMI")
dmi_cols_2 <- c("AHEI_DMI_2", "MSDPS_DMI_2", "hPDI_DMI_2", "DASH_DMI_2",
              "DRRD_DMI_2", "WCRF_DMI_2", "PHDI_DMI_2", "DBI_DMI_2")

score_cols <- c("AHEI_score", "MSDPS_score", "hPDI_score", "DASH_score",
                "DRRD_score", "WCRF_score", "PHDI_score", "DBI_score")



## 2.11 Propagate dietary-intake uncertainty through the DMI ================
.impute_U <- function(W0, T0, direction, kappa = 0.2){
  if (is.na(W0) | is.na(T0)) return(NA_real_)
  d <- tolower(direction)

  if (d == "positive") {
    return(T0 * (1 + kappa))
  } else if (d == "negative") {
    U0 <- T0 * (1 - kappa)
    return(max(U0, 0))
  } else {
    return(NA_real_)
  }
}
.validate_WTU <- function(W, T, U){
  low  <- min(W, U)
  high <- max(W, U)
  is.finite(W) && is.finite(T) && is.finite(U) && (T >= low) && (T <= high)
}

.prepare_ref <- function(ref_table, data_names) {
  ref_use <- ref_table[ref_table$col %in% data_names, , drop = FALSE]
  ref_use$direction <- tolower(ref_use$direction)
  ref_use$sex <- ifelse(is.na(ref_use$sex), "both", ref_use$sex)
  ref_use
}

.DMI_row <- function(row_idx, dt, ref_use, rho, kappa, alpha0, scale_div, use_max_pts, base_kcal) {
  female_i    <- dt$female[row_idx]
  energy_kcal <- dt$energy_kcal[row_idx]

  n <- nrow(ref_use)
  probs <- numeric(n)
  wvec  <- if (isTRUE(use_max_pts)) ref_use$max_pts else rep(1, n)

  for (i in seq_len(n)) {
    r <- ref_use[i,]

    if (r$sex == "F" && female_i != 1) { probs[i] <- NA_real_; next }
    if (r$sex == "M" && female_i != 0) { probs[i] <- NA_real_; next }

    x <- dt[[ r$col ]][row_idx]

    U0i <- if (is.na(r$U0)) {
      .impute_U(r$W0, r$T0, r$direction, kappa)
    } else r$U0

    scale_factor <- if (isTRUE(r$energy_scaled)) energy_kcal / base_kcal else 1
    W <- r$W0 * scale_factor
    T <- r$T0 * scale_factor
    U <- U0i * scale_factor

    low  <- min(W, U)
    high <- max(W, U)
    if (!is.finite(W) || !is.finite(T) || !is.finite(U) ||
        T < low || T > high || is.na(x)) {
      probs[i] <- NA_real_
      next
    }

    xi    <- T
    omega <- max(abs(U - W) / scale_div, 1e-6)
    alpha <- if (r$direction == "positive") -abs(alpha0) else abs(alpha0)

    cdf <- sn::psn(x, xi, omega, alpha)
    p   <- if (r$direction == "positive") cdf else (1 - cdf)
    p   <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    probs[i] <- p
  }

  keep  <- is.finite(probs)
  probs <- probs[keep]
  w     <- wvec[keep]
  if (!length(probs)) return(NA_real_)
  w <- w / sum(w)

  eps   <- 1e-12
  probs <- pmin(pmax(probs, eps), 1 - eps)
  rho0  <- rho

  if (abs(rho0) < 1e-6) {
    return(exp(sum(w * log(probs))))
  } else {
    z <- rho0 * log(probs)
    m <- max(z)
    s <- sum(w * exp(z - m))
    return(exp((log(s) + m) / rho0))
  }
}

calc_DMI_CI_fast <- function(data, ref_table,
                             n_iter = 500,
                             rho = 0.3,
                             kappa = 0.2,
                             alpha0 = 3,
                             scale_div = 4,
                             use_max_pts = TRUE,
                             base_kcal = 2000,
                             rescale_mode = c("pattern","none"),
                             workers = max(1, parallel::detectCores() - 1)) {

  rescale_mode <- match.arg(rescale_mode)
  dt <- as.data.table(data)
  stopifnot(all(c("female","energy_kcal") %in% names(dt)))
  ref_use <- .prepare_ref(ref_table, names(dt))

  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession, workers = workers)

  meas_cols <- grep("(_g|_mg|_pctE)$", names(dt), value = TRUE)
  has_lohi  <- meas_cols[
    paste0(meas_cols, "_lo") %in% names(dt) &
      paste0(meas_cols, "_hi") %in% names(dt)
  ]

  sim_list <- future_lapply(seq_len(n_iter), function(b) {
    dd <- copy(dt)

    for (cn in has_lohi) {
      lo <- dd[[ paste0(cn, "_lo") ]]
      hi <- dd[[ paste0(cn, "_hi") ]]
      dd[[ cn ]] <- runif(nrow(dd), lo, hi)
    }

    DMI_raw <- vapply(seq_len(nrow(dd)), .DMI_row,
                      FUN.VALUE = numeric(1),
                      dt = dd, ref_use = ref_use,
                      rho = rho, kappa = kappa, alpha0 = alpha0,
                      scale_div = scale_div, use_max_pts = use_max_pts, base_kcal = base_kcal)

    ### Rescale the simulated score to the requested range.
    if (rescale_mode == "pattern") {
      rng <- range(DMI_raw, na.rm = TRUE)
      if (is.finite(rng[1]) && diff(rng) > 0) {
        DMI <- (DMI_raw - rng[1]) / (rng[2] - rng[1]) * 100
      } else {
        DMI <- DMI_raw * 100
      }
    } else {
      DMI <- DMI_raw * 100
    }
    DMI
  }, future.seed = TRUE)

  DMI_mat <- do.call(cbind, sim_list)

  DMI_summary <- as.data.frame(dt)[, 0, drop = FALSE]
  DMI_summary$DMI_median <- apply(DMI_mat, 1, median,   na.rm = TRUE)
  DMI_summary$DMI_low    <- apply(DMI_mat, 1, quantile, 0.025, na.rm = TRUE)
  DMI_summary$DMI_high   <- apply(DMI_mat, 1, quantile, 0.975, na.rm = TRUE)
  DMI_summary$DMI_sd     <- apply(DMI_mat, 1, sd,       na.rm = TRUE)

  DMI_summary
}

### 2.11.1 AHEI
DMI_AHEI_UI <- calc_DMI_CI_fast(
  data = AHEI_wide,
  ref_table = ref_AHEI,
  n_iter = 500,
  rho = 0.3,
  kappa = 0.2,
  alpha0 = 3,
  scale_div = 4,
  use_max_pts = TRUE,
  rescale_mode = "none",
  workers = parallel::detectCores() - 1
)
summary(DMI_AHEI_UI$DMI_median)
summary(DMI_AHEI_UI$DMI_low)
summary(DMI_AHEI_UI$DMI_high)

stopifnot(nrow(DMI_AHEI_UI) == nrow(DMI_AHEI_combined))
DMI_AHEI<- dplyr::bind_cols(DMI_AHEI_combined, DMI_AHEI_UI)

### 2.11.2 MSDPS
DMI_MSDPS_UI <- calc_DMI_CI_fast(
  data = MSDPS_wide,
  ref_table = ref_MSDPS,
  n_iter = 500,
  rho = 0.3,
  kappa = 0.2,
  alpha0 = 3,
  scale_div = 4,
  use_max_pts = TRUE,
  rescale_mode = "none",
  workers = parallel::detectCores() - 1
)
summary(DMI_MSDPS_UI$DMI_median)
summary(DMI_MSDPS_UI$DMI_low)
summary(DMI_MSDPS_UI$DMI_high)

stopifnot(nrow(DMI_MSDPS_UI) == nrow(DMI_MSDPS_combined))
DMI_MSDPS<- dplyr::bind_cols(DMI_MSDPS_combined, DMI_MSDPS_UI)

### 2.11.3 hPDI
DMI_hPDI_UI <- calc_DMI_CI_fast(
  data = hPDI_wide,
  ref_table = ref_hPDI,
  n_iter = 500,
  rho = 0.3,
  kappa = 0.2,
  alpha0 = 3,
  scale_div = 4,
  use_max_pts = TRUE,
  rescale_mode = "none",
  workers = parallel::detectCores() - 1
)
summary(DMI_hPDI_UI$DMI_median)
summary(DMI_hPDI_UI$DMI_low)
summary(DMI_hPDI_UI$DMI_high)

stopifnot(nrow(DMI_hPDI_UI) == nrow(DMI_hPDI_combined))
DMI_hPDI <- dplyr::bind_cols(DMI_hPDI_combined, DMI_hPDI_UI)

### 2.11.4 DASH
DMI_DASH_UI <- calc_DMI_CI_fast(
  data = DASH_wide,
  ref_table = ref_DASH,
  n_iter = 500,
  rho = 0.3,
  kappa = 0.2,
  alpha0 = 3,
  scale_div = 4,
  use_max_pts = TRUE,
  rescale_mode = "none",
  workers = parallel::detectCores() - 1
)
summary(DMI_DASH_UI$DMI_median)
summary(DMI_DASH_UI$DMI_low)
summary(DMI_DASH_UI$DMI_high)

stopifnot(nrow(DMI_DASH_UI) == nrow(DMI_DASH_combined))
DMI_DASH <- dplyr::bind_cols(DMI_DASH_combined, DMI_DASH_UI)

### 2.11.5 DRRD
DMI_DRRD_UI <- calc_DMI_CI_fast(
  data = DRRD_wide,
  ref_table = ref_DRRD,
  n_iter = 500,
  rho = 0.3,
  kappa = 0.2,
  alpha0 = 3,
  scale_div = 4,
  use_max_pts = TRUE,
  rescale_mode = "none",
  workers = parallel::detectCores() - 1
)
summary(DMI_DRRD_UI$DMI_median)
summary(DMI_DRRD_UI$DMI_low)
summary(DMI_DRRD_UI$DMI_high)

stopifnot(nrow(DMI_DRRD_UI) == nrow(DMI_DRRD_combined))
DMI_DRRD <- dplyr::bind_cols(DMI_DRRD_combined, DMI_DRRD_UI)

### 2.11.6 WCRF
DMI_WCRF_UI <- calc_DMI_CI_fast(
  data = WCRF_wide,
  ref_table = ref_WCRF,
  n_iter = 500,
  rho = 0.3,
  kappa = 0.2,
  alpha0 = 3,
  scale_div = 4,
  use_max_pts = TRUE,
  rescale_mode = "none",
  workers = parallel::detectCores() - 1
)
summary(DMI_WCRF_UI$DMI_median)
summary(DMI_WCRF_UI$DMI_low)
summary(DMI_WCRF_UI$DMI_high)

stopifnot(nrow(DMI_WCRF_UI) == nrow(DMI_WCRF_combined))
DMI_WCRF <- dplyr::bind_cols(DMI_WCRF_combined, DMI_WCRF_UI)

### 2.11.7 PHDI
DMI_PHDI_UI <- calc_DMI_CI_fast(
  data = PHDI_wide,
  ref_table = ref_PHDI,
  n_iter = 500,
  rho = 0.3,
  kappa = 0.2,
  alpha0 = 3,
  scale_div = 4,
  use_max_pts = TRUE,
  rescale_mode = "none",
  workers = parallel::detectCores() - 1
)
summary(DMI_PHDI_UI$DMI_median)
summary(DMI_PHDI_UI$DMI_low)
summary(DMI_PHDI_UI$DMI_high)

stopifnot(nrow(DMI_PHDI_UI) == nrow(DMI_PHDI_combined))
DMI_PHDI <- dplyr::bind_cols(DMI_PHDI_combined, DMI_PHDI_UI)

### 2.11.8 DBI
DMI_DBI_UI <- calc_DMI_CI_fast(
  data = DBI_wide,
  ref_table = ref_DBI,
  n_iter = 500,
  rho = 0.3,
  kappa = 0.2,
  alpha0 = 3,
  scale_div = 4,
  use_max_pts = TRUE,
  rescale_mode = "none",
  workers = parallel::detectCores() - 1
)
summary(DMI_DBI_UI$DMI_median)
summary(DMI_DBI_UI$DMI_low)
summary(DMI_DBI_UI$DMI_high)

stopifnot(nrow(DMI_DBI_UI) == nrow(DMI_DBI_combined))
DMI_DBI <- dplyr::bind_cols(DMI_DBI_combined, DMI_DBI_UI)

vars <- c(
  "SP.URB.TOTL.IN.ZS",
  "EN.POP.DNST",
  "SH.XPD.CHEX.PC.CD"
)

# 3. Assemble health, demographic and environmental data =========================================
## 3.1 Download World Bank covariates ================
wdi_data <- WDI(country = "all", indicator = vars,
                start = 1990, end = 2018)
selected_years <- c(1990, 1995, 2000, 2005, 2010, 2015, 2018)
wdi_selected <- wdi_data %>%
  filter(year %in% selected_years) %>%
  rename(
    iso3 = iso3c,
    urban_rate = SP.URB.TOTL.IN.ZS,
    pop_density = EN.POP.DNST,
    health_exp_pc = SH.XPD.CHEX.PC.CD
  ) %>%
  select(iso3, year, urban_rate, pop_density, health_exp_pc)

length(unique(wdi_selected$iso3))
sapply(wdi_selected, function(x) sum(is.na(x)))

## 3.2 Read SDI, climate, HALE, LE and population data ================
SDI_data <- read_excel(file.path(GDD_AUX_ROOT, "SDI.xlsx"))

colnames(SDI_data)[1] <- "location"
SDI_data$iso3 <- countrycode(SDI_data$location,
                             origin = "country.name",
                             destination = "iso3c")
colnames(SDI_data)
SDI_long <- SDI_data %>%
  pivot_longer(
    cols = c("1990", "1995", "2000", "2005", "2010", "2015", "2018"),
    names_to = "year",
    values_to = "SDI"
  ) %>%
  mutate(
    year = as.integer(year)
  ) %>%
  filter(!is.na(iso3)) %>%
  select(iso3, location, year, SDI)

length(unique(SDI_data$iso3))

data_qihou <- readRDS(file.path(GDD_AUX_ROOT, "data_qihou.rds"))
colnames(data_qihou)

pm_data <- WDI(
  country = "all",
  indicator = "EN.ATM.PM25.MC.M3",  # Annual mean PM2.5 (micrograms per cubic meter)
  start = 1990, end = 2018
)
pm_selected <- pm_data %>%
  filter(year %in% selected_years) %>%
  rename(
    iso3 = iso3c,
    pm25 = EN.ATM.PM25.MC.M3
  ) %>%
  select(iso3, year, pm25)

HALE_data <- read.csv(file.path(GDD_RAW_ROOT, "GBD2023HALE", "HALE2018.csv"))
LE_data <- read.csv(file.path(GDD_RAW_ROOT, "GBD2023HALE", "LE2018.csv"))

HALE_data$iso3 <- countrycode(HALE_data$location,
                              origin = "country.name",
                              destination = "iso3c")
LE_data$iso3 <- countrycode(LE_data$location,
                              origin = "country.name",
                              destination = "iso3c")

population <- fread(file.path(GDD_AUX_ROOT, "WPP2022_Population1JanuaryByAge5GroupSex_Medium.csv"))
population <- population[Time %in% c(1990, 1995, 2000, 2005, 2010, 2015, 2018)]
population <- population[, c(4, 10, 13, 15, 18, 19, 20)]
colnames(population)
population <- population[ISO3_code != "" & !is.na(ISO3_code)]
colnames(population) <- c("iso3", "country", "year", "age_band", "popmale", "popfemale", "poptotal")
length(unique(population$iso3))
str(population)

hale_raw <- HALE_data %>%
  transmute(
    iso3,
    year,
    sex   = tolower(sex),              # "female"/"male"
    age   = trimws(as.character(age)),
    HALE_val = val,
    HALE_low = lower,
    HALE_high = upper,
    sd = (upper - lower) / (2 * 1.96)
  ) %>%
  mutate(sd = if_else(!is.finite(sd) | sd <= 0, NA_real_, sd))

le_raw <- LE_data %>%
  transmute(
    iso3,
    year,
    sex   = tolower(sex),              # "female"/"male"
    age   = trimws(as.character(age)),
    HALE_val = val,
    HALE_low = lower,
    HALE_high = upper,
    sd = (upper - lower) / (2 * 1.96)
  ) %>%
  mutate(sd = if_else(!is.finite(sd) | sd <= 0, NA_real_, sd))

infant_levels <- c("0-6 days","7-27 days","1-5 months","6-11 months","12-23 months")

map_age_to_band <- function(x){
  x <- trimws(as.character(x))
  n <- length(x)
  out <- rep(NA_character_, n)

  is_infant <- x %in% infant_levels
  out[is_infant] = NA_character_

  is_plus <- grepl("^\\d+\\+\\s*years$", x)
  if (any(is_plus)) {
    lo <- as.numeric(sub("\\+.*","", x[is_plus]))
    out[is_plus] <- ifelse(lo >= 65, "65+", NA_character_)
  }

  m <- stringr::str_match(x, "^(\\d+)[^0-9]+(\\d+)\\s*years$")
  ok <- !is.na(m[,1])
  if (any(ok)) {
    lo <- as.numeric(m[ok,2]); hi <- as.numeric(m[ok,3])
    out[ok] <- dplyr::case_when(
      hi <  2              ~ NA_character_,
      hi <= 19             ~ "2-19",
      lo >= 65             ~ "65+",
      lo >= 20 & hi <= 64  ~ "20-64",
      TRUE                 ~ NA_character_
    )
  }
  out
}

hale_clean <- hale_raw %>%
  mutate(age_band = map_age_to_band(age)) %>%
  filter(!is.na(age_band)) %>%
  select(iso3, year, sex, age, age_band, HALE_val, sd, HALE_low, HALE_high)

le_clean <- le_raw %>%
  mutate(age_band = map_age_to_band(age)) %>%
  filter(!is.na(age_band)) %>%
  select(iso3, year, sex, age, age_band, HALE_val, sd, HALE_low, HALE_high)

population <- population %>%
  mutate(
    iso3  = toupper(trimws(iso3)),
    year  = as.integer(year),
    age_wpp = trimws(as.character(age_band)),
    popmale   = as.numeric(popmale),
    popfemale = as.numeric(popfemale)
  )

### Map abbreviated age labels to the labels used in the health datasets.
wpp_expand <- function(df){
  base <- df %>%
    filter(!age_wpp %in% c("0-4","95-99","100+")) %>%
    mutate(age = paste0(age_wpp, " years"))
  ### Approximate ages 2-4 as three fifths of the population aged 0-4.
  part_24 <- df %>%
    filter(age_wpp == "0-4") %>%
    mutate(age = "2-4 years",
           popmale   = popmale   * 3/5,
           popfemale = popfemale * 3/5)
  part_95p <- df %>%
    filter(age_wpp %in% c("95-99","100+")) %>%
    group_by(iso3, year) %>%
    summarise(popmale = sum(popmale, na.rm = TRUE),
              popfemale = sum(popfemale, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(age = "95+ years")
  bind_rows(base, part_24, part_95p) %>%
    select(iso3, year, age, popmale, popfemale)
}

pop_hale_like <- wpp_expand(population)

pop_long <- pop_hale_like %>%
  pivot_longer(cols = c(popmale, popfemale), names_to = "sex", values_to = "pop_w") %>%
  mutate(sex = ifelse(sex == "popmale", "male", "female"))


hale_w <- hale_clean %>%
  left_join(pop_long, by = c("iso3","year","sex","age"))

le_w <- le_clean %>%
  left_join(pop_long, by = c("iso3", "year", "sex", "age"))

hale_band_ui_weighted <- function(df, n_draws = 1000, workers = max(1, parallel::detectCores()-1)) {
  oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
  future::plan(multisession, workers = workers)

  grp_keys <- c("iso3","year","sex","age_band")

  df2 <- df %>%
    filter(is.finite(HALE_val), is.finite(pop_w), pop_w > 0)

  sims <- future_lapply(seq_len(n_draws), function(b){
    draw <- with(df2, rnorm(nrow(df2), mean = HALE_val,
                            sd = ifelse(is.na(sd) | sd<=0, 0, sd)))
    draw[draw < 0] <- 0
    df2 %>%
      mutate(draw = draw) %>%
      group_by(across(all_of(grp_keys))) %>%
      summarise(HALE_draw = stats::weighted.mean(draw, w = pop_w, na.rm = TRUE),
                .groups = "drop")
  }, future.seed = TRUE)

  sim_all <- bind_rows(Map(function(x,i) mutate(x, .iter = i), sims, seq_along(sims)))

  sim_all %>%
    group_by(across(all_of(grp_keys))) %>%
    summarise(
      HALE_median = median(HALE_draw, na.rm = TRUE),
      HALE_low    = quantile(HALE_draw, 0.025, na.rm = TRUE),
      HALE_high   = quantile(HALE_draw, 0.975, na.rm = TRUE),
      HALE_sd     = sd(HALE_draw, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(female = ifelse(sex == "female", 1L, 0L),
           age_group = age_band) %>%
    select(iso3, year, female, age_group, HALE_median, HALE_low, HALE_high, HALE_sd)
}

HALE_band_UI_2 <- hale_band_ui_weighted(hale_w, n_draws = 300)
LE_band_UI <- hale_band_ui_weighted(le_w, n_draws = 300)

sapply(HALE_band_UI_2, function(x) sum(is.na(x)))
sapply(LE_band_UI, function(x) sum(is.na(x)))

colnames(LE_band_UI)[5:8] <- c("LE_median", "LE_low", "LE_high", "LE_sd")
HEALTH_RDS_DIR <- here::here("outputs", "health_data")
dir.create(HEALTH_RDS_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(HALE_band_UI_2, file.path(HEALTH_RDS_DIR, "HALE_band_UI_2.rds"))
saveRDS(LE_band_UI, file.path(HEALTH_RDS_DIR, "LE_band_UI.rds"))

HALE_band_UI <- HALE_band_UI_2


merge_keys <- c("iso3", "year", "female", "age_group", "energy_kcal")
all_diet_scores <- DMI_AHEI %>%
  select(all_of(merge_keys),
         AHEI_DMI = DMI_median, AHEI_DMI_low = DMI_low, AHEI_DMI_high = DMI_high, AHEI_score) %>%
  left_join(DMI_MSDPS %>%
              select(all_of(merge_keys),
                     MSDPS_DMI = DMI_median, MSDPS_DMI_low = DMI_low, MSDPS_DMI_high = DMI_high, MSDPS_score),
            by = merge_keys) %>%
  left_join(DMI_hPDI %>%
              select(all_of(merge_keys),
                     hPDI_DMI = DMI_median, hPDI_DMI_low = DMI_low, hPDI_DMI_high = DMI_high, hPDI_score),
            by = merge_keys) %>%
  left_join(DMI_DASH %>%
              select(all_of(merge_keys),
                     DASH_DMI = DMI_median, DASH_DMI_low = DMI_low, DASH_DMI_high = DMI_high, DASH_score),
            by = merge_keys) %>%
  left_join(DMI_DRRD %>%
              select(all_of(merge_keys),
                     DRRD_DMI = DMI_median, DRRD_DMI_low = DMI_low, DRRD_DMI_high = DMI_high, DRRD_score),
            by = merge_keys) %>%
  left_join(DMI_WCRF %>%
              select(all_of(merge_keys),
                     WCRF_DMI = DMI_median, WCRF_DMI_low = DMI_low, WCRF_DMI_high = DMI_high, WCRF_score),
            by = merge_keys) %>%
  left_join(DMI_PHDI %>%
              select(all_of(merge_keys),
                     PHDI_DMI = DMI_median, PHDI_DMI_low = DMI_low, PHDI_DMI_high = DMI_high, PHDI_score),
            by = merge_keys) %>%
  left_join(DMI_DBI %>%
              select(all_of(merge_keys),
                     DBI_DMI = DMI_median, DBI_DMI_low = DMI_low, DBI_DMI_high = DMI_high, DBI_score),
            by = merge_keys)

le_df <- LE_band_UI  %>%
  rename(LE = LE_median) %>%
  mutate(
    year = as.numeric(year),
    female = ifelse(tolower(female) %in% c("female","f", "women","1"), 1L,
                    ifelse(tolower(female) %in% c("male","m","men","0"), 0L, NA_integer_)),
    female = as.factor(female),
    age_group = as.factor(age_group)
  ) %>%
  select(iso3, year, female, age_group, LE)
hale_df <- HALE_band_UI  %>%
  rename(HALE = HALE_median) %>%
  mutate(
    year = as.numeric(year),
    female = ifelse(tolower(female) %in% c("female","f", "women","1"), 1L,
                    ifelse(tolower(female) %in% c("male","m","men","0"), 0L, NA_integer_)),
    female = as.factor(female),
    age_group = as.factor(age_group)
  ) %>%
  select(iso3, year, female, age_group, HALE)
sdi_df <- SDI_long %>%
  mutate(
    year = as.numeric(year)
  ) %>%
  select(iso3, year, SDI)
clim_df <- data_qihou %>%
  mutate(
    year = as.numeric(year)
  ) %>%
  select(iso3, year, tmin_mean, tmax_mean, prec_mean, elev)
wdi_selected <- wdi_selected %>%
  mutate(year = as.numeric(year))
pm_selected <- pm_selected %>%
  mutate(year = as.numeric(year))

all_diet_scores <- all_diet_scores %>%
  mutate(
    year = as.numeric(year),
    female = as.factor(female),
    age_group = as.factor(age_group)
  )

## 3.3 Merge all analysis variables ================
panel0 <- all_diet_scores %>%
  left_join(wdi_selected %>% select(iso3, year, urban_rate, pop_density, health_exp_pc), by = c("iso3","year")) %>%
  left_join(pm_selected  %>% select(iso3, year, pm25), by = c("iso3","year")) %>%
  left_join(clim_df,      by = c("iso3","year")) %>%
  left_join(sdi_df,       by = c("iso3","year")) %>%
  left_join(hale_df,      by = c("iso3","year","female","age_group")) %>%
  left_join(le_df,      by = c("iso3","year","female","age_group"))

length(unique(panel0$iso3))

sapply(panel0, function(x) sum(is.na(x)))
sum(duplicated(panel0))

# 4. Perform multiple imputation =================================================================
## 4.1 Prepare variables and verify the prespecified missingness threshold ================
MI_M <- 20L
MI_SEED <- 123L
MI_MAXIT <- 20L
MI_THRESHOLD <- 0.05
MI_DIR <- here("multiple_imputation")
dir.create(MI_DIR, recursive = TRUE, showWarnings = FALSE)

covariate_names <- c(
  "urban_rate", "pop_density", "health_exp_pc", "pm25",
  "tmin_mean", "tmax_mean", "prec_mean", "elev", "SDI"
)

covariate_missingness_data <- panel0 %>%
  select(iso3, year, all_of(covariate_names)) %>%
  distinct()

covariate_missingness <- tibble(
  variable = covariate_names,
  n_total = nrow(covariate_missingness_data),
  n_missing = vapply(covariate_missingness_data[covariate_names], function(x) sum(is.na(x)), integer(1)),
  proportion_missing = vapply(
    covariate_missingness_data[covariate_names],
    function(x) mean(is.na(x)),
    numeric(1)
  )
) %>%
  mutate(below_five_percent = proportion_missing < MI_THRESHOLD)

print(covariate_missingness, n = Inf)
readr::write_csv(covariate_missingness, file.path(MI_DIR, "covariate_missingness.csv"))

if (any(!covariate_missingness$below_five_percent)) {
  warning(
    "The <5% missingness criterion is not met for: ",
    paste(covariate_missingness$variable[!covariate_missingness$below_five_percent], collapse = ", "),
    ". Review the manuscript statement before the final analysis."
  )
}

## 4.2 Attach the observed HALE uncertainty limits ================
bands <- HALE_band_UI %>%
  transmute(
    iso3,
    year = as.numeric(year),
    female = case_when(
      tolower(female) %in% c("female","f","women","1") ~ 1L,
      tolower(female) %in% c("male","m","men","0")     ~ 0L,
      TRUE ~ NA_integer_
    ),
    female = as.factor(female),
    age_group = as.factor(age_group),
    HALE_low,
    HALE_high
  )

panel0_mi <- panel0 %>%
  left_join(bands, by = c("iso3", "year", "female", "age_group")) %>%
  mutate(
    HALE = as.numeric(HALE),
    LE = as.numeric(LE),
    HALE_low = as.numeric(HALE_low),
    HALE_high = as.numeric(HALE_high)
  )

## 4.3 Generate and save 20 completed analysis datasets ================
### Impute at the country-year level because all covariates are shared by sex and age strata.
covariate_panel <- panel0_mi %>%
  select(iso3, year, all_of(covariate_names)) %>%
  distinct()
stopifnot(!anyDuplicated(covariate_panel[c("iso3", "year")]))

imputation_data <- covariate_panel %>%
  select(year, all_of(covariate_names))
imputation_method <- mice::make.method(imputation_data)
imputation_method[] <- ""
imputed_covariates <- covariate_names[colSums(is.na(imputation_data[covariate_names])) > 0]
imputation_method[imputed_covariates] <- "pmm"

### HALE and LE are outcomes, not imputation targets. They remain unchanged in every dataset.
predictor_matrix <- mice::make.predictorMatrix(imputation_data)
predictor_matrix[, ] <- 0
predictor_matrix[imputed_covariates, setdiff(names(imputation_data), imputed_covariates)] <- 1
predictor_matrix[imputed_covariates, imputed_covariates] <- 1
diag(predictor_matrix) <- 0

set.seed(MI_SEED)
mi_object <- mice::mice(
  imputation_data,
  m = MI_M,
  maxit = MI_MAXIT,
  method = imputation_method,
  predictorMatrix = predictor_matrix,
  seed = MI_SEED,
  printFlag = TRUE
)
saveRDS(mi_object, file.path(MI_DIR, "mice_m20.rds"), compress = "xz")

panel0_imputed_list <- lapply(seq_len(MI_M), function(imputation_id) {
  completed_covariates <- mice::complete(mi_object, action = imputation_id) %>%
    select(all_of(covariate_names)) %>%
    bind_cols(covariate_panel %>% select(iso3, year), .) %>%
    distinct(iso3, year, .keep_all = TRUE)

  completed_data <- panel0_mi %>%
    select(-all_of(covariate_names)) %>%
    left_join(completed_covariates, by = c("iso3", "year")) %>%
    mutate(.imp = imputation_id) %>%
    distinct(iso3, year, age_group, female, .keep_all = TRUE)

  output_file <- file.path(MI_DIR, sprintf("panel0_imputed_%02d.rds", imputation_id))
  saveRDS(completed_data, output_file, compress = "xz")
  completed_data
})

panel0_imputed <- panel0_imputed_list[[1]]
print(sapply(panel0_imputed, function(x) sum(is.na(x))))
str(panel0_imputed)

# 5. Construct the HALE/LE sensitivity outcome ===================================================
## 5.1 Read all completed datasets ================
imputed_files <- file.path(MI_DIR, sprintf("panel0_imputed_%02d.rds", seq_len(MI_M)))
stopifnot(all(file.exists(imputed_files)))
panel0_imputed_list <- lapply(imputed_files, readRDS)

LE_for_join <- LE_band_UI %>%
  transmute(
    iso3,
    year = as.numeric(year),
    female = factor(as.character(female), levels = c("0","1")),
    age_group = factor(age_group, levels = levels(panel0_imputed$age_group)),
    LE_low,
    LE_high,
    LE_sd
  )

panel0_combined_list <- lapply(panel0_imputed_list, function(completed_data) {
  completed_data %>%
    left_join(LE_for_join, by = c("iso3", "year", "female", "age_group")) %>%
    mutate(HALE_LE = if_else(is.finite(LE) & LE > 0, HALE / LE, NA_real_))
})

### Use a common country set across all imputations so that Rubin pooling compares like with like.
iso3_drop <- panel0_combined_list %>%
  map(~ filter(.x, HALE_LE > 1) %>% distinct(iso3) %>% pull(iso3)) %>%
  unlist(use.names = FALSE) %>%
  unique()

set.seed(42)
RATIO_BOOT_B <- as.integer(Sys.getenv("GDD_RATIO_BOOT_B", "800"))
n_boot <- RATIO_BOOT_B
dist   <- "normal"
tol    <- 1.01

draw_from_ci <- function(mu, lo, hi, n, dist = c("normal","uniform")) {
  dist <- match.arg(dist)
  if (!is.finite(lo) || !is.finite(hi) || hi <= 0 || lo > hi) return(rep(NA_real_, n))
  if (dist == "uniform") {
    x <- runif(n, min = lo, max = hi)
  } else {
    sd <- (hi - lo) / (2 * 1.96)
    if (!is.finite(mu)) mu <- (lo + hi)/2
    x <- rnorm(n, mean = mu, sd = sd)
    x[x < lo] <- lo; x[x > hi] <- hi
  }
  x
}

boot_ratio_row <- function(h_mu, h_lo, h_hi, l_mu, l_lo, l_hi,
                           n = 1000, dist = "normal", tol = 1.01) {
  if (!all(is.finite(c(h_lo, h_hi, l_lo, l_hi))) || h_hi <= h_lo || l_hi <= l_lo) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }
  h <- draw_from_ci(h_mu, h_lo, h_hi, n, dist)
  l <- draw_from_ci(l_mu, l_lo, l_hi, n, dist)
  h[h < 0]  <- 0
  l[l <= 0] <- 1e-6
  r <- h / l
  r <- ifelse(r <= tol, pmin(r, 1.0), NA_real_)
  c(
    stats::median(r, na.rm = TRUE),
    stats::quantile(r, 0.025, na.rm = TRUE),               # 2.5%
    stats::quantile(r, 0.975, na.rm = TRUE),               # 97.5%
    stats::sd(r, na.rm = TRUE)                             # SD
  )
}

## 5.2 Construct and save the sensitivity outcome in each completed dataset ================
panel_boot_list <- map2(panel0_combined_list, seq_len(MI_M), function(completed_data, imputation_id) {
  panel_boot_i <- completed_data %>%
    filter(!iso3 %in% iso3_drop) %>%
    mutate(HALE_LE_raw = if_else(is.finite(LE) & LE > 0, pmin(HALE / LE, 1.0), NA_real_)) %>%
    rowwise() %>%
    mutate(
      .res = list(
        boot_ratio_row(
          h_mu = HALE, h_lo = HALE_low, h_hi = HALE_high,
          l_mu = LE, l_lo = LE_low, l_hi = LE_high,
          n = n_boot, dist = dist, tol = tol
        )
      ),
      HALE_LE = .res[[1]],
      HALE_LE_low = .res[[2]],
      HALE_LE_high = .res[[3]],
      HALE_LE_sd = .res[[4]]
    ) %>%
    ungroup() %>%
    select(-.res)

  saveRDS(
    panel_boot_i,
    file.path(MI_DIR, sprintf("panel_boot_imputed_%02d.rds", imputation_id)),
    compress = "xz"
  )
  panel_boot_i
})

panel_boot <- panel_boot_list[[1]]
saveRDS(panel_boot, "panel_boot.rds", compress = "xz")
print(summary(panel_boot$HALE_LE))
cat("Common country count across imputations:", n_distinct(panel_boot$iso3), "\n")

# 6. Perform spatial autocorrelation analyses =====================================================
## 6.1 Define shared spatial settings and helper functions ================
colnames(panel0_imputed)

SPATIAL_YEARS <- c(1990, 1995, 2000, 2005, 2010, 2015, 2018)
K_KNN   <- 6
B_BOOT  <- 800
ALPHA_P <- 0.05
SEED    <- 20251022; set.seed(SEED)

DMI_KEYS <- c("AHEI","MSDPS","hPDI","DASH","DRRD","WCRF","PHDI","DBI")
normalize_female <- function(x){
  if (is.numeric(x) || is.integer(x)) return(as.integer(x))
  y <- tolower(as.character(x))
  dplyr::case_when(
    y %in% c("1","female","f","woman","women") ~ 1L,
    y %in% c("0","male","m","man","men")       ~ 0L,
    TRUE ~ NA_integer_
  )
}

build_world_pts <- function(){
  world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
    sf::st_make_valid() %>%
    dplyr::mutate(iso3 = countrycode::countrycode(iso_a3, "iso3c", "iso3c")) %>%
    dplyr::filter(!is.na(iso3), iso_a3 != "-99") %>%
    dplyr::select(iso3, geometry)
  sf::st_point_on_surface(world_sf) %>% dplyr::select(iso3, geometry)
}

build_listw_knn <- function(iso_vec, world_pts, k = K_KNN, k_max = 12){
  pts <- dplyr::filter(world_pts, iso3 %in% iso_vec)
  if (nrow(pts) < 3) return(NULL)
  k <- min(k, nrow(pts)-1); if (k < 1) return(NULL)

  coords <- sf::st_coordinates(pts)
  knn <- spdep::knearneigh(coords, k = k, longlat = TRUE)
  nb  <- spdep::knn2nb(knn)

  while (spdep::n.comp.nb(nb)$nc > 1 && k < min(k_max, nrow(pts)-1)) {
    k <- k + 1
    knn <- spdep::knearneigh(coords, k = k, longlat = TRUE)
    nb  <- spdep::knn2nb(knn)
  }

  nb <- spdep::make.sym.nb(nb)
  spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
}

draw_boot <- function(mu, lo, hi, n = B_BOOT){
  if (is.na(lo) || is.na(hi) || lo >= hi) return(rep(mu, n))
  runif(n, lo, hi)
}

moran_bootstrap_one <- function(df_stratum, value_col, low_col, high_col, world_pts){
  if (nrow(df_stratum) < 3) return(NULL)
  lw <- build_listw_knn(df_stratum$iso3, world_pts, k = K_KNN)
  if (is.null(lw)) return(NULL)

  pts <- world_pts %>% filter(iso3 %in% df_stratum$iso3) %>% st_drop_geometry()
  X   <- pts %>% left_join(df_stratum, by = "iso3")

  nC <- nrow(X)
  Xmat <- map2_dfc(X[[value_col]], seq_len(nC), function(mu, i){
    lo <- X[[low_col]][i]; hi <- X[[high_col]][i]
    tibble::tibble(!!paste0("V", i) := draw_boot(mu, lo, hi, n = B_BOOT))
  })

  I_vals <- numeric(B_BOOT); p_vals <- numeric(B_BOOT)
  for (b in 1:B_BOOT) {
    xb <- as.numeric(Xmat[b, ])
    gt <- tryCatch(moran.test(xb, lw, zero.policy = TRUE), error = function(e) NULL)
    if (is.null(gt)) { I_vals[b] <- NA_real_; p_vals[b] <- NA_real_; next }
    I_vals[b] <- unname(gt$estimate[["Moran I statistic"]])
    p_vals[b] <- gt$p.value
  }

  tibble(
    I_med  = median(I_vals, na.rm = TRUE),
    I_lo   = quantile(I_vals, 0.025, na.rm = TRUE),
    I_hi   = quantile(I_vals, 0.975, na.rm = TRUE),
    p_med  = median(p_vals, na.rm = TRUE),
    p_lo   = quantile(p_vals, 0.025, na.rm = TRUE),
    p_hi   = quantile(p_vals, 0.975, na.rm = TRUE),
    sig_rate  = mean(p_vals < ALPHA_P, na.rm = TRUE),
    n_country = nrow(X)
  )
}

run_global_moran_for_dmi <- function(panel, dmi_key, world_pts){
  v  <- paste0(dmi_key, "_DMI")
  lo <- paste0(dmi_key, "_DMI_low")
  hi <- paste0(dmi_key, "_DMI_high")
  stopifnot(all(c(v,lo,hi) %in% names(panel)))

  strata <- tidyr::expand_grid(female = c(0L,1L),
                               age_band = factor(c("2-19","20-64","65+"),
                                                 levels = c("2-19","20-64","65+")))

  purrr::pmap_dfr(list(strata$female, strata$age_band), function(ff, ab){
    df <- panel %>%
      filter(female == ff, age_band == ab) %>%
      select(iso3, !!v, !!lo, !!hi) %>%
      group_by(iso3) %>%
      summarise(
        !!v  := mean(.data[[v]],  na.rm = TRUE),
        !!lo := mean(.data[[lo]], na.rm = TRUE),
        !!hi := mean(.data[[hi]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(is.finite(.data[[v]])) %>%
      rename(value = !!v, low = !!lo, high = !!hi)

    out <- moran_bootstrap_one(df, "value", "low", "high", world_pts)
    if (is.null(out)) return(tibble(DMI = dmi_key, female = ff, age_band = ab,
                                    I_med = NA_real_, I_lo = NA_real_, I_hi = NA_real_,
                                    p_med = NA_real_, p_lo = NA_real_, p_hi = NA_real_,
                                    sig_rate = NA_real_, n_country = nrow(df)))
    out %>% mutate(DMI = dmi_key, female = ff, age_band = ab)
  })
}

## 6.2 Run and save global Moran's I separately for each study year ================
stopifnot(all(c("iso3", "year", "female", "age_group") %in% names(panel0_imputed)))
world_pts <- build_world_pts()

global_moran_by_year <- setNames(lapply(SPATIAL_YEARS, function(year_i) {
  panel_base_i <- panel0_imputed %>%
    filter(year == year_i) %>%
    mutate(
      age_band = factor(as.character(age_group), levels = c("2-19", "20-64", "65+")),
      female = normalize_female(female)
    ) %>%
    filter(!is.na(age_band), !is.na(female))

  result_i <- lapply(DMI_KEYS, function(pattern) {
    run_global_moran_for_dmi(panel_base_i, pattern, world_pts)
  }) %>%
    bind_rows() %>%
    mutate(year = year_i, female_label = ifelse(female == 1L, "Female", "Male")) %>%
    select(year, DMI, female_label, age_band, n_country, I_med, I_lo, I_hi, p_med, p_lo, p_hi, sig_rate) %>%
    arrange(DMI, female_label, age_band)

  saveRDS(result_i, sprintf("global_moran_%d.rds", year_i), compress = "xz")
  result_i
}), as.character(SPATIAL_YEARS))

res_all <- bind_rows(global_moran_by_year)
print(res_all, n = nrow(res_all))
saveRDS(res_all, "global_moran_all_years.rds", compress = "xz")

## 6.3 Run local Moran's I separately for each study year ================
K_KNN  <- 6
B_BOOT <- 800
ALPHA  <- 0.05
SEED   <- 20251022; set.seed(SEED)

DMI_KEYS <- c("AHEI","MSDPS","hPDI","DASH","DRRD","WCRF","PHDI","DBI")

normalize_female <- function(x){
  if (is.numeric(x) || is.integer(x)) return(as.integer(x))
  y <- tolower(as.character(x))
  dplyr::case_when(
    y %in% c("1","female","f","woman","women") ~ 1L,
    y %in% c("0","male","m","man","men")       ~ 0L,
    TRUE ~ NA_integer_
  )
}
build_world_pts <- function(){
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_make_valid() |>
    dplyr::mutate(iso3 = countrycode::countrycode(iso_a3, "iso3c", "iso3c")) |>
    dplyr::filter(!is.na(iso3), iso_a3 != "-99") |>
    dplyr::select(iso3, geometry) |>
    sf::st_point_on_surface() |>
    dplyr::select(iso3, geometry)
}

build_listw_knn <- function(iso_vec, world_pts, k = K_KNN, k_max = 12){
  pts <- dplyr::filter(world_pts, iso3 %in% iso_vec)
  if (nrow(pts) < 3) return(NULL)
  k <- min(k, nrow(pts)-1); if (k < 1) return(NULL)

  coords <- sf::st_coordinates(pts)
  knn <- spdep::knearneigh(coords, k = k, longlat = TRUE)
  nb  <- spdep::knn2nb(knn)

  while (spdep::n.comp.nb(nb)$nc > 1 && k < min(k_max, nrow(pts)-1)) {
    k <- k + 1
    knn <- spdep::knearneigh(coords, k = k, longlat = TRUE)
    nb  <- spdep::knn2nb(knn)
  }

  nb <- spdep::make.sym.nb(nb)
  spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
}

draw_boot <- function(mu, lo, hi, n = B_BOOT){
  if (is.na(lo) || is.na(hi) || lo >= hi) return(rep(mu, n))
  runif(n, lo, hi)
}

.local_moran_bootstrap_core <- function(df_stratum, world_pts){
  lw <- build_listw_knn(df_stratum$iso3, world_pts)
  if (is.null(lw)) return(NULL)
  pts <- world_pts |> filter(iso3 %in% df_stratum$iso3) |> st_drop_geometry()
  X   <- pts |> left_join(df_stratum, by = "iso3")

  nC <- nrow(X)
  Ii_mat <- matrix(NA_real_, nrow = B_BOOT, ncol = nC)
  p_mat  <- matrix(NA_real_, nrow = B_BOOT, ncol = nC)

  Xmat <- map2_dfc(X$value, seq_len(nC), \(mu, i)
                   tibble::tibble(!!paste0("V", i) := draw_boot(mu, X$low[i], X$high[i], B_BOOT))
  )

  for (b in 1:B_BOOT) {
    xb <- as.numeric(Xmat[b, ])
    loc <- tryCatch(localmoran(xb, lw, zero.policy = TRUE), error = \(e) NULL)
    if (is.null(loc)) next
    Ii_mat[b, ] <- loc[,1];  p_mat[b, ] <- loc[,5]
  }

  tibble(
    iso3 = X$iso3,
    Ii_med = apply(Ii_mat, 2, median,   na.rm = TRUE),
    Ii_lo  = apply(Ii_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    Ii_hi  = apply(Ii_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    p_med  = apply(p_mat,  2, median,   na.rm = TRUE),
    p_lo   = apply(p_mat,  2, quantile, probs = 0.025, na.rm = TRUE),
    p_hi   = apply(p_mat,  2, quantile, probs = 0.975, na.rm = TRUE),
    sig_rate  = apply(p_mat, 2, \(v) mean(v < ALPHA, na.rm = TRUE)),
    n_boot_ok = apply(p_mat, 2, \(v) sum(is.finite(v)))
  )
}

run_local_moran_for_dmi <- function(panel_base, dmi_key, world_pts){
  v  <- paste0(dmi_key, "_DMI")
  lo <- paste0(dmi_key, "_DMI_low")
  hi <- paste0(dmi_key, "_DMI_high")
  stopifnot(all(c(v,lo,hi) %in% names(panel_base)))

  strata <- tidyr::expand_grid(
    female = c(0L,1L),
    age_band = factor(c("2-19","20-64","65+"), levels = c("2-19","20-64","65+"))
  )

  purrr::pmap_dfr(list(strata$female, strata$age_band), \(ff, ab){
    df <- panel_base |>
      filter(female == ff, age_band == ab) |>
      select(iso3, !!v, !!lo, !!hi) |>
      group_by(iso3) |>
      summarise(
        value = mean(.data[[v]],  na.rm = TRUE),
        low   = mean(.data[[lo]], na.rm = TRUE),
        high  = mean(.data[[hi]], na.rm = TRUE),
        .groups = "drop"
      ) |>
      filter(is.finite(value))

    out <- .local_moran_bootstrap_core(df, world_pts)
    if (is.null(out)) return(tibble())
    out |> mutate(DMI = dmi_key, female = ff, age_band = ab)
  })
}

stopifnot(all(c("iso3", "year", "female", "age_group") %in% names(panel0_imputed)))
world_pts <- build_world_pts()

local_moran_by_year <- setNames(lapply(SPATIAL_YEARS, function(year_i) {
  panel_base_i <- panel0_imputed %>%
    filter(year == year_i) %>%
    mutate(
      age_band = factor(as.character(age_group), levels = c("2-19", "20-64", "65+")),
      female = normalize_female(female)
    ) %>%
    filter(!is.na(age_band), !is.na(female))

  result_i <- setNames(
    lapply(DMI_KEYS, function(pattern) {
      run_local_moran_for_dmi(panel_base_i, pattern, world_pts) %>% mutate(year = year_i)
    }),
    DMI_KEYS
  )
  saveRDS(result_i, sprintf("local_moran_bootstrap_%d.rds", year_i), compress = "xz")
  result_i
}), as.character(SPATIAL_YEARS))

saveRDS(local_moran_by_year, "local_moran_bootstrap_all_years.rds", compress = "xz")

# 7. Define the geographically and temporally weighted regression ================================
## 7.1 Load the first completed dataset for function checks ================
panel_boot <- readRDS(panel_boot_files[[1]])

YEARS_AVAILABLE <- c(1990, 1995, 2000, 2005, 2010, 2015, 2018)
LAMBDA <- 0.5
KSI    <- 0.5
KERNEL <- "bisquare"
ADAPT  <- TRUE
SEED   <- 20251022; set.seed(SEED)

DMI_KEYS <- c("AHEI","MSDPS","hPDI","DASH","DRRD","WCRF","PHDI","DBI")

safe_round <- function(x, digits = 3) {
  if (is.numeric(x) && length(x) == 1 && is.finite(x)) round(x, digits) else NA_real_
}

check_and_trim_vars <- function(df, controls, min_n = 30, var_eps = 1e-12) {
  if (nrow(df) < min_n) return(list(ok = FALSE, msg = paste0("too few obs (", nrow(df), ")"),
                                    kept = controls, drop = character(0)))
  vv <- sapply(df[, controls, drop = FALSE], function(z) var(z, na.rm = TRUE))
  bad <- names(vv)[!is.finite(vv) | vv <= var_eps]
  kept <- setdiff(controls, bad)
  list(ok = TRUE, msg = "ok", kept = kept, drop = bad)
}

build_world_pts <- function(){
  world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    st_make_valid() |>
    mutate(iso3 = countrycode(iso_a3, "iso3c", "iso3c")) |>
    filter(!is.na(iso3), iso_a3 != "-99") |>
    select(iso3, geometry)
  st_point_on_surface(world_sf) |>
    mutate(lon = st_coordinates(geometry)[,1],
           lat = st_coordinates(geometry)[,2]) |>
    st_drop_geometry() |>
    select(iso3, lon, lat)
}

gtwr_kernel_weights <- function(dvec, bw, kernel = "bisquare", adaptive = TRUE) {
  dvec <- as.numeric(dvec)

  if (adaptive) {
    nn <- min(max(2L, as.integer(round(bw))), length(dvec))
    h <- sort(dvec, partial = nn)[nn]
    if (!is.finite(h) || h <= 0) {
      h <- max(dvec[dvec > 0], na.rm = TRUE)
    }
  } else {
    h <- bw
  }

  if (!is.finite(h) || h <= 0) return(rep(0, length(dvec)))

  if (kernel == "bisquare") {
    w <- ifelse(dvec < h, (1 - (dvec / h)^2)^2, 0)
  } else if (kernel == "gaussian") {
    w <- exp(-0.5 * (dvec / h)^2)
  } else {
    stop("Currently only 'bisquare' and 'gaussian' kernels are implemented in this diagnostic helper.")
  }

  w[!is.finite(w)] <- 0
  return(w)
}

weighted_r2 <- function(y, X, w) {
  X <- as.matrix(X)
  ok <- is.finite(y) & is.finite(w) & (w > 0)
  if (ncol(X) > 0) ok <- ok & rowSums(!is.finite(X)) == 0

  y <- y[ok]
  X <- X[ok, , drop = FALSE]
  w <- w[ok]

  if (length(y) < (ncol(X) + 1)) return(NA_real_)

  fit <- tryCatch(
    stats::lm.wfit(x = X, y = y, w = w),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)

  ybar <- sum(w * y) / sum(w)
  sse  <- sum(w * fit$residuals^2)
  sst  <- sum(w * (y - ybar)^2)

  if (!is.finite(sst) || sst <= 1e-12) return(NA_real_)

  r2 <- 1 - sse / sst
  r2 <- min(max(r2, 0), 0.999999999)
  return(r2)
}

gtwr_local_collin <- function(dat,
                              pred_vars,
                              stMat,
                              bw,
                              kernel = "bisquare",
                              adaptive = TRUE,
                              save_full_vdp = FALSE,
                              cn_thresh = 30,
                              vif_thresh = 10) {

  X_all <- as.matrix(dat[, pred_vars, drop = FALSE])
  n <- nrow(X_all)
  p <- ncol(X_all)

  cn_vec <- rep(NA_real_, n)

  vif_mat <- matrix(NA_real_, nrow = n, ncol = p,
                    dimnames = list(NULL, paste0("VIF_", pred_vars)))

  vdp_min_mat <- matrix(NA_real_, nrow = n, ncol = p,
                        dimnames = list(NULL, paste0("VDPmin_", pred_vars)))
  vdp_max_mat <- matrix(NA_real_, nrow = n, ncol = p,
                        dimnames = list(NULL, paste0("VDPmax_", pred_vars)))

  full_vdp_list <- if (save_full_vdp) vector("list", n) else NULL

  for (i in seq_len(n)) {
    w <- gtwr_kernel_weights(
      dvec = stMat[i, ],
      bw = bw,
      kernel = kernel,
      adaptive = adaptive
    )

    ok <- is.finite(w) & (w > 0)
    if (sum(ok) < max(8, p + 2)) next

    Xi <- X_all[ok, , drop = FALSE]
    wi <- w[ok]

    mu  <- colSums(Xi * wi) / sum(wi)
    Xc  <- sweep(Xi, 2, mu, "-")
    sdw <- sqrt(colSums(wi * Xc^2) / sum(wi))

    keep <- is.finite(sdw) & (sdw > 1e-10)
    if (sum(keep) < 2) next

    Xi2   <- Xi[, keep, drop = FALSE]
    vars2 <- pred_vars[keep]

    mu2   <- colSums(Xi2 * wi) / sum(wi)
    Xc2   <- sweep(Xi2, 2, mu2, "-")
    sdw2  <- sqrt(colSums(wi * Xc2^2) / sum(wi))
    Z     <- sweep(Xc2, 2, sdw2, "/")

    Zw <- Z * sqrt(wi)
    Cmat <- crossprod(Zw) / sum(wi)

    ee <- tryCatch(eigen(Cmat, symmetric = TRUE), error = function(e) NULL)
    if (is.null(ee)) next

    lam <- pmax(ee$values, 1e-12)
    V   <- ee$vectors

    ### Calculate the Belsley-style condition index.
    cn_vec[i] <- sqrt(max(lam) / min(lam))

    ### Calculate variance-decomposition proportions.
    phi <- sweep(V^2, 2, lam, "/")
    vdp <- phi / rowSums(phi)

    idx_small <- ncol(vdp)

    vdp_min_mat[i, paste0("VDPmin_", vars2)] <- vdp[, idx_small]
    vdp_max_mat[i, paste0("VDPmax_", vars2)] <- apply(vdp, 1, max)

    if (save_full_vdp) {
      full_vdp_list[[i]] <- list(
        vars = vars2,
        eigenvalues = lam,
        VDP = vdp
      )
    }

    ### Calculate local variance-inflation factors.
    if (length(vars2) >= 2) {
      for (j in seq_along(vars2)) {
        yj <- Xi2[, j]
        Xj <- Xi2[, -j, drop = FALSE]

        r2j <- weighted_r2(
          y = yj,
          X = cbind(Intercept = 1, Xj),
          w = wi
        )

        vif_val <- if (!is.finite(r2j)) {
          NA_real_
        } else if ((1 - r2j) <= 1e-10) {
          Inf
        } else {
          1 / (1 - r2j)
        }

        vif_mat[i, paste0("VIF_", vars2[j])] <- vif_val
      }
    }
  }

  point_df <- dat %>%
    dplyr::select(iso3, year) %>%
    dplyr::mutate(local_CN = cn_vec) %>%
    dplyr::bind_cols(as.data.frame(vif_mat)) %>%
    dplyr::bind_cols(as.data.frame(vdp_min_mat)) %>%
    dplyr::bind_cols(as.data.frame(vdp_max_mat))

  vif_cols <- grep("^VIF_", names(point_df), value = TRUE)

  point_df <- point_df %>%
    dplyr::mutate(
      flag_CN_gt30 = local_CN > cn_thresh,
      flag_any_VIF_gt10 = apply(
        as.matrix(dplyr::select(., dplyr::all_of(vif_cols))),
        1,
        function(z) any(z > vif_thresh, na.rm = TRUE)
      )
    )

  summary_df <- tibble::tibble(
    n_points = nrow(point_df),
    median_CN = median(point_df$local_CN, na.rm = TRUE),
    p95_CN    = quantile(point_df$local_CN, 0.95, na.rm = TRUE),
    max_CN    = max(point_df$local_CN, na.rm = TRUE),
    prop_CN_gt30 = mean(point_df$flag_CN_gt30, na.rm = TRUE),
    prop_any_VIF_gt10 = mean(point_df$flag_any_VIF_gt10, na.rm = TRUE)
  )

  out <- list(
    point = point_df,
    summary = summary_df
  )

  if (save_full_vdp) out$full_vdp <- full_vdp_list
  return(out)
}

make_point_panel <- function(df, years = YEARS_AVAILABLE){
  df |>
    filter(year %in% years) |>
    mutate(
      age_band = factor(as.character(age_group), levels = c("2-19","20-64","65+")),
      female   = as.integer(as.character(female))
    ) |>
    filter(!is.na(age_band), !is.na(female)) |>
    group_by(iso3, year, female, age_band) |>
    summarise(
      across(where(is.numeric) & !dplyr::matches("^female$") & !dplyr::matches("(_low|_high)$"),
             ~ mean(.x, na.rm = TRUE)),
      across(dplyr::matches("_low$"),  ~ quantile(.x, 0.025, na.rm = TRUE)),
      across(dplyr::matches("_high$"), ~ quantile(.x, 0.975, na.rm = TRUE)),
      .groups = "drop"
    )
}

zscore_cols <- function(dat, cols){
  dat %>% mutate(across(all_of(cols), ~ as.numeric(scale(.x)), .names = "{.col}"))
}

draw_one <- function(mu, lo, hi, mode = c("uniform","truncnorm")){
  mode <- match.arg(mode)
  if (is.na(lo) || is.na(hi) || lo >= hi) return(mu)
  if (mode == "uniform") runif(1, lo, hi) else {
    sd_est <- (hi - lo) / (2 * 1.96)
    rnorm(1, mean = mu, sd = sd_est)
  }
}

build_formula <- function(outcome, dmi_var, controls){
  rhs <- c(dmi_var, controls)
  as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}

gtwr_for_dmi_all_strata <- function(point_panel,
                                    dmi_key,
                                    years = YEARS_AVAILABLE,
                                    outcome = "HALE",
                                    outcome_low = "HALE_low", outcome_high = "HALE_high",
                                    controls = c("female","SDI","urban_rate","pop_density","health_exp_pc",
                                                 "pm25","tmin_mean","tmax_mean","prec_mean","elev"),
                                    standardize = c("x","none","xy"),
                                    bw_approach = c("CV","AICc"),
                                    parallel_boot = TRUE,
                                    B = 300,
                                    sample_mode = c("uniform","truncnorm"),
                                    fix_bandwidth = TRUE,
                                    lamda = LAMBDA, ksi = KSI, kernel = KERNEL, adaptive = ADAPT,
                                    save_local_collin = TRUE,
                                    save_full_vdp = FALSE,
                                    cn_thresh = 30,
                                    vif_thresh = 10,
                                    verbose = TRUE){

  standardize <- match.arg(standardize)
  bw_approach <- match.arg(bw_approach)
  sample_mode <- match.arg(sample_mode)

  world_pts <- build_world_pts()

  strata <- tidyr::expand_grid(
    age_band = factor(c("2-19","20-64","65+"), levels = c("2-19","20-64","65+"))
  )

  local_list    <- list()
  stratum_list  <- list()
  bw_list       <- list()
  diag_list     <- list()
  bootcoef_list <- list()
  collin_point_list   <- list()
  collin_summary_list <- list()
  collin_full_list    <- list()

  for (i in seq_len(nrow(strata))){
    ab <- strata$age_band[i]

    dmi_var  <- paste0(dmi_key, "_DMI")
    dmi_low  <- paste0(dmi_key, "_DMI_low")
    dmi_high <- paste0(dmi_key, "_DMI_high")

    needed_cols <- c(outcome, outcome_low, outcome_high,
                     dmi_var, dmi_low, dmi_high, controls)
    dat <- point_panel |>
      filter(age_band == ab, year %in% years) |>
      left_join(world_pts, by = "iso3") |>
      select(iso3, year, lon, lat, all_of(needed_cols)) |>
      drop_na(lon, lat, all_of(needed_cols))

    if (nrow(dat) < 30) {
      if (verbose) message("[", dmi_key, "] age=", ab, ": too few obs (", nrow(dat), "). Skip.")
      next
    }

    x_cols <- c(dmi_var, setdiff(controls, "female")); y_cols <- c(outcome)
    dat_std <- dat
    if (standardize %in% c("x","xy")){
      dat_std <- zscore_cols(dat_std, x_cols)
      mu <- mean(dat[[dmi_var]], na.rm=TRUE); sd0 <- sd(dat[[dmi_var]], na.rm=TRUE)
      if (sd0 > 0) {
        dat_std[[dmi_low]]  <- (dat[[dmi_low]]  - mu)/sd0
        dat_std[[dmi_high]] <- (dat[[dmi_high]] - mu)/sd0
      }
    }
    if (standardize == "xy"){
      dat_std <- zscore_cols(dat_std, y_cols)
      mu <- mean(dat[[outcome]], na.rm=TRUE); sd0 <- sd(dat[[outcome]], na.rm=TRUE)
      if (sd0 > 0) {
        dat_std[[outcome_low]]  <- (dat[[outcome_low]]  - mu)/sd0
        dat_std[[outcome_high]] <- (dat[[outcome_high]] - mu)/sd0
      }
    }
    dat_std <- as.data.frame(dat_std, stringsAsFactors = FALSE)

    chk <- check_and_trim_vars(dat_std, controls = controls, min_n = 30, var_eps = 1e-12)
    controls_use <- chk$kept
    if (length(chk$drop) > 0 && verbose) {
      message("  dropped near-constant controls: ", paste(chk$drop, collapse = ", "))
    }

    form      <- build_formula(outcome, dmi_var, controls_use)
    var_names <- c(dmi_var, controls_use)

    ### Convert the analysis data to a spatial points object.
    spdat <- dat_std
    sp::coordinates(spdat) <- ~ lon + lat
    sp::proj4string(spdat) <- sp::CRS("+proj=longlat +datum=WGS84")
    spdat@data <- as.data.frame(spdat@data, stringsAsFactors = FALSE)

    stMat <- GWmodel::st.dist(
      dp.locat = sp::coordinates(spdat),
      obs.tv = spdat@data$year,
      p = 2,
      longlat = TRUE,
      lamda = lamda,
      t.units = "year",
      ksi = ksi
    )

    if (verbose) message("[", dmi_key, "] age=", ab, " : selecting bandwidth (", bw_approach, ") ...")
    bw_base <- tryCatch(
      GWmodel::bw.gtwr(formula = form, data = spdat,
                       obs.tv = spdat@data$year,
                       approach = bw_approach,
                       kernel = kernel, adaptive = adaptive,
                       p = 2, theta = 0, longlat = TRUE,
                       lamda = lamda, t.units = "year",
                       ksi = ksi, st.dMat = stMat),
      error = function(e) NULL
    )
    if (is.null(bw_base)) { if (verbose) message("  bandwidth selection failed. Skip."); next }

    fit_base <- tryCatch(
      GWmodel::gtwr(formula = form, data = spdat,
                    obs.tv = spdat@data$year, reg.tv = spdat@data$year,
                    st.bw = bw_base,
                    kernel = kernel, adaptive = adaptive,
                    p = 2, theta = 0, longlat = TRUE,
                    lamda = lamda, t.units = "year",
                    ksi = ksi, st.dMat = stMat),
      error = function(e) NULL
    )
    if (is.null(fit_base)) { if (verbose) message("  base GTWR fit failed. Skip."); next }

    diag_base <- fit_base$GTW.diagnostic
    if (verbose) {
      message("  AICc=", safe_round(diag_base$AICc, 2), " | R2=", safe_round(diag_base$gw.R2, 3),
              " | R2.adj=", safe_round(diag_base$gwR2.adj, 3))
    }

    if (save_local_collin) {
      if (verbose) message("  computing local collinearity diagnostics ...")

      collin_base <- gtwr_local_collin(
        dat = dat_std,
        pred_vars = c(dmi_var, controls_use),
        stMat = stMat,
        bw = bw_base,
        kernel = kernel,
        adaptive = adaptive,
        save_full_vdp = save_full_vdp,
        cn_thresh = cn_thresh,
        vif_thresh = vif_thresh
      )

      collin_point <- collin_base$point %>%
        dplyr::mutate(DMI = dmi_key, age_band = ab)

      collin_summary <- collin_base$summary %>%
        dplyr::mutate(DMI = dmi_key, age_band = ab)

      if (verbose) {
        message("    median CN = ", round(collin_summary$median_CN, 2),
                " | p95 CN = ", round(collin_summary$p95_CN, 2),
                " | prop(CN>30) = ", round(collin_summary$prop_CN_gt30, 3),
                " | prop(any VIF>10) = ", round(collin_summary$prop_any_VIF_gt10, 3))
      }
    } else {
      collin_point <- NULL
      collin_summary <- NULL
      collin_base <- NULL
    }

    n <- nrow(spdat)
    beta_mat_dmi   <- matrix(NA_real_, nrow = B, ncol = n)
    spacemed_dmi   <- numeric(B)
    coef_store     <- vector("list", B)
    space_med_mat  <- matrix(NA_real_, nrow = B, ncol = length(var_names))
    colnames(space_med_mat) <- var_names

    boot_worker <- function(b){
      tmp <- spdat
      tmp@data <- as.data.frame(tmp@data, stringsAsFactors = FALSE)

      tmp@data[[dmi_var]] <- mapply(
        draw_one,
        mu = dat_std[[dmi_var]],
        lo = dat_std[[paste0(dmi_key, "_DMI_low")]],
        hi = dat_std[[paste0(dmi_key, "_DMI_high")]],
        MoreArgs = list(mode = sample_mode)
      )
      tmp@data[[outcome]] <- mapply(
        draw_one,
        mu = dat_std[[outcome]],
        lo = dat_std[[outcome_low]],
        hi = dat_std[[outcome_high]],
        MoreArgs = list(mode = sample_mode)
      )

      bw_use <- if (fix_bandwidth) bw_base else
        tryCatch(
          GWmodel::bw.gtwr(
            formula = form, data = tmp, obs.tv = tmp@data$year,
            approach = bw_approach, kernel = kernel, adaptive = adaptive,
            p = 2, theta = 0, longlat = TRUE, lamda = lamda, t.units = "year",
            ksi = ksi, st.dMat = stMat
          ),
          error = function(e) NULL
        )
      if (is.null(bw_use)) return(NULL)

      fit_b <- tryCatch(
        GWmodel::gtwr(
          formula = form, data = tmp,
          obs.tv = tmp@data$year, reg.tv = tmp@data$year,
          st.bw = bw_use, kernel = kernel, adaptive = adaptive,
          p = 2, theta = 0, longlat = TRUE, lamda = lamda, t.units = "year",
          ksi = ksi, st.dMat = stMat
        ),
        error = function(e) NULL
      )
      if (is.null(fit_b)) return(NULL)

      coef_mat <- as.matrix(fit_b$SDF@data[, var_names, drop = FALSE])

      dmi_vec <- coef_mat[, dmi_var]
      if (all(!is.finite(dmi_vec))) {
        nm <- grep(paste0("^", dmi_key), names(fit_b$SDF@data), value = TRUE)[1]
        dmi_vec <- as.numeric(fit_b$SDF@data[[nm]])
      }

      return(list(
        coef_mat = coef_mat,
        dmi_vec  = dmi_vec,
        dmi_med  = median(dmi_vec, na.rm = TRUE)
      ))
    }

    if (parallel_boot) {
      boot_res <- future.apply::future_lapply(seq_len(B), boot_worker, future.seed = TRUE)
    } else {
      boot_res <- lapply(seq_len(B), boot_worker)
    }

    ok_idx <- which(vapply(boot_res, function(x) !is.null(x) && is.matrix(x$coef_mat), logical(1)))
    if (length(ok_idx) == 0) {
      if (verbose) message("  no successful bootstrap fits; skip this stratum.")
      next
    }

    for (b in ok_idx) {
      coef_store[[b]]   <- boot_res[[b]]$coef_mat
      beta_mat_dmi[b, ] <- boot_res[[b]]$dmi_vec
      spacemed_dmi[b]   <- boot_res[[b]]$dmi_med
      space_med_mat[b, ] <- apply(boot_res[[b]]$coef_mat, 2, median, na.rm = TRUE)
    }

    coef_store_ok <- boot_res[ok_idx]
    coef_arr <- simplify2array(lapply(coef_store_ok, `[[`, "coef_mat"))

    beta_med_loc <- apply(coef_arr, c(1,2), median,   na.rm = TRUE)
    beta_lo_loc  <- apply(coef_arr, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
    beta_hi_loc  <- apply(coef_arr, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
    beta_var_loc <- apply(coef_arr, c(1,2), var, na.rm = TRUE)

    local_robust_all <- tibble(
      DMI = dmi_key, age_band = ab,
      iso3 = spdat@data$iso3, year = spdat@data$year,
      female = spdat@data$female
    ) %>%
      bind_cols(as.data.frame(beta_med_loc) |> setNames(paste0(colnames(beta_med_loc), "_med"))) %>%
      bind_cols(as.data.frame(beta_lo_loc)  |> setNames(paste0(colnames(beta_lo_loc),  "_lo")))  %>%
      bind_cols(as.data.frame(beta_hi_loc)  |> setNames(paste0(colnames(beta_hi_loc),  "_hi"))) %>%
      bind_cols(as.data.frame(beta_var_loc) |> setNames(paste0(colnames(beta_var_loc), "_within_var"))) %>%
      tidyr::pivot_longer(
        cols = matches("_(med|lo|hi|within_var)$"),
        names_to = c("term",".value"),
        names_pattern = "^(.*)_(med|lo|hi|within_var)$"
      )

    stratum_robust_all <- as_tibble(space_med_mat) %>%
      mutate(.b = row_number()) %>%
      tidyr::pivot_longer(cols = all_of(var_names), names_to = "term", values_to = "space_med") %>%
      group_by(term) %>%
      summarise(
        beta_med = median(space_med, na.rm = TRUE),
        beta_lo  = quantile(space_med, 0.025, na.rm = TRUE),
        beta_hi  = quantile(space_med, 0.975, na.rm = TRUE),
        within_var = var(space_med, na.rm = TRUE),
        B_ok     = sum(is.finite(space_med)),
        .groups  = "drop"
      ) %>%
      mutate(DMI = dmi_key, age_band = ab,
             AICc = safe_round(diag_base$AICc, 2),
             R2   = safe_round(diag_base$gw.R2, 3),
             R2.adj = safe_round(diag_base$gwR2.adj, 3),
             bw   = bw_base)

    key <- paste0(dmi_key, "_", ab)
    local_list[[key]]    <- local_robust_all
    stratum_list[[key]]  <- stratum_robust_all
    bw_list[[key]]       <- bw_base
    diag_list[[key]]     <- diag_base
    bootcoef_list[[key]] <- coef_store

    if (save_local_collin) {
      collin_point_list[[key]]   <- collin_point
      collin_summary_list[[key]] <- collin_summary
      if (save_full_vdp && !is.null(collin_base$full_vdp)) {
        collin_full_list[[key]] <- collin_base$full_vdp
      }
    }
  }

    list(
    local_robust   = bind_rows(local_list),
    stratum_robust = bind_rows(stratum_list),
    bw             = bw_list,
    diag           = diag_list,
    boot_coef      = bootcoef_list,
    local_collin   = dplyr::bind_rows(collin_point_list),
    collin_summary = dplyr::bind_rows(collin_summary_list),
    local_collin_full = collin_full_list
  )
}

# 8. Fit GTWR models in all 20 completed datasets =================================================
## 8.1 Define the repeated-analysis settings ================
GTWR_PATTERNS <- c("AHEI", "MSDPS", "hPDI", "DASH", "DRRD", "WCRF", "PHDI", "DBI")
GTWR_CONTROLS <- c(
  "female", "SDI", "urban_rate", "pop_density", "health_exp_pc",
  "pm25", "tmin_mean", "tmax_mean", "prec_mean", "elev"
)
GTWR_BOOT_B <- as.integer(Sys.getenv("GDD_GTWR_BOOT_B", "500"))
GTWR_WORKERS <- as.integer(Sys.getenv("GDD_GTWR_WORKERS", "6"))
GTWR_MI_DIR <- here("GTWR_MI20")
dir.create(GTWR_MI_DIR, recursive = TRUE, showWarnings = FALSE)

panel_boot_files <- file.path(MI_DIR, sprintf("panel_boot_imputed_%02d.rds", seq_len(MI_M)))
stopifnot(all(file.exists(panel_boot_files)))

## 8.2 Apply the same GTWR function to each imputed dataset ================
plan(multisession, workers = GTWR_WORKERS)

for (imputation_id in seq_len(MI_M)) {
  panel_boot_i <- readRDS(panel_boot_files[[imputation_id]])
  stopifnot(all(c(
    "iso3", "year", "age_group", "female", "HALE", "HALE_low", "HALE_high",
    "HALE_LE", "HALE_LE_low", "HALE_LE_high"
  ) %in% names(panel_boot_i)))
  point_panel_i <- make_point_panel(panel_boot_i, years = YEARS_AVAILABLE)

  for (pattern_name in GTWR_PATTERNS) {
    result_file <- file.path(
      GTWR_MI_DIR,
      sprintf("GTWR_result_%s_imp%02d.rds", pattern_name, imputation_id)
    )

    if (file.exists(result_file)) {
      message("Using existing result: ", result_file)
      next
    }

    message(
      "Fitting ", pattern_name, " in imputation ", imputation_id,
      " of ", MI_M, "."
    )
    result_i <- gtwr_for_dmi_all_strata(
      point_panel_i,
      dmi_key = pattern_name,
      years = YEARS_AVAILABLE,
      outcome = "HALE",
      outcome_low = "HALE_low",
      outcome_high = "HALE_high",
      controls = GTWR_CONTROLS,
      standardize = "x",
      bw_approach = "CV",
      B = GTWR_BOOT_B,
      sample_mode = "truncnorm",
      fix_bandwidth = TRUE,
      lamda = 0.5,
      ksi = 0.5,
      kernel = "bisquare",
      adaptive = TRUE,
      save_local_collin = TRUE,
      save_full_vdp = FALSE,
      cn_thresh = 30,
      vif_thresh = 10,
      verbose = TRUE,
      parallel_boot = TRUE
    )
    saveRDS(result_i, result_file, compress = "xz")
  }
}

plan(sequential)

# 9. Pool GTWR estimates with Rubin's rules =======================================================
## 9.1 Collect the estimates and within-imputation variances ================
result_index <- tidyr::crossing(
  imputation_id = seq_len(MI_M),
  pattern = GTWR_PATTERNS
) %>%
  mutate(
    result_file = file.path(
      GTWR_MI_DIR,
      sprintf("GTWR_result_%s_imp%02d.rds", pattern, imputation_id)
    )
  )
stopifnot(all(file.exists(result_index$result_file)))

local_mi_results <- purrr::pmap_dfr(result_index, function(imputation_id, pattern, result_file) {
  readRDS(result_file)$local_robust %>%
    mutate(.imp = imputation_id, .before = 1)
})

stratum_mi_results <- purrr::pmap_dfr(result_index, function(imputation_id, pattern, result_file) {
  readRDS(result_file)$stratum_robust %>%
    mutate(.imp = imputation_id, .before = 1)
})

saveRDS(local_mi_results, file.path(GTWR_MI_DIR, "GTWR_local_all_imputations.rds"), compress = "xz")
saveRDS(stratum_mi_results, file.path(GTWR_MI_DIR, "GTWR_stratum_all_imputations.rds"), compress = "xz")

## 9.2 Implement Rubin's rules ================
rubin_pool <- function(estimates, within_variances, conf_level = 0.95) {
  valid <- is.finite(estimates) & is.finite(within_variances) & within_variances >= 0
  estimates <- estimates[valid]
  within_variances <- within_variances[valid]
  m_used <- length(estimates)
  if (!m_used) {
    return(tibble(
      estimate = NA_real_, std_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
      df = NA_real_, within_var = NA_real_, between_var = NA_real_, total_var = NA_real_,
      fraction_missing_information = NA_real_, m_used = 0L
    ))
  }

  pooled_estimate <- mean(estimates)
  within_var <- mean(within_variances)
  between_var <- if (m_used > 1) stats::var(estimates) else 0
  total_var <- within_var + (1 + 1 / m_used) * between_var
  relative_increase <- if (within_var > 0) {
    (1 + 1 / m_used) * between_var / within_var
  } else if (between_var > 0) {
    Inf
  } else {
    0
  }
  df <- if (m_used <= 1 || relative_increase == 0) {
    Inf
  } else if (is.infinite(relative_increase)) {
    m_used - 1
  } else {
    (m_used - 1) * (1 + 1 / relative_increase)^2
  }
  critical_value <- if (is.finite(df)) {
    stats::qt(1 - (1 - conf_level) / 2, df = df)
  } else {
    stats::qnorm(1 - (1 - conf_level) / 2)
  }
  standard_error <- sqrt(total_var)

  tibble(
    estimate = pooled_estimate,
    std_error = standard_error,
    conf_low = pooled_estimate - critical_value * standard_error,
    conf_high = pooled_estimate + critical_value * standard_error,
    df = df,
    within_var = within_var,
    between_var = between_var,
    total_var = total_var,
    fraction_missing_information = if (relative_increase == 0) {
      0
    } else if (is.infinite(relative_increase)) {
      1
    } else {
      (relative_increase + 2 / (df + 3)) / (relative_increase + 1)
    },
    m_used = m_used
  )
}

## 9.3 Pool local and stratum-level GTWR coefficients ================
pooled_local_results <- local_mi_results %>%
  group_by(DMI, age_band, iso3, year, female, term) %>%
  group_modify(~ rubin_pool(.x$med, .x$within_var)) %>%
  ungroup()

pooled_stratum_results <- stratum_mi_results %>%
  group_by(DMI, age_band, term) %>%
  group_modify(~ rubin_pool(.x$beta_med, .x$within_var)) %>%
  ungroup()

diagnostic_summary <- stratum_mi_results %>%
  summarise(
    bandwidth_median = median(bw, na.rm = TRUE),
    bandwidth_min = min(bw, na.rm = TRUE),
    bandwidth_max = max(bw, na.rm = TRUE),
    AICc_median = median(AICc, na.rm = TRUE),
    R2_median = median(R2, na.rm = TRUE),
    adjusted_R2_median = median(R2.adj, na.rm = TRUE),
    .by = c(DMI, age_band)
  )

saveRDS(pooled_local_results, file.path(GTWR_MI_DIR, "GTWR_local_Rubin_pooled.rds"), compress = "xz")
saveRDS(pooled_stratum_results, file.path(GTWR_MI_DIR, "GTWR_stratum_Rubin_pooled.rds"), compress = "xz")
readr::write_csv(pooled_local_results, file.path(GTWR_MI_DIR, "GTWR_local_Rubin_pooled.csv"))
readr::write_csv(pooled_stratum_results, file.path(GTWR_MI_DIR, "GTWR_stratum_Rubin_pooled.csv"))
readr::write_csv(diagnostic_summary, file.path(GTWR_MI_DIR, "GTWR_diagnostics_across_imputations.csv"))

cat("All GTWR models and Rubin pooling steps are complete.\n")
cat("Results were saved under:", GTWR_MI_DIR, "\n")
cat("Pipeline completed:", format(Sys.time()), "\n")
}, finally = {
  if (sink.number(type = "message") != prior_message_sink) sink(type = "message")
  if (sink.number(type = "output") > prior_output_sinks) sink(type = "output")
  if (isOpen(pipeline_log_connection)) close(pipeline_log_connection)
})
