# Creedenzymatic Analysis

library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(creedenzymatic)

process_creedenzymatic <-
  function(krsa_path, peptide_path) {
    krsa_data <- readr::read_csv(krsa_path, show_col_types = FALSE) |>
      select(Kinase, Score = AvgZ) |>
      read_krsa(trns = "abs", sort = "desc")

    peptide_data <-
      readr::read_csv(peptide_path, show_col_types = FALSE) |>
      select(Peptide, Score = totalMeanLFC)

    kea3_data <-
      read_kea(
        peptide_data,
        sort = "asc",
        trns = "abs",
        method = "MeanRank",
        lib = "kinase-substrate"
      )

    ptmsea_data <-
      read_ptmsea(peptide_data)

    combined <- combine_tools(
      KRSA_df = krsa_data,
      KEA3_df = kea3_data,
      PTM_SEA_df = ptmsea_data
    )

    combined
  }

krsa_files <- list.files("results", "krsa", full.names = TRUE) |>
  keep(~ str_detect(.x, "pat9")) |>
  sort()

peptide_files <- list.files("results", "dpp", full.names = TRUE) |>
  keep(~ str_detect(.x, "pat9")) |>
  sort()

result_names_a <- str_extract(krsa_files, "results/(.*)_ptk-krsa_table_full_PLNB\\-(.*)\\.csv", 1)
result_names_b <- str_extract(krsa_files, "results/(.*)_ptk-krsa_table_full_PLNB\\-(.*)\\.csv", 2)

result_names <- str_c(result_names_a, result_names_b, sep = "_")

result <-
  list(
    krsa_path = krsa_files,
    peptide_path = peptide_files
  ) |>
  pmap(process_creedenzymatic) |>
  set_names(result_names) |>
  imap_dfr(~ write_csv(
    .x,
    str_glue("results/{.y}_creedenzymatic.csv")
  ), .id = "Comparison")
