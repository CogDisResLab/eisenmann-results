# Creedenzymatic Analysis

library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(creedenzymatic)

process_creedenzymatic <-
  function(krsa_path, uka_path, peptide_path) {
    krsa_data <- readr::read_csv(krsa_path, show_col_types = FALSE) |>
      select(Kinase, Score = AvgZ) |>
      read_krsa(trns = "abs", sort = "desc")

    uka_data <- readr::read_tsv(uka_path, show_col_types = FALSE) |>
      dplyr::select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
      read_uka(trns = "abs", sort = "desc")

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
      UKA_df = uka_data,
      KEA3_df = kea3_data,
      PTM_SEA_df = ptmsea_data
    )

    combined
  }

krsa_files <- list.files("results", "krsa", full.names = TRUE) |>
  keep(~ str_detect(.x, "Pat27|Pat31")) |>
  sort()

uka_files <- c(
  "kinome_data/PTK/UKA/Pat27/Summaryresults 20231025-1512.txt",
  "kinome_data/STK/UKA/Pat27/Summaryresults 20231024-1346.txt",
  "kinome_data/PTK/UKA/Pat31/Summaryresults 20231025-1527.txt",
  "kinome_data/STK/UKA/Pat31/Summaryresults 20231024-1356.txt"
)

peptide_files <- list.files("results", "dpp", full.names = TRUE) |>
  keep(~ str_detect(.x, "Pat27|Pat31")) |>
  sort()

result_names <- str_extract(krsa_files, "full_(.*)\\.csv", 1)

result <-
  list(
    krsa_path = krsa_files,
    uka_path = uka_files,
    peptide_path = peptide_files
  ) |>
  pmap(process_creedenzymatic) |>
  set_names(result_names) |>
  imap_dfr(~ write_csv(
    .x,
    str_glue("results/{.y}_creedenzymatic.csv")
  ), .id = "Comparison")
