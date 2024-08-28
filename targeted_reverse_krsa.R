# Generate reverse KRSA with a targeted list of kinases

library(tidyverse)
library(KRSA)

dpp_files <- list.files("results", "dpp", full.names = TRUE) |>
  set_names(~ str_remove(basename(.x), fixed("_dpp.csv"))) |>
  keep(~ str_detect(.x, fixed("Pat27")) | str_detect(.x, fixed("Pat31")))

reverse_krsa <- function(file) {
  data <- read_csv(file)

  plot <- krsa_reverse_krsa_plot(
    chipCov = KRSA_coverage_PTK_PamChip_86402_v1,
    lfc_table = data,
    kinases = c("EGFR", "VEGFR", "ABL", "PDGFR", "FGFR"),
    lfc_thr = 0.2
  )
}

plots <- map(dpp_files, reverse_krsa) |>
  imap(
    ~ ggsave(.x, file = file.path("results", paste0(.y, "_targeted_reverse_krsa.png")))
  )



mean_reverse_krsa_values <- function(file) {
  data <- read_csv(file) |>
    left_join(
      {
        KRSA_coverage_PTK_PamChip_86402_v1 |> rename(Kinase = Kin, Peptide = Substrates)
      },
      by = "Peptide"
    ) |>
    select(Kinase, Peptide, totalMeanLFC) |>
    distinct() |>
    nest(.by = Kinase) |>
    mutate(mean = map_dbl(data, ~ mean(.x[["totalMeanLFC"]], na.rm = TRUE)))
}

mean_values <- map(dpp_files, mean_reverse_krsa_values) |>
  imap(
    ~ write_csv(.x, file = file.path("results", paste0(.y, "_mean_reverse_krsa.csv")))
  )
