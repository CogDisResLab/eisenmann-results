# Reshape the results of the Creedenzymatic analysis to a format that can be used annotating
# "favorite" kinases

library(tidyverse)

result_files <- list.files("results", "creedenzymatic", full.names = TRUE) |>
  set_names(~ str_remove(basename(.x), fixed("_creedenzymatic.csv"))) |>
  keep(~ str_detect(.x, fixed("Pat27")) | str_detect(.x, fixed("Pat31")))

reshape_data <- function(file_path) {
  read_csv(file_path) |>
    select(hgnc_symbol, Perc, Method) |>
    mutate(
      Method = str_remove(Method, fixed("-"))
    ) |>
    pivot_wider(
      names_from = Method,
      values_from = Perc, values_fn = unique
    ) |>
    mutate(
      AverageScore = (PTMSEA + KEA3 + UKA + KRSA) / 4L,
      Presence = sign(PTMSEA) + sign(KEA3) +
        sign(UKA) + sign(KRSA)
    ) |>
    arrange(AverageScore)
}

results <- map(result_files, reshape_data) |>
  imap(
    ~ write_csv(.x, file = file.path("results", paste0(.y, "_reshaped.csv")))
  )
