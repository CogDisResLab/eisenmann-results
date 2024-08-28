# Collect and compile all the signal data for all four states

library(tidyverse)

data_files <- list.files("results", "signal", full.names = TRUE) |>
  set_names(~ basename(.x) |> str_remove(".csv"))

combined_data <- data_files |>
  map(read_csv) |>
  bind_rows(.id = "filename") |>
  mutate(Chip = str_extract(filename, "[SP]TK")) |>
  nest(.by = Chip) |>
  mutate(
    data = map(
      data,
      ~ select(
        .x,
        SampleName, Peptide, slope
      )
    ),
    pivoted = map(
      data,
      ~ pivot_wider(
        .x,
        names_from = SampleName,
        values_from = slope,
        values_fill = NA
      )
    ),
    written = map2(
      pivoted, Chip,
      ~ write_csv(.x,
        file = file.path("results", paste0(.y, "_collated_signal_data.csv"))
      )
    )
  )
