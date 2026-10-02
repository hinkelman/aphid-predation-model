# Data import and cleaning ------------------------------------------------
#
# Raw files in data/ are treated as read-only. These functions return tidy
# tibbles with snake_case names and documented fixes applied.

data_path <- function(...) file.path("data", ...)

# Several CSVs were saved with classic Mac (CR-only) line endings, which
# read.csv/readr misread as a single line. Normalize line endings first.
read_mac_csv <- function(file, ...) {
  txt <- readChar(file, file.size(file), useBytes = TRUE)
  txt <- gsub("\r\n?", "\n", txt)
  readr::read_csv(I(txt), show_col_types = FALSE, ...)
}

#' Aphid clip-cage life tables (Dec 2009 - Jan 2010)
#'
#' One row per aphid per day. `day` is age in days (day 1 = first day as L1).
#' Fixes applied:
#' * Bean start dates recorded as 2010-12-15 are a typo for 2009-12-15; dates
#'   are rebuilt from the true start date and `day`.
#' * `status2` is renamed `died` (1 = death observed on last day, 0 = censored).
#' * Pea aphid 23 lived 24 days, reached adulthood, and never reproduced; it is
#'   excluded as an abnormal individual (`drop_sterile = FALSE` keeps it).
read_clip_cage <- function(drop_sterile = TRUE) {
  raw <- suppressWarnings(suppressMessages(
    readxl::read_excel(data_path("ClipCageData.xlsx"))
  ))

  raw |>
    dplyr::select(Block:Status2) |>
    janitor::clean_names() |>
    dplyr::rename(died = status2) |>
    dplyr::mutate(
      aphid = factor(aphid, levels = c("pea", "bean")),
      id = as.integer(id),
      day = as.integer(day),
      date = as.Date("2009-12-15") + day - 1L,
      stage = factor(stage, levels = c("L1", "L2", "L3", "L4", "A")),
      offspring = as.integer(offspring)
    ) |>
    dplyr::filter(!(drop_sterile & aphid == "pea" & id == 23L)) |>
    dplyr::arrange(aphid, id, day)
}

#' One row per aphid: lifespan (days) and fate.
clip_cage_lifespans <- function(clip_cage = read_clip_cage()) {
  clip_cage |>
    dplyr::summarise(
      lifespan = max(day),
      died = dplyr::first(died),
      total_offspring = sum(offspring, na.rm = TRUE),
      .by = c(aphid, id)
    )
}
