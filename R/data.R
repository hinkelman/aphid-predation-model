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

#' Single-encounter behavior trials (2008; dissertation Ch. 2)
#'
#' 4th-instar larvae (`l4_day` = 1 or 2 days into the instar) starved for
#' `starve` hours, given one adult bean aphid or a size-matched pea aphid on
#' a fava leaf. Times in minutes:
#' * `handle`: from securing the aphid until moving away from the feeding
#'   site (`handle_end` = 1 if observed).
#' * `post_handle`: from the end of handling until leaving the leaf
#'   (`left` = 1 if observed). One trial has post_handle = 0.
read_behavior_trials <- function() {
  read_mac_csv(data_path("handle_depart_move.csv")) |>
    dplyr::transmute(
      id = ID,
      aphid = factor(Aphid, levels = c("pea", "bean")),
      l4_day = Age,
      starve = Starve,
      mass = Mass * 1000, # mg
      handle = Handle,
      handle_end = HandleEnd,
      partial = Partial,
      post_handle = Depart,
      left = Left,
      brush = Brush,
      harass = Harass,
      leaf_area = LeafArea
    )
}

#' Handling time by aphid age (2005)
#'
#' Larvae starved ~2-4 h were offered one aphid of known age (days). Handling
#' ended when the larva moved away from the feeding site (a censored trial is
#' annotated "finished eating, but had not moved away from feeding site").
#' `rejected` = 1 when the larva would not attack or eat the aphid; those rows
#' have no handling time.
read_handle_age <- function() {
  read_mac_csv(data_path("HandleAge.csv")) |>
    dplyr::transmute(
      id = ID,
      aphid = factor(tolower(Aphid.Species), levels = c("pea", "bean")),
      aphid_age = Aphid.Age,
      starve = Starve,
      handle = Handle,
      handle_end = 1L - Censored,
      partial = Partial,
      rejected = Rejected
    )
}
