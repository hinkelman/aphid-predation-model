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
#' Movement over the tracked part of post-handling time (up to 90 min; NA for
#' 13 untracked videos): `move_speed` (mm/s while moving) and `activity`
#' (proportion of tracked time moving).
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
      leaf_area = LeafArea,
      move_speed = MoveSpeed,
      activity = Activity / 100
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

#' Larval diet experiment (Sep 2008; dissertation Ch. 1)
#'
#' Larvae reared individually from hatching on bean (B), pea (P) or a 50:50
#' mix (M) of size-matched aphids (~0.9 mg), supplied daily ad libitum.
#' Returns one row per larva per check interval k = (check k, check k + 1]:
#' * `stage`: stage at check k; `next_stage`: stage at check k + 1
#' * `pea`, `bean`: aphids killed during the interval (the "eaten" counts
#'   recorded at check k + 1; dead aphids with no piercing are excluded)
#' * `died`: larva found dead at check k + 1
#' * `censored`: larva removed (disabled) at check k + 1
#' The final record of each larva (death, adult, or removal) has no interval.
#' Intake is 0 in non-feeding intervals (pre-pupa, pupa), which were not
#' recorded, and NA for a few recording gaps (mixed-diet larvae 13 and 14).
read_diet_experiment <- function() {
  develop <- read_mac_csv(data_path("develop.csv")) |>
    janitor::clean_names() |>
    dplyr::select(id, diet, day, stage, disabled, dead)

  kills <- read_mac_csv(data_path("feed.csv")) |>
    janitor::clean_names() |>
    dplyr::filter(status == "eaten") |>
    dplyr::mutate(aphid = tolower(aphid), day = day - 1) |> # recorded at k + 1
    tidyr::pivot_wider(id_cols = c(id, day), names_from = aphid, values_from = count)

  develop |>
    dplyr::arrange(id, day) |>
    dplyr::mutate(
      next_stage = dplyr::lead(stage),
      died = dplyr::lead(dead),
      censored = dplyr::lead(disabled),
      .by = id
    ) |>
    dplyr::filter(!is.na(next_stage)) |>
    dplyr::left_join(kills, by = c("id", "day")) |>
    dplyr::mutate(
      # Non-feeding intervals were not recorded: the pre-pupa at the end of
      # L4 (interval ending in pupation) and pupae. Kills there are 0.
      non_feeding = stage %in% c("Pupa", "Adult") | (stage == "L4" & next_stage == "Pupa"),
      pea = dplyr::if_else(non_feeding & is.na(pea), 0, pea),
      bean = dplyr::if_else(non_feeding & is.na(bean), 0, bean),
      # single-species diets: the absent species is 0, not missing
      pea = dplyr::if_else(diet == "B" & !is.na(bean), 0, pea),
      bean = dplyr::if_else(diet == "P" & !is.na(pea), 0, bean),
      diet = factor(diet, levels = c("P", "M", "B")),
      stage = factor(stage, levels = c("L1", "L2", "L3", "L4", "Pupa", "Adult")),
      next_stage = factor(next_stage, levels = levels(stage))
    ) |>
    dplyr::select(id, diet, day, stage, next_stage, pea, bean, died, censored)
}

#' Diet timing experiment: 4th instars switched to bean aphids (B) or no food
#' (S) on day 1-4 of the instar (`timing`); controls (P) stayed on pea.
#' `l4_mass` is mass (mg) at the switch; `adult_mass` in mg.
read_diet_timing <- function() {
  suppressWarnings(suppressMessages(readxl::read_excel(data_path("DietTimingData.xlsx")))) |>
    dplyr::select(ID:Eclosed) |>
    janitor::clean_names() |>
    dplyr::mutate(
      dplyr::across(c(timing, fourth, pupa, develop, trt_to_pupa, trt_to_adult, t_weight, a_weight),
                    \(x) suppressWarnings(as.numeric(x))),
      trt = factor(trt, levels = c("P", "B", "S")),
      l4_mass = t_weight * 1000,
      adult_mass = a_weight * 1000
    ) |>
    dplyr::select(id, block, trt, timing, l4_days = fourth, pupa_days = pupa, develop,
                  trt_to_pupa, l4_mass, adult_mass, sex, pupated, eclosed)
}
