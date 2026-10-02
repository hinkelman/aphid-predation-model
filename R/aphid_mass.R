# Aphid mass at age ---------------------------------------------------------
#
# No aphid weights were recorded, so mass at age is a simple assumed model:
# exponential growth from neonate mass at birth to adult mass at maturity,
# then constant. Anchors:
#   * adult apterous mass: pea 3.8 mg, bean 0.9 mg (dissertation Ch. 1 & 2)
#   * pea neonate 0.15 mg (literature value for 1st instars; source not yet
#     verified - see docs/model-design.md)
#   * bean neonate: no value found; assume the same neonate:adult ratio as pea
#   * age at maturity: median age at first reproduction in the clip cages
# Treat these as provisional and revisit if better data turn up.

aphid_mass_anchors <- function(clip_cage = read_clip_cage()) {
  maturity <- clip_cage |>
    dplyr::filter(offspring > 0) |>
    dplyr::summarise(first_repro = min(day), .by = c(aphid, id)) |>
    dplyr::summarise(age_mature = median(first_repro) - 0.5, .by = aphid)

  tibble::tibble(
    aphid = factor(c("pea", "bean"), levels = c("pea", "bean")),
    mass_adult = c(3.8, 0.9),
    mass_neonate = c(0.15, 0.9 * 0.15 / 3.8)
  ) |>
    dplyr::left_join(maturity, by = "aphid") |>
    dplyr::mutate(growth_rate = log(mass_adult / mass_neonate) / age_mature)
}

#' Mass (mg) of aphids of a given age (days); vectorized over `age`.
aphid_mass <- function(age, anchors) {
  pmin(anchors$mass_neonate * exp(anchors$growth_rate * age), anchors$mass_adult)
}

#' Age (days) at which an aphid reaches `mass` (mg); inverse of aphid_mass().
aphid_age_at_mass <- function(mass, anchors) {
  log(mass / anchors$mass_neonate) / anchors$growth_rate
}
