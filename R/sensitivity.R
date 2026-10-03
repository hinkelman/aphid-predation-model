# Sensitivity analysis helpers ----------------------------------------------
#
# The design and elementary effects come from sensitivity::morris(); these
# helpers define the free parameters' ranges and map design points onto the
# simulation's parameter list. Factors with `log = TRUE` are sampled on a
# log10 scale: their design columns hold log10(value).

free_parameter_ranges <- function() {
  tibble::tribble(
    ~factor,              ~low,   ~high,  ~log,
    "detection_width",    2,      5,      FALSE, # mm, L4
    "plant_area",         400,    1600,   TRUE,  # cm^2 searchable surface
    "colony_min_area",    2,      40,     TRUE,  # cm^2
    "colony_density",     1,      10,     TRUE,  # mg aphid per cm^2 of colony
    "prey_size_ratio",    1,      4,      TRUE,  # prey vs larval length scale
    "capture_pea",        0.02,   0.3,    TRUE,  # ~1/72 for adults on alfalfa
    "capture_bean",       0.5,    1.0,    FALSE,
    "leave_scaling",      0.25,   4,      TRUE,  # multiplier on area-scaled leave times
    "plant_path",         450,    1800,   TRUE,  # mm between plants
    "digestion_rate",     2,      8,      TRUE,  # per day
    "learn_meals",        5,      100,    TRUE,
    "starvation_scale",   0.5,    2,      TRUE,  # multiplier on starvation medians
    "plant_biomass",      10,     40,     TRUE,  # g fresh; sets crowding density
    "competition_alpha",  0,      1,      FALSE, # weight of other-species aphids
    "density_lag",        0,      3,      FALSE, # days, rounded to whole days
    "hatch_day",          2,      5,      FALSE  # day the larva hatches (scenario)
  )
}

#' Lower and upper bounds for sensitivity::morris() (log10 for log factors).
design_bounds <- function(ranges = free_parameter_ranges()) {
  list(
    binf = ifelse(ranges$log, log10(ranges$low), ranges$low),
    bsup = ifelse(ranges$log, log10(ranges$high), ranges$high)
  )
}

#' Convert one row of the design matrix to parameter values (natural units).
design_values <- function(x, ranges = free_parameter_ranges()) {
  x <- unlist(x)[ranges$factor]
  setNames(ifelse(ranges$log, 10^x, x), ranges$factor)
}

#' Apply factor values (natural units) to a parameter list. `hatch_day` is a
#' scenario setting, applied by the caller.
apply_factors <- function(params, values) {
  v <- as.list(values)
  direct <- c("detection_width", "plant_area", "colony_min_area", "colony_density",
              "prey_size_ratio", "leave_scaling", "plant_path", "digestion_rate",
              "learn_meals", "plant_biomass", "competition_alpha")
  for (nm in intersect(names(v), direct)) params[[nm]] <- v[[nm]]
  if (!is.null(v$density_lag)) params$density_lag <- as.integer(round(v$density_lag))
  if (!is.null(v$capture_pea)) params$capture[["pea"]] <- v$capture_pea
  if (!is.null(v$capture_bean)) params$capture[["bean"]] <- v$capture_bean
  if (!is.null(v$starvation_scale)) params$starvation_median <- params$starvation_median * v$starvation_scale
  params
}

#' Tidy Morris statistics from a told sensitivity::morris object with a
#' multi-output response (x$ee is trajectories x factors x outputs).
morris_table <- function(mo) {
  ee <- mo$ee
  if (length(dim(ee)) == 2) ee <- array(ee, c(dim(ee), 1), dimnames = c(dimnames(ee), list("y")))
  stat <- function(f) {
    m <- apply(ee, c(2, 3), f)
    tibble::as_tibble(m, rownames = "factor") |>
      tidyr::pivot_longer(-factor, names_to = "output")
  }
  stat(mean) |>
    dplyr::rename(mu = value) |>
    dplyr::left_join(dplyr::rename(stat(\(e) mean(abs(e))), mu_star = value), by = c("factor", "output")) |>
    dplyr::left_join(dplyr::rename(stat(sd), sigma = value), by = c("factor", "output"))
}
