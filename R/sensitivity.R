# Sensitivity analysis helpers ----------------------------------------------
#
# The design and elementary effects come from sensitivity::morris(); these
# helpers define the free parameters' ranges and map design points onto the
# simulation's parameter list. Factors with `log = TRUE` are sampled on a
# log10 scale: their design columns hold log10(value).

free_parameter_ranges <- function() {
  tibble::tribble(
    ~factor,              ~low,   ~high,  ~log,
    "encounter_rate",     0.001,  0.02,   TRUE,  # per aphid per minute
    "capture_pea",        0.1,    0.6,    FALSE,
    "capture_bean",       0.5,    1.0,    FALSE,
    "digestion_rate",     2,      8,      TRUE,  # per day
    "learn_meals",        5,      100,    TRUE,
    "giving_up_time",     30,     480,    TRUE,  # minutes
    "travel_time",        15,     240,    TRUE,  # minutes
    "starvation_scale",   0.5,    2,      TRUE,  # multiplier on starvation medians
    "aphid_capacity",     1000,   4000,   TRUE   # mg adult-mass equivalents
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

#' Apply factor values (natural units) to a parameter list.
apply_factors <- function(params, values) {
  v <- as.list(values)
  if (!is.null(v$encounter_rate)) { params$search_rate <- v$encounter_rate; params$plant_area <- 1 }
  if (!is.null(v$capture_pea)) params$capture[["pea"]] <- v$capture_pea
  if (!is.null(v$capture_bean)) params$capture[["bean"]] <- v$capture_bean
  if (!is.null(v$digestion_rate)) params$digestion_rate <- v$digestion_rate
  if (!is.null(v$learn_meals)) params$learn_meals <- v$learn_meals
  if (!is.null(v$giving_up_time)) params$giving_up_time <- v$giving_up_time
  if (!is.null(v$travel_time)) params$travel_time <- v$travel_time
  if (!is.null(v$starvation_scale)) params$starvation_median <- params$starvation_median * v$starvation_scale
  if (!is.null(v$aphid_capacity)) params$aphid_capacity <- v$aphid_capacity
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
