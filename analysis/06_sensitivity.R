# Sensitivity analysis of free parameters (Morris screening) ---------------
#
# Which free (unmeasured) parameters drive the predator's suppression of pea
# aphids and the indirect effects of bean aphids on pea aphids?
#
# Scenarios (two fava plants; aphids on plant 1 | plant 2, 2 founding adults
# each; one L1 larva hatches on plant 1 on `hatch_day` (a factor); 35-day
# runs):
#   pp  = Pea | pea
#   pb  = Pea | bean
#   mix = Pea + bean | none
# Outputs per parameter set (means over replicate runs):
#   supp_*    1 - pea aphid-days on plant 1 with predator / without
#   ie_separate = supp_pb - supp_pp   (bean vs pea aphids on the other plant)
#   ie_shared   = supp_mix - supp_pp  (bean aphids on the same plant; includes
#                                      competition between aphid species)
#   pupate_*  proportion of larvae pupating
#
# Morris screening with sensitivity::morris(): 15 factors, 4 levels, grid jump
# 2, 10 trajectories = 160 parameter sets. Log-range factors are sampled on a
# log10 scale; elementary effects are scaled by factor range (scale = TRUE).
# 20 replicate runs per scenario per set, with common random numbers (same
# seeds in every set). Results are cached per set in output/sensitivity/, so
# the script can be stopped and resumed. ~40 min on 11 cores.
#
# Outputs: output/sensitivity/*.rds, output/figures/sens_*.png

library(dplyr)
library(tidyr)
library(ggplot2)
library(sensitivity)

for (f in c("R/data.R", "R/aphid_demography.R", "R/aphid_events.R", "R/aphid_mass.R",
            "R/predator_behavior.R", "R/predator_development.R", "R/params.R",
            "R/simulation.R", "R/scenarios.R", "R/sensitivity.R", "R/plot_theme.R")) source(f)

# Environment overrides (for quick smoke tests): SENS_OUT, SENS_R, SENS_REPS
out_dir <- Sys.getenv("SENS_OUT", "output/sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
n_cores <- max(1L, parallel::detectCores() - 1L)

base_params <- default_params()
stopifnot(is.finite(base_params$bean_handling_experienced))
ranges <- free_parameter_ranges()

run_days <- 35
reps <- as.integer(Sys.getenv("SENS_REPS", "20"))
n_trajectories <- as.integer(Sys.getenv("SENS_R", "10"))
baseline_reps <- 2L * reps
make_scenarios <- function(hatch_day = 3) list(
  pp = two_plant_scenario(c(pea = 2), c(pea = 2), run_days = run_days, predator_day = hatch_day),
  pb = two_plant_scenario(c(pea = 2), c(bean = 2), run_days = run_days, predator_day = hatch_day),
  mix = two_plant_scenario(c(pea = 2, bean = 2), c(pea = 0), run_days = run_days, predator_day = hatch_day)
)
scenarios <- make_scenarios()
no_predator <- purrr::map(scenarios, \(s) { s$predator <- NULL; s })

focal_pea_days <- function(run) {
  sum(run$census$n[run$census$plant == 1 & run$census$species == "pea"])
}

# Design ----------------------------------------------------------------------

design_file <- file.path(out_dir, "morris_design.rds")
if (!file.exists(design_file)) {
  set.seed(1991)
  bounds <- design_bounds(ranges)
  mo <- morris(
    model = NULL, factors = ranges$factor, r = n_trajectories,
    design = list(type = "oat", levels = 4, grid.jump = 2),
    binf = bounds$binf, bsup = bounds$bsup, scale = TRUE
  )
  saveRDS(mo, design_file)
}
mo <- readRDS(design_file)
design <- as_tibble(mo$X) |>
  mutate(across(everything(), \(x) round(x, 8)), set = row_number(), .before = 1)

# No-predator baselines: depend only on aphid capacity and density lag ------

baseline_file <- file.path(out_dir, "baselines.rds")
if (!file.exists(baseline_file)) {
  tasks <- design |>
    distinct(aphid_capacity, density_lag) |>
    crossing(scenario = names(scenarios), rep = seq_len(baseline_reps))
  pea_days <- parallel::mclapply(seq_len(nrow(tasks)), \(i) {
    tk <- tasks[i, ]
    p <- apply_factors(base_params, c(aphid_capacity = 10^tk$aphid_capacity, density_lag = 10^tk$density_lag))
    focal_pea_days(simulate_two_plants(p, no_predator[[tk$scenario]], seed = 5e5 + i))
  }, mc.cores = n_cores)
  baselines <- tasks |>
    mutate(pea_days = unlist(pea_days)) |>
    summarise(baseline = mean(pea_days), .by = c(aphid_capacity, density_lag, scenario))
  saveRDS(baselines, baseline_file)
}
baselines <- readRDS(baseline_file)

# Parameter sets --------------------------------------------------------------

run_set <- function(row) {
  values <- design_values(row[ranges$factor], ranges)
  p <- apply_factors(base_params, values)
  scenarios <- make_scenarios(values[["hatch_day"]])
  tasks <- crossing(scenario = names(scenarios), rep = seq_len(reps))
  res <- parallel::mclapply(seq_len(nrow(tasks)), \(i) {
    tk <- tasks[i, ]
    # common random numbers: the seed depends only on scenario and replicate
    seed <- 1e4 * match(tk$scenario, names(scenarios)) + tk$rep
    r <- simulate_two_plants(p, scenarios[[tk$scenario]], seed = seed)
    tibble(pea_days = focal_pea_days(r), fate = r$predator$fate, fate_day = r$predator$fate_day)
  }, mc.cores = n_cores)
  bind_cols(tasks, bind_rows(res))
}

for (i in seq_len(nrow(design))) {
  f <- file.path(out_dir, sprintf("set_%03d.rds", i))
  if (file.exists(f)) next
  t0 <- Sys.time()
  saveRDS(run_set(design[i, ]), f)
  message(sprintf("set %d/%d done in %.0f s", i, nrow(design), as.numeric(Sys.time() - t0, units = "secs")))
}

# Outputs per set ---------------------------------------------------------------

set_outputs <- design |>
  mutate(runs = purrr::map(set, \(i) readRDS(file.path(out_dir, sprintf("set_%03d.rds", i))))) |>
  select(set, aphid_capacity, density_lag, runs) |>
  unnest(runs) |>
  left_join(baselines, by = c("aphid_capacity", "density_lag", "scenario")) |>
  summarise(
    supp = 1 - mean(pea_days) / first(baseline),
    pupate = mean(fate == "pupated"),
    .by = c(set, scenario)
  ) |>
  pivot_wider(names_from = scenario, values_from = c(supp, pupate)) |>
  mutate(ie_separate = supp_pb - supp_pp, ie_shared = supp_mix - supp_pp) |>
  arrange(set)

outputs <- c("supp_pp", "supp_pb", "supp_mix", "ie_separate", "ie_shared",
             "pupate_pp", "pupate_pb", "pupate_mix")
stopifnot(!anyNA(set_outputs[outputs]))
tell(mo, as.matrix(set_outputs[outputs]))
effects <- morris_table(mo)
results <- design |> left_join(set_outputs, by = "set")
saveRDS(list(morris = mo, results = results, effects = effects), file.path(out_dir, "morris.rds"))

effects |>
  filter(output %in% c("supp_pp", "ie_separate", "ie_shared")) |>
  arrange(output, desc(mu_star)) |>
  print(n = 30)

# Figures ---------------------------------------------------------------------

output_labels <- c(
  supp_pp = "Suppression: Pea | pea", supp_pb = "Suppression: Pea | bean",
  supp_mix = "Suppression: Pea + bean | none",
  ie_separate = "Indirect effect: bean on other plant", ie_shared = "Indirect effect: bean on same plant",
  pupate_pp = "Pupation: Pea | pea", pupate_pb = "Pupation: Pea | bean", pupate_mix = "Pupation: Pea + bean | none"
)
factor_labels <- c(
  detection_width = "Detection width", plant_area = "Plant area",
  colony_min_area = "Colony minimum area", colony_density = "Colony packing density",
  prey_size_ratio = "Prey size limit", capture_pea = "Capture: pea", capture_bean = "Capture: bean",
  leave_scaling = "Leave-time scaling", plant_path = "Path between plants",
  digestion_rate = "Digestion rate", learn_meals = "Bean learning (meals)",
  starvation_scale = "Starvation tolerance", aphid_capacity = "Aphid capacity",
  density_lag = "Density-dependence lag", hatch_day = "Hatch day"
)

p_rank <- effects |>
  filter(output %in% c("supp_pp", "ie_separate", "ie_shared", "pupate_pb")) |>
  mutate(
    output = factor(output, levels = c("supp_pp", "ie_separate", "ie_shared", "pupate_pb"), labels = output_labels[c("supp_pp", "ie_separate", "ie_shared", "pupate_pb")]),
    # order factors by mu_star within each facet
    factor = reorder(paste(factor_labels[factor], output, sep = "___"), mu_star)
  ) |>
  ggplot(aes(mu_star, factor)) +
  geom_segment(aes(x = 0, xend = mu_star, yend = factor), colour = "grey75", linewidth = 2.5) +
  geom_point(aes(x = sigma), shape = 4, size = 2, colour = "grey30") +
  facet_wrap(~output, scales = "free", nrow = 1, labeller = label_wrap_gen(28)) +
  scale_y_discrete(labels = \(x) sub("___.*$", "", x)) +
  labs(
    x = "Influence (bar: μ*, mean |elementary effect|; ×: σ, interactions/nonlinearity)", y = NULL,
    title = "Which free parameters matter?",
    subtitle = sprintf("Morris screening (sensitivity::morris): %d trajectories × %d factors, %d runs per scenario per set. Effects scaled by factor range",
                       n_trajectories, nrow(ranges), reps)
  ) +
  theme_model()
fig_dir <- if (out_dir == "output/sensitivity") "output/figures" else out_dir
ggsave(file.path(fig_dir, "sens_morris_ranking.png"), p_rank, width = 13, height = 4.5, dpi = 200)

p_spread <- results |>
  select(set, all_of(c("supp_pp", "ie_separate", "ie_shared", "pupate_pp", "pupate_pb", "pupate_mix"))) |>
  distinct() |>
  pivot_longer(-set, names_to = "output") |>
  mutate(output = factor(output, levels = names(output_labels), labels = output_labels)) |>
  ggplot(aes(value, output)) +
  geom_vline(xintercept = 0, colour = "grey75") +
  geom_point(position = position_jitter(height = 0.15, width = 0, seed = 1), alpha = 0.5, size = 1.2, colour = "#2a78d6") +
  labs(
    x = "Value across the 100 parameter sets", y = NULL,
    title = "Range of outcomes across free-parameter space",
    subtitle = "Indirect effect < 0: bean aphids reduce how much the larva suppresses pea aphids on plant 1"
  ) +
  theme_model()
ggsave(file.path(fig_dir, "sens_output_spread.png"), p_spread, width = 8, height = 4.5, dpi = 200)
