# Two-plant simulation: calibration, validation, first experiment ---------
#
# 1. Calibrates how fast bean-experienced larvae handle bean aphids so that
#    simulated bean-diet vials match the 2008 kill rate of 4th instars.
# 2. Checks the full simulation (vial mode) against the diet experiment.
# 3. Runs a first two-plant experiment on indirect effects of bean aphids on
#    pea aphids. Free parameters (search, capture, travel, giving-up time,
#    starvation, learning rate) are at provisional defaults, so these results
#    are illustrative, not conclusions.
#
# Outputs:
#   output/params/simulation.rds
#   output/figures/sim_*.png

library(dplyr)
library(tidyr)
library(ggplot2)

for (f in c("R/data.R", "R/aphid_demography.R", "R/aphid_events.R", "R/aphid_mass.R",
            "R/predator_behavior.R", "R/predator_development.R", "R/params.R",
            "R/simulation.R", "R/scenarios.R", "R/plot_theme.R")) source(f)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
n_cores <- max(1L, parallel::detectCores() - 1L)
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)

run_many <- function(params, scenario, reps, summarize = identity) {
  parallel::mclapply(seq_len(reps), \(i) summarize(simulate_two_plants(params, scenario)),
                     mc.cores = n_cores, mc.set.seed = TRUE)
}

params <- default_params()
diet_intervals <- read_diet_experiment()
diet_summary <- read_mac_csv(data_path("diet_summary.csv")) |> janitor::clean_names()

diet_labels <- c(P = "Pea diet", M = "Mixed diet", B = "Bean diet")
diet_colours <- c(P = species_colours[["pea"]], M = "#8c6bb1", B = species_colours[["bean"]])

vial_runs <- function(diet, m, reps) {
  p <- vial_params(params)
  p$bean_handling_experienced <- m
  run_many(p, vial_scenario(diet, p), reps, \(r) summarize_vial_run(r, p)) |>
    bind_rows(.id = "rep") |>
    mutate(diet = diet, m = m)
}

# 1. Calibrate experienced bean handling -------------------------------------

obs_rates <- fit_max_intake(diet_intervals) # kills/day by diet and instar
target_bean_l4 <- obs_rates |> filter(diet == "B", stage == "L4") |> pull(kills_per_day)

calib <- tibble(m = seq(0.3, 0.8, by = 0.1)) |>
  mutate(runs = purrr::map(m, \(m) vial_runs("B", m, reps = 200))) |>
  mutate(l4_rate = purrr::map_dbl(runs, \(r) {
    r |> filter(stage == "L4", !is.na(days), days > 0) |> summarise(mean(kills / days)) |> pull()
  }))
calib |> select(m, l4_rate)

# kill rate falls with m; interpolate the m that gives the observed rate
bean_handling_experienced <- approx(calib$l4_rate, calib$m, xout = target_bean_l4)$y
bean_handling_experienced
params$bean_handling_experienced <- bean_handling_experienced
saveRDS(list(bean_handling_experienced = bean_handling_experienced, calibration = select(calib, m, l4_rate)),
        "output/params/simulation.rds")

p_calib <- ggplot(calib, aes(m, l4_rate)) +
  geom_hline(yintercept = target_bean_l4, colour = "grey60", linetype = "22") +
  geom_line(colour = species_colours[["bean"]]) +
  geom_point(colour = species_colours[["bean"]], size = 2) +
  annotate("text", x = max(calib$m), y = target_bean_l4, label = "Observed (2008 vials)",
           hjust = 1, vjust = -0.6, size = 3, colour = "grey40") +
  labs(
    x = "Bean handling multiplier for experienced larvae (m)",
    y = "L4 bean kills per day",
    title = "Calibrating faster bean handling with experience",
    subtitle = sprintf("Simulated bean-diet vials; m = %.2f matches the observed rate", bean_handling_experienced)
  ) +
  theme_model()
ggsave("output/figures/sim_bean_handling_calibration.png", p_calib, width = 6, height = 4, dpi = 200)

# 2. Validate against the diet experiment ------------------------------------

validation <- c("P", "M", "B") |>
  purrr::map(\(d) vial_runs(d, bean_handling_experienced, reps = 400)) |>
  bind_rows()

pupation_check <- validation |>
  distinct(diet, rep, fate) |>
  summarise(sim = mean(fate == "pupated"), .by = diet) |>
  left_join(diet_summary |> summarise(obs = mean(pupate), n = n(), .by = diet), by = "diet")
pupation_check

obs_stage <- diet_intervals |>
  filter(stage %in% larval_stages) |>
  summarise(days = n(), kills = sum(pea + bean), .by = c(id, diet, stage)) |>
  left_join(select(diet_summary, id, pupate), by = "id") |>
  filter(pupate == 1) |>
  mutate(
    # observed L4 includes the pre-pupa; simulated L4 'days' end at the pre-pupa
    days = if_else(stage == "L4", days - params$prepupa_days, as.numeric(days)),
    source = "Observed"
  )

sim_stage <- validation |>
  group_by(diet, rep) |>
  filter(any(fate == "pupated")) |>
  ungroup() |>
  mutate(source = "Simulated")

stage_check <- bind_rows(
  select(obs_stage, diet, stage, days, kills, source),
  select(sim_stage, diet, stage, days, kills, source)
) |>
  filter(!is.na(days), days > 0, !is.na(kills)) |>
  mutate(rate = kills / days) |>
  pivot_longer(c(days, rate), names_to = "measure", values_to = "value") |>
  summarise(mean = mean(value), lo = quantile(value, 0.1), hi = quantile(value, 0.9),
            .by = c(diet, stage, source, measure)) |>
  mutate(
    diet = factor(diet, levels = c("P", "M", "B")),
    measure = factor(measure, levels = c("days", "rate"),
                     labels = c("Days in instar (feeding)", "Kills per day"))
  )

p_validation <- ggplot(stage_check, aes(stage, mean, colour = diet, shape = source)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.35) +
  facet_grid(measure ~ diet, scales = "free_y", switch = "y",
             labeller = labeller(diet = diet_labels)) +
  scale_colour_manual(values = diet_colours, guide = "none") +
  scale_shape_manual(values = c(Observed = 16, Simulated = 1)) +
  labs(
    x = NULL, y = NULL, shape = NULL,
    title = "Full simulation in vial mode vs the 2008 diet experiment",
    subtitle = paste0(
      "Larvae that pupated; points = mean, bars = 10th–90th percentiles. Pupation (obs / sim): ",
      paste(sprintf("%s %.0f%% / %.0f%%", c(P = "pea", M = "mixed", B = "bean")[pupation_check$diet],
                    100 * pupation_check$obs, 100 * pupation_check$sim), collapse = "; ")
    )
  ) +
  theme_model() +
  theme(strip.placement = "outside")
ggsave("output/figures/sim_vial_validation.png", p_validation, width = 9, height = 5.5, dpi = 200)

# 3. First two-plant experiment ----------------------------------------------
#
# Two fava bean plants; aphid species differ only in where they start. The
# predator (one L1 Hippodamia convergens) starts on plant 1, which always
# holds pea aphids. Does what is on plant 2 - or alongside the pea aphids on
# plant 1 - change how much the predator suppresses pea aphids on plant 1?
# Effects are measured against the same scenario without the predator.
# Scenario labels: aphids on plant 1 | aphids on plant 2.

scenarios <- tibble(
  scenario = c("Pea | none", "Pea | pea", "Pea | bean", "Pea + bean | none"),
  plant1 = list(c(pea = 5), c(pea = 5), c(pea = 5), c(pea = 5, bean = 5)),
  plant2 = list(c(pea = 0), c(pea = 5), c(bean = 5), c(pea = 0))
)

reps <- 150
experiment <- scenarios |>
  crossing(predator = c(TRUE, FALSE)) |>
  mutate(runs = purrr::pmap(list(plant1, plant2, predator), \(p1, p2, pred) {
    sc <- two_plant_scenario(p1, p2, predator = pred, run_days = 14)
    run_many(params, sc, reps, \(r) list(
      census = r$census,
      fate = r$predator$fate,
      fate_day = r$predator$fate_day,
      meals = r$meals |> count(plant, species)
    ))
  }))

census <- experiment |>
  mutate(census = purrr::map(runs, \(rs) bind_rows(purrr::map(rs, "census"), .id = "rep"))) |>
  select(scenario, predator, census) |>
  unnest(census)

focal <- census |>
  filter(plant == 1, species == "pea") |>
  summarise(median = median(n), lo = quantile(n, 0.25), hi = quantile(n, 0.75),
            .by = c(scenario, predator, day)) |>
  mutate(scenario = factor(scenario, levels = scenarios$scenario))

p_focal <- ggplot(focal, aes(day, median, linetype = predator)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = predator), fill = species_colours[["pea"]], alpha = 0.15) +
  geom_line(colour = species_colours[["pea"]], linewidth = 0.7) +
  facet_wrap(~scenario, nrow = 1) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"),
                        labels = c(`TRUE` = "With predator", `FALSE` = "No predator")) +
  labs(
    x = "Day", y = "Pea aphids on plant 1 (log scale)", linetype = NULL,
    title = "Pea aphids on the predator's starting plant",
    subtitle = "Two fava plants. Panels: aphid species on plant 1 | plant 2 (5 adults each). Median and IQR of 150 runs"
  ) +
  theme_model()
ggsave("output/figures/sim_focal_pea.png", p_focal, width = 10, height = 4, dpi = 200)

# Suppression of pea on plant 1 at day 14 relative to no predator
suppression <- census |>
  filter(plant == 1, species == "pea", day == 14) |>
  summarise(mean_n = mean(n), .by = c(scenario, predator)) |>
  pivot_wider(names_from = predator, values_from = mean_n, names_prefix = "pred_") |>
  mutate(suppression = 1 - pred_TRUE / pred_FALSE)

predator_outcomes <- experiment |>
  filter(predator) |>
  mutate(
    fates = purrr::map(runs, \(rs) tibble(fate = purrr::map_chr(rs, "fate"))),
    meals = purrr::map(runs, \(rs) bind_rows(purrr::map(rs, "meals"), .id = "rep"))
  )

fates <- predator_outcomes |>
  select(scenario, fates) |>
  unnest(fates) |>
  count(scenario, fate) |>
  mutate(prop = n / sum(n), .by = scenario)

diet <- predator_outcomes |>
  select(scenario, meals) |>
  unnest(meals) |>
  summarise(n = sum(n) / reps, .by = c(scenario, species)) |>
  pivot_wider(names_from = species, values_from = n, values_fill = 0)

experiment_summary <- suppression |>
  left_join(diet, by = "scenario") |>
  left_join(fates |> filter(fate == "pupated") |> select(scenario, pupated = prop), by = "scenario")
experiment_summary

saveRDS(list(census = census, fates = fates, summary = experiment_summary),
        "output/params/two_plant_experiment.rds")
