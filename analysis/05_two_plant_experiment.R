# Two-plant experiment ----------------------------------------------------
#
# One plant species (two fava plants), two aphid species (pea, bean), one
# predator (a single Hippodamia convergens larva). Aphid density dependence:
# births thinned by max(0, 1 - density / capacity), density in adult-mass
# equivalents per plant (both species combined).
#
# 1. Aphid dynamics without predators: plateaus and competition on a shared
#    plant.
# 2. Indirect effects of bean aphids on pea aphids via the predator: the
#    larva starts as an L1 on plant 1, which always holds pea aphids; what is
#    on plant 2, or alongside the pea aphids on plant 1, varies. Effects are
#    measured against the same scenario without the predator.
#
# Free parameters (search, capture, travel, giving-up time, starvation,
# learning rate, aphid capacity) are at provisional defaults, so results are
# illustrative until the sensitivity analysis is done.
#
# Requires output/params from analysis/01-04.
# Outputs: output/figures/exp_*.png, output/params/two_plant_experiment.rds

library(dplyr)
library(tidyr)
library(ggplot2)

for (f in c("R/data.R", "R/aphid_demography.R", "R/aphid_events.R", "R/aphid_mass.R",
            "R/predator_behavior.R", "R/predator_development.R", "R/params.R",
            "R/simulation.R", "R/scenarios.R", "R/plot_theme.R")) source(f)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
n_cores <- max(1L, parallel::detectCores() - 1L)
RNGkind("L'Ecuyer-CMRG")
set.seed(4040)

params <- default_params()
stopifnot(is.finite(params$bean_handling_experienced))

run_many <- function(scenario, reps, summarize) {
  parallel::mclapply(seq_len(reps), \(i) summarize(simulate_two_plants(params, scenario)),
                     mc.cores = n_cores, mc.set.seed = TRUE)
}

census_of <- function(runs) bind_rows(purrr::map(runs, "census"), .id = "rep")

# 1. Aphid dynamics without predators ----------------------------------------

no_pred <- tibble(
  setting = c("Pea aphids alone", "Bean aphids alone", "Pea + bean aphids together"),
  plant1 = list(c(pea = 5), c(bean = 5), c(pea = 5, bean = 5))
) |>
  mutate(census = purrr::map(plant1, \(p1) {
    sc <- two_plant_scenario(p1, c(pea = 0), predator = FALSE, run_days = 60)
    census_of(run_many(sc, 24, \(r) list(census = r$census)))
  })) |>
  select(setting, census) |>
  unnest(census) |>
  filter(plant == 1, n > 0 | day == 0)

capacity_lines <- tibble(species = c("pea", "bean"),
                         n = params$aphid_capacity / params$mass$mass_adult[match(c("pea", "bean"), params$mass$aphid)])

p_no_pred <- no_pred |>
  summarise(median = median(n), lo = quantile(n, 0.1), hi = quantile(n, 0.9), .by = c(setting, species, day)) |>
  mutate(setting = factor(setting, levels = c("Pea aphids alone", "Bean aphids alone", "Pea + bean aphids together"))) |>
  filter(!(setting == "Pea aphids alone" & species == "bean"), !(setting == "Bean aphids alone" & species == "pea")) |>
  ggplot(aes(day, median, colour = species, fill = species)) +
  geom_hline(data = capacity_lines, aes(yintercept = n, colour = species), linetype = "22", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~setting, nrow = 1) +
  scale_colour_manual(values = species_colours, labels = species_labels) +
  scale_fill_manual(values = species_colours, labels = species_labels) +
  labs(
    x = "Day", y = "Aphids on the plant", colour = NULL, fill = NULL,
    title = "Aphid dynamics on one fava plant without predators",
    subtitle = sprintf("Capacity %s mg adult-mass equivalents (dashed: capacity in aphids of one species). Median, 10th–90th percentiles of 24 runs",
                       format(params$aphid_capacity, big.mark = ","))
  ) +
  theme_model()
ggsave("output/figures/exp_aphids_no_predator.png", p_no_pred, width = 10, height = 4, dpi = 200)

# 2. Indirect effects via the predator -------------------------------------

run_days <- 40
reps <- 100

scenarios <- tibble(
  scenario = c("Pea | none", "Pea | pea", "Pea | bean", "Pea + bean | none"),
  plant1 = list(c(pea = 5), c(pea = 5), c(pea = 5), c(pea = 5, bean = 5)),
  plant2 = list(c(pea = 0), c(pea = 5), c(bean = 5), c(pea = 0))
)

# Simulations take ~15 min on 11 cores; cached so summaries and figures can
# be revised without rerunning. Set rerun <- TRUE after changing the model.
rerun <- FALSE
cache <- "output/params/two_plant_experiment_runs.rds"
if (rerun || !file.exists(cache)) {
  experiment <- scenarios |>
    crossing(predator = c(TRUE, FALSE)) |>
    mutate(runs = purrr::pmap(list(plant1, plant2, predator), \(p1, p2, pred) {
      sc <- two_plant_scenario(p1, p2, predator = pred, run_days = run_days)
      run_many(sc, reps, \(r) list(
        census = r$census,
        fate = r$predator$fate,
        fate_day = r$predator$fate_day,
        meals = r$meals |> count(plant, species)
      ))
    }))
  saveRDS(experiment, cache)
} else {
  experiment <- readRDS(cache)
}

census <- experiment |>
  mutate(census = purrr::map(runs, census_of)) |>
  select(scenario, predator, census) |>
  unnest(census) |>
  mutate(scenario = factor(scenario, levels = scenarios$scenario))

p_focal <- census |>
  filter(plant == 1, species == "pea") |>
  summarise(median = median(n), lo = quantile(n, 0.25), hi = quantile(n, 0.75),
            .by = c(scenario, predator, day)) |>
  ggplot(aes(day, median, linetype = predator)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = predator), fill = species_colours[["pea"]], alpha = 0.15) +
  geom_line(colour = species_colours[["pea"]], linewidth = 0.7) +
  facet_wrap(~scenario, nrow = 1) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"),
                        labels = c(`TRUE` = "With predator", `FALSE` = "No predator")) +
  labs(
    x = "Day", y = "Pea aphids on plant 1", linetype = NULL,
    title = "Pea aphids on the predator's starting plant",
    subtitle = sprintf("Two fava plants. Panels: aphid species on plant 1 | plant 2 (5 adults each). Median and IQR of %d runs", reps)
  ) +
  theme_model()
ggsave("output/figures/exp_focal_pea.png", p_focal, width = 10, height = 4, dpi = 200)

# Suppression: 1 - (pea aphid-days on plant 1 with predator) / (mean without),
# per run with a predator, over the whole run.
aphid_days <- census |>
  filter(plant == 1, species == "pea") |>
  summarise(aphid_days = sum(n), .by = c(scenario, predator, rep))

suppression <- aphid_days |>
  filter(predator) |>
  left_join(
    aphid_days |> filter(!predator) |> summarise(baseline = mean(aphid_days), .by = scenario),
    by = "scenario"
  ) |>
  mutate(suppression = 1 - aphid_days / baseline)

predator_outcomes <- experiment |>
  filter(predator) |>
  mutate(out = purrr::map(runs, \(rs) tibble(
    rep = as.character(seq_along(rs)),
    fate = purrr::map_chr(rs, "fate"),
    fate_day = purrr::map_dbl(rs, "fate_day"),
    pea_eaten = purrr::map_dbl(rs, \(r) sum(r$meals$n[r$meals$species == "pea"])),
    bean_eaten = purrr::map_dbl(rs, \(r) sum(r$meals$n[r$meals$species == "bean"]))
  ))) |>
  select(scenario, out) |>
  unnest(out) |>
  mutate(scenario = factor(scenario, levels = scenarios$scenario))

summary_table <- suppression |>
  left_join(predator_outcomes, by = c("scenario", "rep")) |>
  summarise(
    suppression_se = sd(suppression) / sqrt(n()),
    suppression = mean(suppression),
    pupated = mean(fate == "pupated"),
    died = mean(fate %in% c("died", "starved", "failed_pupation")),
    pupation_day = median(fate_day[fate == "pupated"]),
    pea_eaten = mean(pea_eaten),
    bean_eaten = mean(bean_eaten),
    .by = scenario
  ) |>
  arrange(scenario)
summary_table

p_suppression <- suppression |>
  ggplot(aes(scenario, suppression)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_violin(fill = species_colours[["pea"]], colour = NA, alpha = 0.2) +
  stat_summary(fun.data = \(x) mean_se(x, mult = 1.96), colour = species_colours[["pea"]], size = 0.4) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    x = "Aphid species on plant 1 | plant 2", y = "Reduction in pea aphid-days on plant 1",
    title = "How much one larva suppresses pea aphids on its starting plant",
    subtitle = sprintf("Over %d days, relative to the same scenario without a predator. Point: mean ± 95%% CI; shape: distribution across runs", run_days)
  ) +
  theme_model()
ggsave("output/figures/exp_suppression.png", p_suppression, width = 7, height = 4.5, dpi = 200)

saveRDS(list(census = census, suppression = suppression, predator = predator_outcomes,
             summary = summary_table, params = params),
        "output/params/two_plant_experiment.rds")
