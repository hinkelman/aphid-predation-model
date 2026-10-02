# Two-plant simulation: calibration and validation -------------------------
#
# 1. Calibrates how fast bean-experienced larvae handle bean aphids so that
#    simulated bean-diet vials match the 2008 kill rate of 4th instars.
# 2. Checks the full simulation (vial mode) against the diet experiment.
# The two-plant experiments are in analysis/05.
#
# Outputs:
#   output/params/simulation.rds
#   output/figures/sim_bean_handling_calibration.png, sim_vial_validation.png

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
