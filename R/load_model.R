# Load the full simulation model ---------------------------------------------
#
# Sources every R/ file the simulation scripts (analysis/04-06) need: data
# readers, fitting and sampling helpers, parameters, the simulation itself,
# scenarios, sensitivity helpers, and plot conventions. Run from the project
# root.

for (f in c("R/data.R", "R/aphid_demography.R", "R/aphid_events.R", "R/aphid_mass.R",
            "R/predator_behavior.R", "R/predator_development.R", "R/params.R",
            "R/simulation.R", "R/scenarios.R", "R/sensitivity.R", "R/plot_theme.R")) {
  source(f)
}
rm(f)
