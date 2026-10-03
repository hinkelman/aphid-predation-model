# Aphid–ladybird predation model

An individual-based, discrete-event simulation of one convergent ladybird
larva (*Hippodamia convergens*) foraging on two fava bean plants (*Vicia faba*)
colonised by pea aphids (*Acyrthosiphon pisum*) and black bean aphids
(*Aphis fabae*). The model is parameterised from laboratory experiments run in
2005–2010 (see `data/`) and asks how the two aphid species affect each other
through their shared predator.

## Executive summary

### Motivation

Pea aphids and bean aphids differ sharply as prey for *H. convergens* larvae:

- **Pea aphids are high-quality prey that are hard to catch.** They escape by
  walking away or dropping off the plant.
- **Bean aphids are poor, possibly toxic, prey that are easy to catch.**
  - Only 13% of larvae reared on bean aphids reached adulthood, against 70% on
    pea aphids.
  - After eating a single bean aphid, a larva handles prey for much longer,
    moves less, and stays on the plant longer.

The question is how these differences shape **predator-mediated indirect
effects** of bean aphids on pea aphids. Does the presence of bean aphids raise
pea aphids' predation risk (shared doom) or lower it (associational refuge)?
Several mechanisms could act:

- **Short term (behaviour):** long bean meals, lethargy afterwards, and staying
  longer on plants after bean meals all change how the predator's time and
  attention are divided between the prey.
- **Longer term (predator performance):** slower development, higher mortality
  and failed pupation on bean aphids change how much predation the larva
  delivers over its life.

### Methods

**Data.** Six experiments, re-analysed from the original data and the
author's dissertation (Ch. 1–2):

- clip-cage life tables of both aphid species on fava;
- larval development and daily kills on bean, pea and mixed diets;
- diet switches and starvation in the 4th instar;
- video-tracked single-prey trials: handling time, leaving the leaf, and
  activity after a meal;
- handling time by aphid age.

All predator experiments used size-matched prey and pea-reared larvae.

**Aphids** (`R/aphid_*.R`):

- Lifespans follow a Weibull distribution; at each event, the remaining
  lifespan is drawn given the aphid's current age.
- Births follow the age-specific birth rate from the clip cages, with a
  per-aphid fecundity multiplier for individual variation.
- Crowding (adapted from the ALMaSS aphid model, Thomsen et al. 2024) raises
  aphid mortality, with a 1-day lag, and triggers the production of winged
  aphids that leave the plant. Crowding is measured per gram of plant, and
  aphids of the other species count with weight α.
- Colonies start from 10 founders whose ages follow the species' stable age
  distribution.
- Aphid mass at age is a provisional growth curve.

**Predator behaviour** (`R/predator_behavior.R`), from the 2008 trials:

- Handling time depends on aphid species × hunger, and scales with aphid age.
- Bean handling gets faster with experience, calibrated to the kill rates in
  the 2008 diet-experiment vials.
- Leaving a plant: time to leave a leaf depends on aphid species × hunger.
  Leaving a whole plant means giving up on many leaves in turn, so the time is
  summed over leaves.
- Activity drops after a bean meal (lethargy).
- Search rate = walking speed × activity × detection width, scaled by body
  size.
- Larvae search inside colonies, whose area depends on aphid numbers and on
  species-specific packing density.
- Capture depends on aphid species and on prey size relative to larval size.
- Intake is limited by a digesting gut.

**Predator development and survival** (`R/predator_development.R`), fitted to
each larva's own observed daily kills:

- **Food thresholds per instar** are 9, 17, 31 and 100 units, where one unit
  is one size-matched pea aphid.
- **A bean kill is worth 0.50 units on its own, but only 0.15 when the larva
  is also eating pea.**
- **Mortality:**
  - a background death rate for each stage;
  - in the first instar, a bean-related death risk that pea in the diet
    reduces;
  - a chance of failing to pupate that rises with lifetime bean intake.
- Starved 4th instars that have eaten more than a critical amount still
  pupate.

**Simulation** (`R/simulation.R`):

- Continuous time, in minutes. The single larva and the aphid populations run
  together as one event-driven simulation.
- Default setup: two medium fava plants (800 cm² searchable, about 20 g), 0.9 m
  apart. Each colony is founded by 10 aphids, and the larva hatches on day 3.
  Runs last 35 days.

**Validation and analysis:**

- Simulated larvae on ad-lib diets in vials reproduce the observed pupation
  rates (pea 70% vs 77%, mixed 53% vs 52%, bean 13% vs 13%). They also match
  instar durations and kill rates. Simulated aphid cohorts reproduce the
  clip-cage survival and fecundity curves.
- The free parameters were screened with a Morris sensitivity analysis
  (`sensitivity::morris`): 17 parameters, 180 parameter sets, 3 scenarios, and
  20 replicate runs each with shared random seeds.

### Key results

All four scenarios start with pea aphids on plant 1, where the larva hatches.
The **predator's impact on pea aphids** is the reduction in pea aphid-days on
plant 1 relative to the same scenario without a predator. The **benefit to
pea aphids from bean aphids** compares that impact with the Pea | pea
scenario, where the alternative prey is another pea colony; positive means
bean aphids reduce pea aphids' predation. Results are illustrative: several
behavioural parameters have no direct data.

| Aphids on plant 1 \| plant 2 | Predator impact on pea aphids | Benefit to pea from bean aphids | Larvae pupating | Days to pupation |
|---|---|---|---|---|
| Pea \| none | 33% | — | 60% | 13.7 |
| Pea \| pea | 42% | (reference) | 71% | 14.0 |
| Pea \| bean | 32% | +0.11 | 65% | 14.3 |
| Pea + bean \| none | −6% (pea aphids *gain*) | +0.48 | 54% | 17.8 |

(Experiment `analysis/05`, 100 runs per scenario. Sensitivity analysis
`analysis/06`, 180 parameter sets.)

1. **Bean aphids on the larva's plant act as an associational refuge for pea
   aphids.** They reduce pea aphids' predation in 98% of parameter sets (median
   benefit +0.16; 10th–90th percentile +0.03 to +0.55). Three mechanisms
   contribute:
   - **Dilution and distraction (short term).** The larva spends search time,
     handling time and gut capacity on easily caught bean aphids. Long,
     lethargic bean meals take it out of action, and pea aphids eaten fall
     from 192 to 161 per larva.
   - **Reduced predator performance (longer term).** On a mixed diet the larva
     develops more slowly (17.8 vs 14.0 days to pupation), spends longer in the
     small, low-intake instars, and more larvae die or fail to pupate.
   - **Predator-mediated competitive release.** The two aphids partly compete
     for the plant. By eating bean aphids, the larva frees resources for pea
     aphids, which can end up more numerous with the predator than without
     it. The sensitivity analysis suggests this is a minor contributor: the
     strength of competition between the aphids (α) had little influence.
   - **Shared doom loses to the refuge.** Bean meals make larvae stay on a
     plant longer, which keeps the predator near pea aphids. The model
     includes this retention, but dilution and reduced predator performance
     outweigh it.
2. **Bean aphids on a neighbouring plant give at most a weak refuge.** Their
   benefit to pea aphids is positive in 74% of parameter sets but small
   (median +0.02). A larva seldom leaves a pea colony it is still exploiting,
   so it rarely experiences both prey. In this model the spatial scale of the
   indirect effect is set by predator movement, and it is mostly within a
   plant.
3. **The refuge is strongest when the predator is good at catching pea
   aphids.** Pea capture success is the most influential parameter for every
   output: it sets how much predation there is for bean aphids to divert.
   Bean-specific parameters (capture, colony packing, learning rate) barely
   matter.
4. **Timing matters.** One larva strongly depresses a pea colony only if it
   arrives while the colony is small; hatch day is the second-most
   influential parameter.
5. **Consuming bean aphids carries the costs seen in the lab:**
   - bean food is worth about half as much as pea, and less when pea is also
     eaten;
   - mortality concentrates in the first instar;
   - pupation fails more often.

   These costs drive the longer-term part of the refuge, through reduced
   predator performance.
6. **The 2005 and 2008 handling-time experiments disagreed** on how long larvae
   took to handle bean aphids (about 15–20 vs 70–130 min). The difference is
   mostly explained by starvation: 2008 larvae were starved for 2–24 h,
   against about 2–4 h in 2005, and bean handling lengthens with hunger.

### Main uncertainties

- **Pea capture success per encounter**, the parameter that most strongly
  scales the refuge. The default is 0.1. The only
  quantitative evidence comes from adult beetles, which consumed 1 of 72 pea
  aphids encountered in field arenas (Nelson & Rosenheim 2006).
- **Colony packing densities.** There are no published values; the defaults
  are bounded by body size.
- **Leaf-to-plant scaling of leaving times.** This is an assumption.
- **Aphid crowding parameters** (competition α, plant biomass, lag). ALMaSS
  values are field-scale.
- **Model scope:**
  - one larva within one generation: there is no predator reproduction or
    aggregation, so the model captures effects through the predator's
    behaviour and individual performance, not changes in predator numbers;
  - no cannibalism or intraguild predation;
  - no non-lethal effects of disturbance on pea aphids (dropping costs
    feeding time);
  - no handling-time data for instars other than the 4th;
  - a provisional mass-at-age curve.

Details, decisions and sources are in [docs/model-design.md](docs/model-design.md).

## Running the analyses

Requires R ≥ 4.4 with dplyr, tidyr, ggplot2, purrr, readxl, janitor, survival,
MASS and sensitivity. Run the scripts with `Rscript` from the project root.
Scripts 04–06 fork with `parallel::mclapply()`, which IDE consoles such as
Positron may block.

```bash
Rscript analysis/01_aphid_demography.R       # aphid lifespan, births, frailty
Rscript analysis/02_predator_behavior.R      # handling, leaving, movement, aphid mass
Rscript analysis/03_predator_development.R   # thresholds, mortality, pupation (~1 min)
Rscript analysis/04_simulation_calibration.R # experienced bean handling; vial validation
Rscript analysis/05_two_plant_experiment.R   # two-plant experiment (~5 min)
Rscript analysis/06_sensitivity.R            # Morris sensitivity analysis (~1 h)
```

Fitted parameters, cached runs and figures are written to `output/` (not
tracked). Scripts 05 and 06 reuse cached simulations. Delete
`output/params/two_plant_experiment_runs.rds` or `output/sensitivity/` to
rerun them after changing the model.

## Repository layout

```
R/          data import, model fits, samplers, simulation engine, scenarios, sensitivity helpers
analysis/   numbered scripts: fit, calibrate, validate, experiment, sensitivity
data/       raw data from the 2005–2010 experiments (read-only)
docs/       model design document: data interpretation, decisions, sources
```
