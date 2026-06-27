# Open-Source Sea Turtle Population Viability Analysis & Trend Toolkit

This repository provides an open-access, stochastic simulation framework for sea turtle population assessments. The underlying statistical methods transition traditional deterministic evaluations into state-space modeling environments, directly supporting robust Endangered Species Act (ESA) Section 7 consultation frameworks and conservation benefit accounting.

The methodologies embedded here are based on the state-space architectures developed by **Martin et al. (2020)** and the Conway-Maxwell-Poisson (CMP) take distribution adaptations defined by **Siders et al. (2023)**. This version has been expanded to support a **generalized multi-mode inversion**, allowing users to assess either traditional anthropogenic losses ("take") or population recovery impacts ("benefits") with strict mathematical sign parity.

---

## Prerequisites and Installation

### 1. System Dependency: JAGS
The trend estimation engine relies on **JAGS (Just Another Gibbs Sampler)** via the `jagsUI` package. You must have the core JAGS binary application installed on your operating system before running the R code:
* [Download and Install JAGS](https://sourceforge.net/projects/mcmc-jags/)

### 2. R Packages
Execute the following block in your R console to install the necessary processing, simulation, multivariate distribution drawing, and visualization libraries:

```R
install.packages(c("dplyr", "tidyr", "ggplot2", "mvtnorm", "truncnorm", "jagsUI"))
```

---

## Repository Directory

### Core Model Files
* **`singleUQ.txt`**: A multivariate state-space trend model written automatically to your workspace by the engine. It estimates a single, underlying, regionally shared trend ($U$) across an arbitrary number of nesting beaches while isolating process variance ($Q$) and site-specific observation errors ($R_j$).

### Execution Scripts
* **`sea_turtle_pva.R`**: The primary, self-contained pipeline script. It features a multi-mode projection loop that dynamically scales filters to any timeline or number of beaches. The script automates a JAGS state-space trend model, executes a von Bertalanffy growth inversion to track lifecycle structures, and implements a bivariate size-mortality framework to map individual reproductive value metrics ($ANE$) into forward-looking policy options.

---

## Output Routing Ledger

All generated plots and data metrics are systematically dumped into a local `/output` folder automatically generated in your working directory:

| Filename | Type | Description |
| :--- | :--- | :--- |
| `table_1_abundance_matrix.csv` | CSV Data | Cleaned historical nesting matrix formatted for JAGS. |
| `table_2_jags_trend_scorecard.csv` | CSV Data | Posterior summaries for baseline parameters ($U$ and $Q$). |
| `table_3_threshold_collapse_risk.csv` | CSV Data | Calculated probability of the stock dropping below $50\%$, $25\%$, or $12.5\%$ thresholds. |
| `table_4_fishery_ledger.csv` | CSV Data | Comparative metrics contrasting observed vs. counterfactual histories. |
| `figure_1_jags_model_fit.png` | Image Plot | Multivariate state-space model fits over the historical data. |
| `figure_2_pva_forward_abundance_trajectories.png` | Image Plot | Forward forecasts assuming completely unimpacted natural trends. |
| `figure_3_pva_scenario_overlay.png` | Image Plot | Adaptive scenario comparison displaying true history vs. management trajectories. |
| `figure_4_isolated_demographic_difference.png` | Image Plot | The isolated, standalone impact profile ($\Delta$ Total Population) across the timeline. |

---

## How to Adapt This Toolkit for Your Own Data

This framework can be applied to any marine megafauna population with nesting time series and known interaction/addition histories. To adapt the pipeline to your target dataset, configure these three control parameters at the top of `sea_turtle_pva.R`:

### 1. Set the Objective Mode & Species Target
Choose whether you want to calculate fishery takes or conservation additions, and name your target species run identifier:

```R
target_species <- "dc"   # Your target species/population tag code

# --- GLOBAL MODEL MODE SWITCH ---
# Set to "take" to model fishery removals (Original Gold Standard)
# Set to "benefit" to model conservation additions (The Inverse Framework)
model_mode     <- "benefit" 
```

### 2. File Path and Data Ingestion Configuration
Switch the user execution toggle to read your custom files instead of the automatic mock workspace generator:

```R
use_mock_data <- FALSE
```

Once disabled, uncomment the tracking parameters directly below it and point them to your data sheets:

```R
nest_file  <- "path/to/your/nest_counts_time_series.csv"
safe_file  <- "path/to/your/total_annual_expanded_take_or_benefit.csv"
obs_file   <- "path/to/your/individual_observer_measurements.csv"
```
* *Data Sheet Formatting Constraint:* Your `nest_file` must contain a distinct column tracking the calendar census year (e.g., `Year` or `Season_Year`). All remaining columns will be automatically parsed as separate, individual index beaches—protecting the trend engine against single-site dimension collapses.

### 3. Species-Specific Biological Constants
Calibrate the life-history modifiers within the script's biological parameter block to align with your target species' published traits:

```R
sp_clutch_freq <- 5.5      # Average nests deposited per female per season
sp_remig_int   <- 3.06     # Average years between breeding seasons
sp_pf          <- 0.73     # Sex ratio modifier (proportion female)

# von Bertalanffy growth function parameters for size-to-age inversion
sp_linf        <- 142.7    # Asymptotic maximum straight carapace length (cm)
sp_k           <- 0.2262   # Curvature growth speed coefficient
sp_tknot       <- -0.17    # Theoretical age at zero length
```

---

## Citations
* **Martin et al. 2020** (NOAA Tech Memo NMFS-PIFSC-101): *Evaluation of the population-level impacts of incidental take in the Hawaii deep-set longline fishery on sea turtles.*
* **Siders et al. 2023** (Update to Tech Memo NMFS-PIFSC-101): *Applying the Conway-Maxwell-Poisson distribution to model under- and over-dispersed protected species interactions.*
