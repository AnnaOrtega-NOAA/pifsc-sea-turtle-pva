# =======================================================================================
# Sea turtle population viability analysis validation engine (~North Pacific loggerheads)
# Timeline envelope: 1985 - 2015 (31-year baseline matrix)
# Core authors: Summer L. Martin, Tomoharu Eguchi, and Zachary Siders
# Collaborators and Reviewers: Robert Ahrens, Todd Jones, Milani Chaloupka, Larry Crowder
# Adapted by: Anna Ortega (CIMAR)
# Description: Provides a monolithic, self-contained analytical pipeline that integrates 
#              state-space trend estimation with stochastic population viability analysis 
#              (PVA) to evaluate fishery-induced demographic risks. The engine 
#              implements a multivariate state-space model in JAGS to determine a shared 
#              regional growth rate across nesting beaches, applies a 
#              Conway-Maxwell-Poisson distribution to generate over-dispersed interaction 
#              frequencies, and utilizes a von Bertalanffy growth function 
#              inversion to map turtle lengths to lifecycle age trajectories. Individual 
#              reproductive value losses are tracked via annual nester equivalents (ANE) 
#              across a 100-year forward projection horizon, enabling a comparative 
#              counterfactual assessment of take versus no-take population scenarios. 
#              Features an integrated control switch to seamlessly toggle between automated 
#              mock data generation and custom user-supplied datasets.
# =======================================================================================

# Load necessary libraries for data cleaning, restructuring, graphics, and distributions
library(dplyr)      # Used for data manipulation, cleaning, and filtering pipelines
library(tidyr)      # Used for restructuring matrices (pivoting from wide to long formats)
library(ggplot2)    # Used for creating publication-grade predictive trajectories and graphics
library(mvtnorm)    # Used for simulating multivariate normal distributions for length and mortality
library(truncnorm)  # Used for handling truncated normal boundaries in demographic simulations
library(jagsUI)     # Used as the interface to execute the JAGS MCMC state-space engine

# ==============================================================================
# 0. User execution controls and directory setup
# ==============================================================================

# User configuration toggle: Set to TRUE to run a test pipeline using fake data.
# Set to FALSE when you are ready to point the script toward your own local CSV files.
use_mock_data <- TRUE 

# Define the target species identifier code used across generated files
target_species <- "cc"

# Dynamically establish a relative output directory pathway within the current working folder
out_dir <- paste0(getwd(), "/output")

# Check if the output folder exists on the machine; if not, generate it recursively
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Output initialization success messages directly to the R console
message("Initialization successful for ", target_species, " (1985-2015). Storing in: ", out_dir)

# ==============================================================================
# 1. Data routing switch (Mock generation vs. user file configuration)
# ==============================================================================

if (use_mock_data) {
  
  # --- Step A: Define local file names for the automated mock workspace ---
  nest_file  <- "combined_loggerhead_nests.csv" # Saved name for mock nesting counts
  safe_file  <- "cc_safe_rep.csv"               # Saved name for mock total take estimates
  obs_file   <- "cc_dsll_interactions.csv"      # Saved name for mock observer measurements
  
  message("Alert: Generating complete mock environment (nesting, safe, and observer data)...")
  set.seed(2026) # Lock the random number generator seed for exact table stability
  
  # Establish the temporal bounds of the 31-year historical timeline matrix
  baseline_years <- 1985:2015
  n_years_baseline <- length(baseline_years)
  
  # 1. Mock nesting data matrix (Simulating fluctuating annual counts for three index beaches)
  mock_nest_data <- data.frame(
    Season_Year = baseline_years,
    Maehama = round(pmax(100, rnorm(n_years_baseline, mean = 800, sd = 250) + (baseline_years - 1985) * 10)),
    Inakahama = round(pmax(80, rnorm(n_years_baseline, mean = 600, sd = 180) + (baseline_years - 1985) * 5)),
    Yotsusehama = round(pmax(20, rnorm(n_years_baseline, mean = 200, sd = 70)))
  )
  write.csv(mock_nest_data, nest_file, row.names = FALSE) # Export the matrix to the local folder
  
  # 2. Mock safe reporting data (Simulating estimated total expanding interaction numbers)
  mock_safe_data <- data.frame(
    Year = baseline_years,
    Species = "Loggerhead",
    CC = round(runif(n_years_baseline, min = 10, max = 35))
  )
  write.csv(mock_safe_data, safe_file, row.names = FALSE) # Export total estimates to local folder
  
  # 3. Mock observer interaction records (Simulating measured sample sizes for length and mortality)
  mock_sample_size <- 65
  mock_obs_years <- sample(2002:2015, mock_sample_size, replace = TRUE) # Sample years within active limits
  
  mock_obs_data <- data.frame(
    Year = sort(mock_obs_years),
    Species = "Loggerhead",
    SCL_raw = round(rnorm(mock_sample_size, mean = 69.0, sd = 7.5), 1), # Mean length from published text
    M_mu = round(runif(mock_sample_size, min = 0.15, max = 0.85), 3)     # Simulated individual capture mortality rates
  )
  write.csv(mock_obs_data, obs_file, row.names = FALSE) # Export individual observations to local folder
  
  message("Mock data factory execution complete. All local CSV dependencies written to workspace.")
  
} else {
  
  # --- Step B: User data input configuration ---
  # To run your own custom datasets, change 'use_mock_data' to FALSE and uncomment the lines below.
  # Modify the file paths inside the quotation marks to match your local CSV files exactly.
  
  # nest_file  <- "path/to/your/nest_counts_time_series.csv"
  # safe_file  <- "path/to/your/total_annual_expanded_take.csv"
  # obs_file   <- "path/to/your/individual_observer_measurements.csv"
  
  # Note on file column structures required for downstream code:
  # 1. nest_file requirements: Must contain a column representing the year and columns for each beach.
  # 2. safe_file requirements: Must contain a 'Year' column and a column specifying total model take.
  # 3. obs_file requirements: Must contain a 'Year' column, a carapace length column, and a mortality rate column.
  
}

# ==============================================================================
# 2. Species biology profiles and constants
# ==============================================================================

# Species baseline parameters configured for North Pacific Loggerheads
sp_clutch_freq <- 4.6      # Average number of nests laid by a female in a single nesting year
sp_remig_int   <- 3.3      # Average number of years between active nesting seasons for an adult female
sp_ane_pj      <- 0.80     # Baseline annual survival probability for a juvenile sea turtle
sp_ane_pa      <- 0.895    # Baseline annual survival probability for a mature adult sea turtle
sp_linf        <- 86.9     # von Bertalanffy growth model asymptotic maximum length in centimeters
sp_k           <- 0.09     # von Bertalanffy growth coefficient determining curvature speed
sp_tknot       <- -2.467   # Theoretical age parameter at which size is exactly zero
sp_lmat        <- 86.9 * 0.99 # Straight Carapace Length proxy triggering mature ogive maturation
sp_sig_mat     <- 6.34     # Standard deviation scale parameter determining maturity ogive logistic slope
sp_max_age     <- 60.0     # Hard biological upper boundary ceiling for survival lifespan
sp_pf          <- 0.65     # Population sex ratio scale factor representing the proportion of female hatchlings

# Conway-Maxwell-Poisson (CMP) parameters for stochastic future take generation
cmp_mu         <- 15.0     # Central expected value anchor for future longline fishery interaction attempts
cmp_nu         <- 0.1      # Dispersion coefficient (values < 1 define heavy over-dispersion models)
beach_cols     <- c("Maehama", "Inakahama", "Yotsusehama") # Names of the indexing nesting columns to track

n_sims  <- 1000   # Total count of Monte Carlo iterations used to build error matrices
horizon <- 100    # Number of sequential years to project forward into the future PVA timeline

# ==============================================================================
# 3. Mathematical conversion functions
# ==============================================================================

# Standard inverse-logit mapping function to scale real values back onto a 0-1 probability scale
inv_logit <- function(x) { 
  exp(x) / (1 + exp(x)) 
}

# Safe-bounded logit transformation protecting extreme boundaries from infinite results
safe_logit <- function(x) { 
  x_safe <- pmax(0.001, pmin(0.999, x)) # Constrain values strictly between 0.001 and 0.999
  log(x_safe / (1 - x_safe))            # Execute standard log-odds math transformation
}

# Recursive algorithm computing the normalizing constant Z for the Conway-Maxwell-Poisson layout
compute_CMP_constant <- function(Lambda, Nu, Mu, Tol, Max, Log=TRUE, Type="Z"){
  # Direct large-sample asymptotic calculation approximation to save memory and processing time
  if( (!is.na(Lambda) & Lambda > 10^Nu) | (!is.na(Mu) & Mu^Nu > 10^Nu) ){
    if(Type=="Z"){ ln_Const = Nu*Lambda^(1/Nu) - ((Nu-1)/(2*Nu))*log(Lambda) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu) }
    if(Type=="S"){ ln_Const = Nu*Mu - ((Nu-1)/(2))*log(Mu) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu) }
  } else {
    # Sequential series expansion looping until probability changes drop below user tolerance threshold
    Const = rep(0,Max+1); Index = 1; Const[Index] = 1
    while( Const[Index]/Const[1] > Tol ){
      if(Type=="Z") Const[Index+1] = Const[Index] * ( Lambda / Index^Nu )
      if(Type=="S") Const[Index+1] = Const[Index] * ( Mu / Index )^Nu; Index = Index + 1
    }
    ln_Const = log(sum(Const)) # Return log-scale sum of expansion series
  }
  if(Log) return(ln_Const) else return(exp(ln_Const))
}

# Custom density function computing log-likelihoods under a Conway-Maxwell-Poisson model
dCMP <- function( x, lambda, mu, nu, log=TRUE, tol=0.01, iter.max=200 ){
  if(missing(mu) & !missing(lambda)) loglike = x*log(lambda) - nu*lfactorial(x) - compute_CMP_constant(Lambda=lambda, Nu=nu, Mu=NA, Tol=tol, Max=iter.max, Log=TRUE, Type="Z")
  if(!missing(mu) & missing(lambda)) loglike = nu*x*log(mu) - nu*lfactorial(x) - compute_CMP_constant(Lambda=NA, Nu=nu, Mu=mu, Tol=tol, Max=iter.max, Log=TRUE, Type="S")
  if(log) return(loglike) else return(exp(loglike))
}

# Stochastic random sampler generating discrete integer counts under CMP dispersion
rCMP <- function( n, lambda, mu, nu, tol=0.01, x_max=200 ){
  loglike_x = rep(NA, x_max+1) # Pre-allocate vector for sample range density probabilities
  for( x in 0:x_max ){
    # Populate the distribution arrays with calculated point probabilities
    if(missing(mu) & !missing(lambda)) loglike_x[x+1] = dCMP( x=x, lambda=lambda, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
    if(!missing(mu) & missing(lambda)) loglike_x[x+1] = dCMP( x=x, mu=mu, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
  }
  # Draw a random observation using the built-in probability weight distribution array
  return(sample(x=0:x_max, size=n, replace=TRUE, prob=exp(loglike_x)))
}

# ==============================================================================
# 4. Self-contained JAGS model generator
# ==============================================================================
# Writes the multivariate state-space model structure to disk as a text file for compiler ingestion
cat("model{
   A[1] <- 0; # Set the baseline offset for primary indexing beach to zero anchor
   for(j in 2:n.timeseries) { A[j] ~ dnorm(a_mean, 1/(a_sd^2)); } # Estimate beach offset deviations
   U ~ dnorm(u_mean, 1/(u_sd^2)); # Shared log-scale regional population trend parameter
   tauQ ~ dgamma(q_alpha, q_beta); # Gamma prior on process error precision
   Q <- 1/tauQ; # Convert precision parameter into process error variance
   X[1] ~ dnorm(x0_mean, 1/(x0_sd^2)); # Initial baseline state tracking node assignment
   for(v in 1:1){ X0[v] <- X[1] - U; } # Back-calculate initial index state tracking point
   for(tt in 2:n.yrs) {
      predX[tt] <- X[tt-1] + U; # Project latent state expectation forward using trend parameter
      X[tt] ~ dnorm(predX[tt], tauQ); # Latent population state vector subject to process error
   }
   for(j in 1:n.timeseries) {
      tauR[j] ~ dgamma(r_alpha, r_beta); # Site precision observation error distributions
      R[j] <- 1/tauR[j]; # Convert observation precision indices into actual variance scales
      for(tt in 1:n.yrs) {
          predY[j,tt] <- Z[j,1] * X[tt] + A[j]; # Reconstruct specific beach log-scale expectation tracks
          Y[j,tt] ~ dnorm(predY[j,tt], tauR[j]); # Evaluate true observations against observation error
      }
   }
}", file = "singleUQ.txt")

# ==============================================================================
# Part 1 and 2: Data ingestion and matrix shaping (1985 - 2015 bounds)
# ==============================================================================
message("--- Phase 1 and 2: Building trend matrices for ", target_species, " (1985-2015) ---")

# Import the local working data files from the project folder paths
df_raw  <- read.csv(nest_file, stringsAsFactors = FALSE)
obs_df  <- read.csv(obs_file, stringsAsFactors = FALSE)
safe_df <- read.csv(safe_file, stringsAsFactors = FALSE)

# Dynamically locate the specific text column tracking the calendar census year
year_col <- names(df_raw)[grep("year", names(df_raw), ignore.case = TRUE)[1]]

# Filter, scale, and clean the census data matrix for JAGS engine consumption
abund_df <- df_raw %>%
  rename(Year = !!sym(year_col)) %>% # Relabel the dynamic year column header name to a standard label
  filter(Year >= 1985 & Year <= 2015) %>% # Restrain dataset strictly to 31 baseline years
  select(Year, all_of(beach_cols)) %>% # Drop extra metadata columns outside tracking beach arrays
  group_by(Year) %>% # Collapse any duplicate row entries within shared census years
  summarise(across(all_of(beach_cols), ~sum(as.numeric(.), na.rm=TRUE)), .groups="drop") %>%
  tidyr::pivot_longer(cols = all_of(beach_cols), names_to = "Site", values_to = "Annual_Nesters") %>% # Pivot long
  mutate(Annual_Nesters = Annual_Nesters / sp_clutch_freq) %>% # Convert raw nest counts into active nesting females
  mutate(Annual_Nesters = ifelse(Annual_Nesters <= 0, NA_real_, Annual_Nesters)) # Replace zero values with true structural NAs

# Save the structured historical matrix cleanly into your relative output folder path
write.csv(abund_df, paste0(out_dir, "/table_1_abundance_matrix.csv"), row.names = FALSE)

# ==============================================================================
# Trend engine execution
# ==============================================================================
message("--- Running original singleUQ trend model ---")

# Define execution function wrapping the JAGS interface loop
run_trend_model <- function(df_abund, beaches) {
  # Reconstruct structural arrays into wide dimensions for row indexing transposition
  df_wide <- df_abund %>% select(Year, Site, Annual_Nesters) %>% pivot_wider(names_from=Site, values_from=Annual_Nesters) %>% arrange(Year)
  d_mat <- t(log(df_wide[, beaches])) # Transpose log-transformed counts into a site-by-year data matrix
  
  n_ts <- length(beaches) # Calculate total tracking beach columns
  Z_m <- matrix(0, n_ts+1, n_ts) # Pre-allocate indicator routing allocation mapping array
  for(i in 1:n_ts) Z_m[i,1] <- 1 # Route all site pathways back to single latent population tracking trend
  
  # Build the formal data list structure explicitly demanded by the compiled text script
  j_data <- list(
    Y = rbind(d_mat, NA), n.yrs = 31, n.timeseries = n_ts, Z = Z_m,
    a_mean=0, a_sd=4, u_mean=0, u_sd=0.5,
    q_alpha=0.01, q_beta=0.01, r_alpha=0.01, r_beta=0.01,
    x0_mean=mean(d_mat, na.rm=TRUE), x0_sd=10
  )
  # Execute MCMC chains sequentially inside the native engine shell
  jagsUI::jags(
    data = j_data, inits = NULL,
    parameters.to.save = c("A", "U", "Q", "R", "X0", "X"),
    model.file = "singleUQ.txt", n.chains = 3, n.burnin = 25000, n.thin = 10, n.iter = 50000, parallel = FALSE, verbose = TRUE
  )
}

# Run the unadjusted baseline execution track
trend_base <- run_trend_model(abund_df, beach_cols)

# Extract latent population state tracks from posterior vectors
X_post <- trend_base$sims.list$X
total_pop_post <- exp(X_post) # Initialize total population array scaled back from log formats
for(j in 2:length(beach_cols)) {
  # Add individual beach counts back into pool by scaling deviations via regional offset posteriors
  total_pop_post = total_pop_post + exp(X_post + trend_base$sims.list$A[, j])
}
total_pop_post <- total_pop_post * sp_remig_int # Scale active nesting population back into absolute adult cohort population size
start_take_draws <- total_pop_post[, 31] # Extract all draws matching terminal baseline evaluation index (Year 2015)

# ------------------------------------------------------------------------------
# Generate baseline figure 1 layout graphics
# ------------------------------------------------------------------------------
# Extract summary statistics to recreate classic base layout visual distribution curves
thedata <- abund_df %>% select(Year, Site, Annual_Nesters) %>% pivot_wider(names_from=Site, values_from=Annual_Nesters) %>% arrange(Year)
raw_sums <- rowSums(thedata[, beach_cols], na.rm = TRUE); raw_sums[raw_sums == 0] <- NA; log_obs <- log(raw_sums)

X_total <- apply(cbind(trend_base$sims.list$X0, trend_base$sims.list$X), 2, function(v) rowSums(apply(trend_base$sims.list$A, 2, function(x) exp(v + x))))
X0_total <- rowSums(apply(trend_base$sims.list$A, 2, function(x) exp(trend_base$sims.list$X0 + x)))
X_total_med <- apply(log(X_total), 2, median); X_q <- apply(log(X_total), 2, quantile, probs = c(0.025, 0.5, 0.975))

yrange_vals <- c(log_obs, X_q); yrange <- range(yrange_vals[is.finite(yrange_vals)], na.rm = TRUE)
den.X0 <- density(log(X0_total), adj = 2); den.X0$y2 <- (den.X0$y / max(den.X0$y)); q.X0 <- quantile(log(X0_total), probs = c(0.025, 0.975)); xid <- sapply(q.X0, function(x) { which.min(abs(x - den.X0$x)) })
den.N0 <- density(log(X_total[, ncol(X_total)]), adj = 2); den.N0$y2 <- (den.N0$y / max(den.N0$y)); q.N0 <- quantile(log(X_total[, ncol(X_total)]), probs = c(0.025, 0.975)); xid.n0 <- sapply(q.N0, function(x) { which.min(abs(x - den.N0$x)) })

# Open drawing device canvas pathway and construct graphic elements
png(filename = paste0(out_dir, "/figure_1_jags_model_fit.png"), width = 7, height = 4.62, units = "in", res = 300)
layout(matrix(1:2, 1, 2), width = c(1, 0.23))
par(mar = c(4, 4, 1, 1))

plot(thedata$Year, log_obs, pch = 16, type = "n", las = 1, ylab = "log(Annual Nesters)", xlab = "Season", ylim = yrange, xlim = c(min(thedata$Year) - 1.1, max(thedata$Year) + 0.25), xaxs = 'i')
polygon(c(c(min(thedata$Year) - 1, thedata$Year), rev(c(min(thedata$Year) - 1, thedata$Year))), c(X_q[1, ], rev(X_q[3, ])), col = 'grey85', border = NA)
lines(c(min(thedata$Year) - 1, thedata$Year), X_q[2, ], lwd = 3, col = 'gray50')
polygon(x = c(rep((thedata$Year[1] - 1), length(xid[1]:xid[2])), rev(den.X0$y2[xid[1]:xid[2]] + (thedata$Year[1] - 1))), y = c(den.X0$x[xid[1]:xid[2]], rev(den.X0$x[xid[1]:xid[2]])), lwd = 2, border = NA, col = adjustcolor("dodgerblue3", alpha.f = 0.3))
lines(den.X0$y2[xid[1]:xid[2]] + (thedata$Year[1] - 1), den.X0$x[xid[1]:xid[2]], lwd = 2, col = 'dodgerblue3')
polygon(x = c(rep(thedata$Year[nrow(thedata)], length(xid.n0[1]:xid.n0[2])), rev(-den.N0$y2[xid.n0[1]:xid.n0[2]] + thedata$Year[nrow(thedata)])), y = c(den.N0$x[xid.n0[1]:xid.n0[2]], rev(den.N0$x[xid.n0[1]:xid.n0[2]])), lwd = 2, border = NA, col = adjustcolor("darkorchid3", alpha.f = 0.3))
lines(-den.N0$y2[xid.n0[1]:xid.n0[2]] + thedata$Year[nrow(thedata)], den.N0$x[xid.n0[1]:xid.n0[2]], lwd = 2, col = 'darkorchid3')
points(thedata$Year, log_obs, pch = 16)
points(c(min(thedata$Year) - 1, thedata$Year), X_total_med, pch = 16, col = c('dodgerblue3', rep('red', length(X_total_med) - 2), "darkorchid3"))

par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = c(expression(sum(N[list(obs,j)], j, "")), expression(sum(T[j] + a[j], j, "")), "Median r", "95% r", expression(paste(T[0], " (95%CI)")), expression(paste(N[final], " (95% CI)"))), pch = c(16, 16, NA, 15, NA, NA), lwd = c(NA, NA, 3, NA, 2, 2), pt.cex = c(1, 1, NA, 3, NA, NA), col = c("black", "red", "gray50", "gray85", "dodgerblue3", "darkorchid3"), bty = "n", y.intersp = 0.9, xpd = NA)
dev.off()

# Save primary unadjusted baseline summary performance markers into a clean scorecard file
table_jags_df <- data.frame(Framework = "JAGS Baseline", Strategy = "Regional Shared Trend", LogScale_U = mean(trend_base$sims.list$U), GrowthRate_Pct = (exp(mean(trend_base$sims.list$U))-1)*100, Process_Variance_Q = mean(trend_base$sims.list$Q))
write.csv(table_jags_df, paste0(out_dir, "/table_2_jags_trend_scorecard.csv"), row.names = FALSE)

# ==============================================================================
# Part 3: Headless PVA forward forecasts
# ==============================================================================
message("--- Phase 3: Projecting unconstrained PVA matrix ---")

# Establish future projection calendar framework timeline series indices
proj_years <- (max(abund_df$Year) + 1):(max(abund_df$Year) + horizon)
sim_reg_matrix <- matrix(NA, nrow = n_sims, ncol = horizon) # Pre-allocate tracking array

# Draw random matching distribution pairings directly out of baseline posterior distribution arrays
r_reg <- sample(trend_base$sims.list$U, n_sims, replace = TRUE)
q_reg <- sample(trend_base$sims.list$Q, n_sims, replace = TRUE)
n0_reg <- sample(start_take_draws, n_sims, replace = TRUE)

# Execute forward iterative population viability projections
pb_pva <- txtProgressBar(min = 0, max = n_sims, style = 3) # Open interactive loop loading bar
for(s in 1:n_sims) {
  current <- n0_reg[s] # Extract initial base population size draw allocation
  for(t in 1:horizon) {
    # Step forward annually subject to log-growth expectation rates and process error variance
    current <- rnorm(1, current * exp(r_reg[s]), sqrt(q_reg[s]))
    sim_reg_matrix[s,t] <- current # Save calculation into trajectory matrix row positions
  }
  if(any(sim_reg_matrix[s,] < 0)) {
    # If values cross below zero point, lock all remaining years to absolute zero extinction state
    extinction_yr <- min(which(sim_reg_matrix[s,] < 0))
    sim_reg_matrix[s, extinction_yr:horizon] <- 0
  }
  setTxtProgressBar(pb_pva, s) # Step progress bar indicators forward
}
close(pb_pva) # Close tracking graphics tool

# Summarize tracking arrays and calculate predictive confidence boundaries
pva_summary <- data.frame(Year = proj_years, Median = apply(sim_reg_matrix, 2, median), Lower = apply(sim_reg_matrix, 2, function(x) quantile(x, 0.025)), Upper = apply(sim_reg_matrix, 2, function(x) quantile(x, 0.975)))
p_pva_traj <- ggplot(pva_summary, aes(x = Year, y = Median)) + geom_line(color = "purple4", linewidth = 1.2) + geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "purple", alpha = 0.15) + theme_bw() + labs(title = paste(target_species, "PVA forward abundance trajectory (Zero future take)"), y = "Total Population")
ggsave(paste0(out_dir, "/figure_2_pva_forward_abundance_trajectories.png"), p_pva_traj, width = 8, height = 5)

# Calculate long-term population collapse crash probability indices across tracking profiles
reg_thresh <- data.frame(Metric = paste(target_species, "Pooled Regional Cohort"), Risk_Below_50 = track_risk(sim_reg_matrix, n0_reg, 0.50), Risk_Below_25 = track_risk(sim_reg_matrix, n0_reg, 0.25), Risk_Below_12.5 = track_risk(sim_reg_matrix, n0_reg, 0.125))
write.csv(reg_thresh, paste0(out_dir, "/table_3_threshold_collapse_risk.csv"), row.names = FALSE)

# ==============================================================================
# Part 4: ESA incidental take model (Reconstruction bound to 1985-2015)
# ==============================================================================
message("--- Phase 4: Reconstructing historical take (2023 methods) ---")

# Strict filtering of file tracking bounds to keep analysis isolated to the 2015 limits
safe_df <- safe_df %>% filter(Year <= 2015)
obs_df  <- obs_df %>% filter(Year <= 2015)

# Standardize variable column text strings across dynamic custom entry formats
if ("Point.Estimate" %in% names(safe_df)) { safe_df$Total_Est <- as.numeric(safe_df$Point.Estimate)
} else if ("CC" %in% names(safe_df)) { safe_df$Total_Est <- as.numeric(safe_df$CC)
} else { safe_df$Total_Est <- 16.5 }

if ("SCL_raw" %in% names(obs_df)) { obs_df$Length <- as.numeric(obs_df$SCL_raw)
} else if ("Len" %in% names(obs_df)) { obs_df$Length <- as.numeric(obs_df$Len)
} else { obs_df$Length <- NA_real_ }

if ("M_mu" %in% names(obs_df)) { obs_df$Mortality <- as.numeric(obs_df$M_mu)
} else if (all(c("M.low", "M.high") %in% names(obs_df))) { obs_df$Mortality <- rowMeans(obs_df[, c("M.low", "M.high")], na.rm = TRUE)
} else { obs_df$Mortality <- NA_real_ }

# Calculate empirical average mortality parameters to handle missing data fields
emp_mean_mort <- mean(obs_df$Mortality, na.rm = TRUE)
if (is.nan(emp_mean_mort) | is.na(emp_mean_mort)) emp_mean_mort <- 0.535
obs_df$Mortality[is.na(obs_df$Mortality) | is.nan(obs_df$Mortality)] <- emp_mean_mort

# Join files together into unified data tracking matrices
merged_data <- obs_df %>% left_join(safe_df, by = "Year")
merged_data$Total_Est[is.na(merged_data$Total_Est)] <- 16.5

# Build predictive relationship profile between turtle length distributions and total annual capture rates
len_lm <- lm(log(Length) ~ Total_Est, data = merged_data)
sim_beta0 <- coef(len_lm)[1]; sim_beta1 <- ifelse(is.na(coef(len_lm)[2]), 0, coef(len_lm)[2])
sim_sigma_L <- summary(len_lm)$sigma; if(is.na(sim_sigma_L) || sim_sigma_L == 0) sim_sigma_L <- 0.338
sim_mu0 <- mean(safe_logit(merged_data$Mortality), na.rm=TRUE); sim_sigma_D <- sd(safe_logit(merged_data$Mortality), na.rm=TRUE)
sim_rho <- cor(log(merged_data$Length), safe_logit(merged_data$Mortality), use="complete.obs"); if(is.na(sim_rho)) sim_rho <- -0.51

# Construct formal variance-covariance array mapping length and mortality parameters together
sim_cov <- matrix(c(sim_sigma_L^2, sim_sigma_L * sim_sigma_D * sim_rho, sim_sigma_L * sim_sigma_D * sim_rho, sim_sigma_D^2), 2, 2)
emp_mean_len <- mean(obs_df$Length, na.rm = TRUE); if(is.nan(emp_mean_len)) emp_mean_len <- 117.5

sp_years <- sort(unique(safe_df$Year))
all_sp_turtles <- data.frame() # Clear baseline storage dataframe file

# Step through historical years to reconstruct the full characteristics of unobserved interactions
for (target_yr in sp_years) {
  obs_yr <- obs_df[which(obs_df$Year == target_yr), ]
  total_est <- sum(safe_df$Total_Est[safe_df$Year == target_yr], na.rm = TRUE)
  n_unobs <- max(0, round(total_est) - nrow(obs_yr)) # Quantify unobserved interactions missing measurements
  
  if (nrow(obs_yr) > 0) {
    # Save the existing empirical records directly to the master ledger file
    all_sp_turtles <- rbind(all_sp_turtles, data.frame(CaptureYear = target_yr, Length = coalesce(obs_yr$Length, emp_mean_len), Mortality = coalesce(obs_yr$Mortality, emp_mean_mort)))
  }
  if (n_unobs > 0) {
    # Generate joint draws for unobserved individual attributes using the covariance model matrices
    mu_l <- sim_beta0 + sim_beta1 * total_est
    draws <- mvtnorm::rmvnorm(n_unobs, mean = c(mu_l, sim_mu0), sigma = sim_cov)
    all_sp_turtles <- rbind(all_sp_turtles, data.frame(CaptureYear = target_yr, Length = exp(draws[,1]), Mortality = inv_logit(draws[,2])))
  }
}

sp_ledger <- data.frame() # Pre-allocate ledger file for forward projection tracks
# Run the historical age inversion and forward survival projection loops for every individual turtle
for (i in 1:nrow(all_sp_turtles)) {
  c_year <- all_sp_turtles$CaptureYear[i]
  # Invert von Bertalanffy model back from carapace length metrics to calculate absolute individual age estimates
  age_start <- sp_tknot - (1 / sp_k) * log(1 - (pmin(all_sp_turtles$Length[i], sp_linf - 0.1) / sp_linf))
  if(is.nan(age_start) || is.na(age_start)) age_start <- sp_max_age - 2 # Apply safety ceiling buffers
  
  future_years <- seq(c_year, 2015) # Generate timeline windows from capture point up to baseline envelope edge
  l_y <- length(future_years); if(l_y < 1) next
  
  ages_traj <- seq(age_start, length.out = l_y, by = 1) # Advance individual ages linearly
  lens_traj <- sp_linf * (1 - exp(-sp_k * (ages_traj - sp_tknot))) # Recalculate lengths tracking growth paths
  p_mat <- ifelse(lens_traj >= 0.99 * sp_linf, 1.0, 1.0 / (1.0 + exp(-(lens_traj - sp_lmat) / sp_sig_mat))) # Compute maturity probabilities
  
  p_binom <- rbinom(l_y, 1, p_mat) # Evaluate maturity state changes stochastically via binomial coin-flips
  if(any(p_binom == 1)) p_binom[min(which(p_binom == 1)):l_y] <- 1 # Lock turtle into absorbing mature state
  
  # Apply lifecycle matrix multipliers matching mature vs juvenile survival states
  surv_vector <- cumprod((1 - p_binom) * sp_ane_pj + p_binom * sp_ane_pa)
  # Calculate absolute reproductive contribution penalization value metrics (Annual Nester Equivalents)
  sr <- surv_vector * p_binom * sp_pf * all_sp_turtles$Mortality[i] * (1 / sp_remig_int)
  sp_ledger <- rbind(sp_ledger, data.frame(CalendarYear = future_years, ANE = sr))
}

# Calculate regional proportional composition layout weights across beaches
site_props <- abund_df %>% group_by(Site) %>% summarise(mean_n = mean(Annual_Nesters, na.rm = TRUE), .groups = "drop") %>% mutate(prop = mean_n / sum(mean_n, na.rm = TRUE))
final_historical_ane <- sp_ledger %>% group_by(CalendarYear) %>% summarise(Total_Cumulative_ANE = sum(ANE, na.rm=TRUE), .groups = "drop")

# Reconstruct a counterfactual adjusted history matrix by adding lost reproductive counts back into dataset
abund_adj <- abund_df %>%
  left_join(final_historical_ane, by = c("Year" = "CalendarYear")) %>%
  mutate(Total_Cumulative_ANE = ifelse(is.na(Total_Cumulative_ANE), 0, Total_Cumulative_ANE)) %>%
  left_join(site_props, by = "Site") %>%
  mutate(Annual_Nesters = Annual_Nesters + (Total_Cumulative_ANE * prop)) %>% # Add lost turtles back to background tracks
  select(-Total_Cumulative_ANE, -prop, -mean_n)

message("--- Running trend engine on fishery adjusted data ---")
trend_adj <- run_trend_model(abund_adj, beach_cols) # Re-fit state-space model to counterfactual matrix

# Extract counterfactual adult population size draws matching terminal baseline year (2015)
X_post_adj <- trend_adj$sims.list$X; total_pop_post_adj <- exp(X_post_adj)
for(j in 2:length(beach_cols)) { total_pop_post_adj <- total_pop_post_adj + exp(X_post_adj + trend_adj$sims.list$A[, j]) }
total_pop_post_adj <- total_pop_post_adj * sp_remig_int; start_notake_draws <- total_pop_post_adj[, 31]

# Re-sample model posterior vectors to configure dual comparative projection tracks
r_notake <- sample(trend_adj$sims.list$U, n_sims, replace = TRUE)
q_notake <- sample(trend_adj$sims.list$Q, n_sims, replace = TRUE)
r_take   <- sample(trend_base$sims.list$U, n_sims, replace = TRUE)
q_take   <- sample(trend_base$sims.list$Q, n_sims, replace = TRUE)
n0_notake <- sample(start_notake_draws, n_sims, replace = TRUE)

sim_notake <- matrix(NA, n_sims, horizon); sim_take <- matrix(NA, n_sims, horizon)

# Execute final comparative dual-scenario forward projection matrices (Take vs No-Take)
message("--- Projecting dual fishery scenarios ---")
pb_dual <- txtProgressBar(min = 0, max = n_sims, style = 3)
for(i in 1:n_sims) {
  curr_nt <- n0_notake[i] # Counterfactual trajectory initial population anchor track
  curr_t <- n0_reg[i]     # True observed trajectory initial population anchor track
  take_ane_vector <- rep(0, horizon) # Array tracking distributed forward reproductive losses over time
  
  for(y in 1:horizon) {
    # Advance counterfactual un-impacted profile track forward into next timeline step
    curr_nt <- rnorm(1, curr_nt * exp(r_notake[i]), sqrt(q_notake[i])); sim_notake[i,y] <- curr_nt
    
    # Draw stochastically simulated future interactions from CMP distribution models
    takes <- rCMP(1, mu = cmp_mu, nu = cmp_nu)
    if (takes > 0) {
      # Reconstruct lengths and mortality profiles across future interaction events using multivariate systems
      mu_l <- sim_beta0 + sim_beta1 * takes
      draws <- mvtnorm::rmvnorm(takes, mean = c(mu_l, sim_mu0), sigma = sim_cov)
      for(t in 1:takes) {
        # Execute demographic age conversions and track forward lifecycle paths into future projection horizons
        age_t <- sp_tknot - (1 / sp_k) * log(1 - (pmin(exp(draws[t,1]), sp_linf - 0.1) / sp_linf))
        if(is.nan(age_t) || is.na(age_t)) age_t <- sp_max_age - 2
        
        rem_yrs <- y:horizon; l_y <- length(rem_yrs)
        ages_traj <- seq(age_t, length.out = l_y, by = 1)
        lens_traj <- ifelse(ages_traj >= sp_max_age, rep(sp_linf * 0.99, l_y), sp_linf * (1 - exp(-sp_k * (ages_traj - sp_tknot))))
        p_mat <- ifelse(lens_traj >= 0.99 * sp_linf, 1.0, 1.0 / (1.0 + exp(-(lens_traj - sp_lmat) / sp_sig_mat)))
        p_binom <- rbinom(l_y, 1, p_mat); if(any(p_binom == 1)) p_binom[min(which(p_binom == 1)):l_y] <- 1
        
        surv_vector <- cumprod((1 - p_binom) * sp_ane_pj + p_binom * sp_ane_pa)
        sr <- surv_vector * p_binom * sp_pf * inv_logit(draws[t,2]) * (1 / sp_remig_int)
        
        # Distribute annual reproductive value loss dynamically over remaining future years
        take_ane_vector[y:horizon] <- take_ane_vector[y:horizon] + sr
      }
    }
    # Subtract realized dynamic loss profiles from baseline trend trajectories
    curr_t <- rnorm(1, (curr_t - take_ane_vector[y]) * exp(r_take[i]), sqrt(q_take[i])); sim_take[i,y] <- curr_t
  }
  # Extinction checks ensuring tracking channels drop cleanly to zero rather than negative numeric space
  if(any(sim_notake[i,] < 0)) { extinction_yr <- min(which(sim_notake[i,] < 0)); sim_notake[i, extinction_yr:horizon] <- 0 }
  if(any(sim_take[i,] < 0)) { extinction_yr <- min(which(sim_take[i,] < 0)); sim_take[i, extinction_yr:horizon] <- 0 }
  setTxtProgressBar(pb_dual, i)
}
close(pb_dual)

# Save final graphic figures and structural logs into relative output subdirectory
yrs_proj <- (max(abund_df$Year) + 1):(max(abund_df$Year) + horizon)
df_notake <- calc_summary(sim_notake) %>% mutate(Scenario = "No Take (Nj + F)")
df_take   <- calc_summary(sim_take) %>% mutate(Scenario = "Take (Nj)")
df_diff   <- calc_summary(sim_notake - sim_take)

p_overlay <- ggplot(bind_rows(df_notake, df_take), aes(x = Year, y = Median, color = Scenario, fill = Scenario)) + geom_line(linewidth = 1.2) + geom_ribbon(aes(ymin = L95, ymax = U95), alpha = 0.15, color = NA) + scale_color_manual(values = c("No Take (Nj + F)" = "chartreuse4", "Take (Nj)" = "dodgerblue4")) + scale_fill_manual(values = c("No Take (Nj + F)" = "chartreuse3", "Take (Nj)" = "dodgerblue3")) + theme_bw() + labs(title = paste("2023 methods validation: PVA projection overlay (", target_species, ")", sep=""), y = "Total Population")
ggsave(paste0(out_dir, "/figure_3_take_vs_notake_baseline.png"), p_overlay, width = 8, height = 5)

p_diff <- ggplot(df_diff, aes(x = Year, y = Median)) + geom_line(color = "darkorchid4", linewidth = 1.2) + geom_ribbon(aes(ymin = L95, ymax = U95), fill = "darkorchid", alpha = 0.25) + theme_bw() + labs(title = paste("Isolated deep-set fleet cumulative penalty profile (", target_species, ")", sep=""), y = "\u0394 Total Population")
ggsave(paste0(out_dir, "/figure_4_isolated_fishery_impact_difference.png"), p_diff, width = 8, height = 4)

write.csv(data.frame(Baseline = c("True State (Nj)", "Fishery-Adjusted (Nj + F)"), Trend_r = c(mean(trend_base$sims.list$U), mean(trend_adj$sims.list$U)), Process_Variance_Q = c(mean(trend_base$sims.list$Q), mean(trend_adj$sims.list$Q))), paste0(out_dir, "/table_4_fishery_ledger.csv"), row.names = FALSE)

# ==============================================================================
# Automated text generation console display execution
# ==============================================================================
cat("\n==============================================================================\n")
cat("               Loggerhead extraction engine: Output for 1985-2015 data        \n")
cat("==============================================================================\n\n")

pct_change_base <- (exp(mean(trend_base$sims.list$U)) - 1) * 100; pct_cri_base <- (exp(quantile(trend_base$sims.list$U, probs = c(0.025, 0.975))) - 1) * 100
diff_matrix_cc <- sim_notake - sim_take

cat("### Draft text generation ###\n\n")
cat(sprintf("In the updated analysis, the estimated long-term trend for North Pacific loggerheads was %.2f%% (95%% CrI: %.2f%% to %.2f%%).\n", pct_change_base, pct_cri_base[1], pct_cri_base[2]))
cat(sprintf("Current annual nester abundance was estimated at %.1f adults in the active cohort (95%% CrI: %.1f to %.1f).\n\n", median(start_take_draws), quantile(start_take_draws, 0.025), quantile(start_take_draws, 0.975)))
cat(sprintf("Differences between take and no-take scenarios were %.2f adult female equivalents (95%% CrI: %.2f to %.2f) over the short-term projection period (10 years) and %.2f (95%% CrI: %.2f to %.2f) over the 100-year projection period.\n", median(diff_matrix_cc[, 10]), quantile(diff_matrix_cc[, 10], 0.025), quantile(diff_matrix_cc[, 10], 0.975), median(diff_matrix_cc[, 100]), quantile(diff_matrix_cc[, 100], 0.025), quantile(diff_matrix_cc[, 100], 0.975)))
cat("==============================================================================\n")