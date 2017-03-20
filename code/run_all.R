## Description
## -----------
##
## R script to run all of the analyses described in 
## `../docs/ipm_extinction_*.pdf`.
##
##
## WARNING: With the default number of simulations and cores used in the 
## analysis the `run_all.R` script may take >24 hours to run.  
## If you want to speed this up in order to play around with how the scripts 
## are working (but, of course, at the cost of not getting a good summary of 
## the stochastic simulations!), make the following changes to reduce the 
## number of simulations.

## 1. `sensitivity_analysis.R`: Line 50 -> Change to SIMS = 10
## 2. `sensitivity_analysis_bayesian.R`:  Line 50 -> Change to SIMS = 10
## 3. `zdecay_plots.R`: Line 30 -> Change to sims = 10
## 4. `stochastic_multiseason_simulation.R`: Line 40 -> Change to SIMS = 10


# Will need to change this working directory for your particular machine
cwd = "~/Repos/density_dependent_ipm_for_github/code"
setwd(cwd)

# Check for required packages and install them if necessary
list.of.packages <- c("fields", "grid", "popbio", "reshape2", "plyr", 
                      "AICcmodavg", "GGally", "MASS", "RColorBrewer", "abind",
                      "boot", "data.table", "ggplot2", "gridExtra", "lme4", 
                      "nlme", "parallel", "plyr", "rethinking", "rpart",
                      "rpart.plot", "rstan", "scales", "tree", "yaml")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

############ Begin the analysis ######################

# Fit all the vital rate functions and transmission parameters
print("Running Bayesian analysis on vital rate functions")
setwd("stat_analysis/")
source("ipm_vital_rates_bayesian.R")
setwd("..")
print("Completed Bayesian parameter analysis of vital rate functions.")

# Compute R0
print("Performing R0 calculation for different densities and temperatures...")
setwd("model_analysis")
source("R0_environmental.R")
setwd("..")
print("Completed R0 calculations")

# Generate the extinction curves with the hybrid model. Also get the example
# trajectories given in Figure 3 in `../docs/ipm_extinction_manuscript.pdf`
print("Simulating hybrid model to generate extinction curves...")
setwd("model_analysis/")
source("stochastic_multiseason_simulation.R")
source("example_trajectories.R")
source("zdecay_plots.R")
setwd("..")
print("Completed simulation of hybrid model")

# Perform the sensitivity analysis
print("Performing sensitivity analysis with lognormal perturbations...")
setwd("model_analysis/")
source("sensitivity_analysis.R")
print("Performing sensitivity analysis with posterior perturbations...")
source("sensitivity_analysis_bayesian.R")
setwd("..")
print("Completed sensitivity analysis")

print("Analysis complete")

## See stat_analysis/manuscript_sensitivity_and_plots.ipynb for IPython
## notebook containing the code to make the plots displayed in 
## `../docs/ipm_extinction_*.pdf`. 
