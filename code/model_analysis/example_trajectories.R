## Description
## -----------
## 
## Script runs hybrid model simulations with different transmission functions
## and plots the resulting trajectories.  All plots are saved into the results
## folder.
##
## Author: Mark Wilber
## 
## Dependencies
## ------------
## 1. IPM_functions_for_R.R
## 2. model_parameter.yml
## 3. ../../results/bayesian_parameter_estimates.rds 

source("IPM_functions_for_R.R")
library(yaml)

# Load in IPM parameters
pfull = readRDS("../../results/bayesian_parameter_estimates.rds")
pother = yaml.load_file("model_parameters.yml")

# NOTE: You can adjust the various parameters in pother (see model_parameter.yml)
# to see how various assumptions affect the simulation trajectories.

# Strings specifying the types of transmission functions to use.
# see `get_infection_probability` function in IPM_functions_for_R.R
# for a full list of transmission functions to explore.
trans_types = c(
            "zoospore_pool_infection_const",
            "zoospore_pool_infection_fd",
            "zoospore_pool_infection_dd")

# Loop through transmission functions, run simulation, plot and save
# results
for(trans in trans_types){

    pother$ipm$transmission_type = trans
    # pfull$dd_intercept = 0.001
    pother$sim$years = 20

    res = multiseason_simulation(pfull, pother, 1, 1, sims=5)
    plot_simulation(res, pfull, pother)
    dev.copy(pdf, paste("../../results/trajectories_", trans, ".pdf", sep=""))
    dev.off()
}




