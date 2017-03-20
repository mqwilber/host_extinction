## Description
## -----------
##
## Script that performs stochastic simulations of R. muscosa-Bd IPM model.
## Assuming no heterogeneity in susceptibility, this script simulates a
## temperature-dependent IPM with three different transmission functions
## parameterized from data.  For each transmission function, it computes an 
## extinction curve over 30 years.  It also uses a multitype branching process
## approximation to compute the sensitivity of the extinction curves to the 
## various vital rate functions in the IPM.  These results are the basis for
## figure 4 in the manuscript.
##
## Author: Mark Wilber

# Source the IPM functions
source("IPM_functions_for_R.R")
library(yaml)
library(scales)
library(abind)
library(ggplot2)
library(parallel)

# Load in parameters for the IPM model 
pfull = readRDS("../../results/bayesian_parameter_estimates.rds")#readRDS("../../results/IPM_parameters_no_26_linear.rds")
pother = yaml.load_file("model_parameters.yml")
pother$sim$years= 25

# Best fitting transmission functions
transmission_fxns = list(
            zoospore_pool_infection_const="zoospore_pool_infection_const",
            zoospore_pool_infection_fd="zoospore_pool_infection_fd",
            zoospore_pool_infection_dd="zoospore_pool_infection_dd")

# Parameters 
standing_var = 0.001 # Variation in growth slope
resist_bins = 1 # Number of resistance classes (1 if there is no heterogeneity)
SIMS = 500 # Number of simulations 
elas_vals = list(elasphi=1,
                 elasg0=1,
                 elasg=1,
                 elass=1,
                 elasr=1,
                 elasz=1)
delta = 0.001

sim_function = function(trans, elas_vals, method){
    # Runs multiseason simulation and computes extinction metrics
    # "annual": Uses annual multitype branching process approximation
    # "daily": Uses daily multitype branching process approximation
    # "full": Uses the full stochastic simulation

    pother$ipm$transmission_type = trans

    if(method == "annual"){ # Annual branching process

        prob_extinct = extinction_prob_annual(pfull, pother, standing_var, 
                                resist_bins, elas_vals=elas_vals)
        means = prob_extinct
        years = 0:(length(prob_extinct) - 1)

    } else if(method == "daily"){ # Daily branching process

        prob_extinct = extinction_prob_daily(pfull, pother, standing_var, 
                                resist_bins, elas_vals=elas_vals)
        years = prob_extinct$time
        means = prob_extinct$probs

    } else { # Full stochastic simulation

        sim_res = multiseason_simulation(pfull, pother, standing_var, resist_bins, 
                            elas_vals=elas_vals, sims=SIMS)
        extinction_times = simulation_extinction_times(sim_res)
        means = apply(extinction_times, 2, mean)
        years = (0:(length(means) - 1)) * pother$sim$step_length / pother$sim$days_per_year
    }

    return(list(years=years, means=means))
}

## This code generates the extinction curves for the full density-dependent
## simulation

# Run full stochastic simulation if necessary
fname_full = "../../results/extinction_simulations_full_bayes.rds"
if(length(Sys.glob(fname_full)) == 0){

    anal_results_full = mclapply(transmission_fxns, sim_function, elas_vals, "full", 
                                    mc.cores=3)
    saveRDS(anal_results_full, fname_full)
}


## This code assumes a density-independent situation so that we can use
## a branching process approximation to look at the local sensitivity of 
## extinction risk to perturbations in the vital rate functions.

# Check if base simulation file exists. If so don't run the base simulation
fname_annual = "../../results/extinction_simulations_annual_bayes.rds"
if(length(Sys.glob(fname_annual)) == 0){

    # Multiprocess loop for increased speed
    anal_results = mclapply(transmission_fxns, sim_function, elas_vals, "annual", 
                                    mc.cores=3)

    # Save simulation results as rds
    saveRDS(anal_results, fname_annual)

} else{
    anal_results = readRDS(fname_annual)
}

# Compute the fitted parameter AUC extinction estimates
base_aucs = list()
for(i in 1:length(transmission_fxns)){

    trans = transmission_fxns[[i]]
    base_aucs[[trans]] = pracma::trapz(anal_results[[trans]]$years, 
                                            anal_results[[trans]]$means)

}

## Compute the elasticities of the AUC estimates to the vital rate fxns ###

vitalrate_elasfxn = function(j, elas_vals, fitted_auc, trans, delta){
    # Function to parallelize the AUC elasticity analysis
    # 

    vital_fxn = names(elas_vals)[j]
    updated_evals = elas_vals
    updated_evals[[vital_fxn]] = elas_vals[[vital_fxn]] - delta # Perturb
    tsim_res = sim_function(trans, updated_evals, "annual") # Simulate with perturbation

    elas_auc = pracma::trapz(tsim_res$years, tsim_res$means) # Calc AUC
    telas_val = (elas_auc - fitted_auc) / (delta * fitted_auc) # Calc elas
    return(telas_val)

}

auc_elasticity = list()
for(i in 1:length(transmission_fxns)){ # Loop through transmission functions
    
    trans = transmission_fxns[[i]]
    fitted_auc = base_aucs[[trans]]
    auc_elasticity[[trans]] = list()
    print(paste("Working on", trans))

    # Perturb each vital rate fxn in turn. Using parallel processing
    auc_elasticity[[trans]] = mclapply(seq_along(elas_vals), 
                                vitalrate_elasfxn, elas_vals, fitted_auc,
                                trans, delta, mc.cores=3)

    # Rename the elasticity functions
    names(auc_elasticity[[trans]]) = names(elas_vals)

}

# Save the sensitivity results
saveRDS(auc_elasticity, "../../results/auc_elasticity_analysis_bayes.rds")



