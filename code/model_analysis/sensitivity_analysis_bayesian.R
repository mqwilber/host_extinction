## Description
## -----------
##
## Perform a sensitivity analysis in which we are drawing from the posterior
## distributions of the vital rate functions.  This is accounting for the
## correlation between parameter values and the experimental uncertainty. This 
## allows some parameters to have larger ranges of uncertainty.

## The other approach that is taken in `sensitivity_analysis.R` simply perturbs
## all the parameters the same amount, but doesn't account for correlation between
## the parameters.  The advantage the approach in `sensitivity_analysis.R` is 
## that we can assume all parameters are independent to look at each parameter
## independently. This can help us isolate effects.  On the other hand it does
## ignore correlation and uncertainty in the parameter estimates.

## Author: Mark Wilber

source("IPM_functions_for_R.R")
library(yaml)
library(scales)
library(abind)
library(ggplot2)
library(parallel)

# Load in IPM parameters
pfull = readRDS("../../results/bayesian_parameter_estimates.rds") #readRDS("../../results/IPM_parameters_no_26_linear.rds")
pother = yaml.load_file("model_parameters.yml")
pother$sim$years = 8
pother$sim$bins = 30

# List of parameters to perturb in pfull
param_names = c("surv_int", "surv_slope",
                "growth_int", "growth_temp", "growth_size", 
                "growth_sigma2", "growth_sigma2_exp",
                "clump_int", "clump_temp", "clump_sigma2", "clump_sigma2_exp",
                "loss_int", "loss_temp", "loss_size", 
                "zfd_intercept", "zfd_slope",
                "zdd_intercept", "zdd_slope",
                "zconst_intercept")


transmission_fxns = list("zoospore_pool_infection_const",
            "zoospore_pool_infection_fd",
            "zoospore_pool_infection_dd")

# Drop the necessary values depending on the transmission function
drop_vals = list(15:18, 17:19, c(15, 16, 19))

# Simulation results
SIMS = 1000

# Arrays to hold results
param_matrix = data.frame(matrix(NA, nrow=SIMS, ncol=length(param_names)))
colnames(param_matrix) = param_names
extinction_matrix = array(NA, dim=SIMS)

standing_var = 0.001
resist_bins = 1

extinction_sens = function(i, sigma){

    # Perturb the parameters
    print(paste("Beginning simulation", i))
    pfull_up = bayesian_params(pfull) # Bayesian draw

    sim_res = multiseason_simulation(pfull_up, pother, standing_var, resist_bins,
                         sims=1, output=FALSE)
    ext_times = simulation_extinction_times(sim_res)

    return(list(params=unlist(pfull_up[param_names]), 
                ext=!all(!ext_times)))

}

print("Beginning analysis...")

for(j in 1:length(transmission_fxns)){

    trans = transmission_fxns[[j]]
    print(paste("Working on", trans))

    # Set the transmission function
    pother$ipm$transmission_type = trans

    res = mclapply(1:SIMS, extinction_sens, sigma, mc.cores=4)

    # # Unpack results and save them to a csv
    X = do.call(rbind, lapply(res, function(x) x$params))
    X = X[, -drop_vals[[j]]]
    y = unlist(lapply(res, function(x) x$ext))

    dat = as.data.frame(X)
    dat$y = as.numeric(y)
    write.csv(dat, paste("../../results/sens_results_", trans, 
                                            "_bayesian.csv", sep=""))

}

print("Completed analysis")
