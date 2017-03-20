## Description
## -----------
##
## Script performs a sensitivity analysis on the R muscosa-Bd model.  Using
## the methods outlined in Sobie 2009 and Harper et al. 2011, we perform 
## a proportional, global sensitivity analysis of the the model to the 
## within-season, disease-related parameters. We perturb all parameters 
## simultaneously with a random draw from a lognormal distribution centered at
## 1 and then simulate the model for 8 years and record whether or not
## the population went extinction in this time frame. We then save the perturbed
## parameters and the extinction response from each simulation and analyze it
## using various classification methods to identify parameter importance.
##
## Author: Mark Wilber

# Source the necessary functions
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

# Depending on the transmission function, drop some of th parameters that are
# not relevant
drop_vals = list(15:18, 17:19, c(15, 16, 19))

# Number of simulations
SIMS = 1000

# Arrays to hold results
param_matrix = data.frame(matrix(NA, nrow=SIMS, ncol=length(param_names)))
colnames(param_matrix) = param_names
extinction_matrix = array(NA, dim=SIMS)

# No variation in resistance
standing_var = 0.001
resist_bins = 1

extinction_sens = function(i, sigma){

    # Perturb the parameters
    print(paste("Beginning simulation", i))
    pfull_up = perturb_params(pfull, param_names, sigma=sigma)

    sim_res = multiseason_simulation(pfull_up, pother, standing_var, resist_bins,
                         sims=1, output=FALSE)
    ext_times = simulation_extinction_times(sim_res)

    return(list(params=unlist(pfull_up[param_names]), 
                ext=!all(!ext_times)))

}

# Different sigmas to check consistency of sensitivity analysis
sigmas = c(0.3) 

# Uncomment to try sensitivity analysis with a range of different sigmas
# sigmas = c(0.1, 0.2, 0.3, 0.4, 0.5)

for(sigma in sigmas){

    print(paste('Starting sigma', sigma))

    for(j in 1:length(transmission_fxns)){

        trans = transmission_fxns[[j]]
        print(paste("Working on", trans))

        # Set the transmission function
        pother$ipm$transmission_type = trans

        res = mclapply(1:SIMS, extinction_sens, sigma, mc.cores=3)

        # # Unpack results and save them to a csv
        X = do.call(rbind, lapply(res, function(x) x$params))
        X = X[, -drop_vals[[j]]]
        y = unlist(lapply(res, function(x) x$ext))

        dat = as.data.frame(X)
        dat$y = as.numeric(y)
        write.csv(dat, paste("../../results/sens_results_", trans, 
                                                "_", sigma, ".csv", sep=""))

    }

}