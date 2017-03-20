## Description
## -----------
##
## Dynamics and extinction plots for when zoospore decay rate in the
## environment is high (survival probability is low).
##
## Simulating the dynamical predictions if zoospores decay very quickly
## and extracting the resulting extinction curve.
##
## Author: Mark Wilber
##
## Dependencies
## ------------
## 1. IPM_functions_for_R.R
## 2. model_parameter.yml
## 3. ../../results/bayesian_parameter_estimates.rds

# Source necessary simulation functions
source("IPM_functions_for_R.R")
library(yaml)

# Load in IPM parameters
pfull = readRDS("../../results/bayesian_parameter_estimates.rds")
pother = yaml.load_file("model_parameters.yml")
pother$sim$bins = 30 # For speed
pother$sim$years= 25
pother$ipm$transmission_type = "zoospore_pool_infection_dd"

# Make dynamics plot
sims = 500
res = multiseason_simulation(pfull, pother, 1, 1, sims=sims, zsurv=0.05)

# Plot some simulations

plot_it = FALSE
if(plot_it){
    # Plot trajectories with high zoospore decay

    plot_simulation(res[20:30], pfull, pother)
    dev.copy(pdf, paste("../../results/zoospore_decay_trajectories.pdf", sep=""))
    dev.off()
}


# Make extinction plot
ext_vals = simulation_extinction_times(res)
ext = apply(ext_vals, 2, mean)
steps = with(pother$sim, (days_per_year * years) / step_length)
time = with(pother$sim, (c(0, 1:steps) * step_length) / days_per_year)
dat = data.frame(years=time, means=ext)
dat$zoospore_decay = "Constant, high rate of decay"

# Save the results
write.csv(dat, "../../results/zsurv_high_traj.csv", row.names=FALSE)
 

# Load an plot the other simulation results
# full_sim=readRDS("../../results/extinction_simulations_full_bayes_pres.rds")
# dat_zp = data.frame(full_sim$zoospore_pool_infection_fd)
# dat_fd = data.frame(full_sim$zoospore_pool_infection_dd)
# dat_zp$zoospore_decay = "Empirically estimated decay"

# full_dat = rbind(dat, dat_zp)

# # Plot and save zoospore decay data
# ggplot(full_dat, aes(x=years, y=means)) + geom_line(aes(color=zoospore_decay)) + 
#                         theme_bw() + theme(legend.position=c(.8, .1)) + 
#                         xlab("Time (years)") +
#                         ylab("Extinction probability")
# ggsave("../../results/compare_extinction_zdecay.pdf")





