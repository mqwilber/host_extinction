## Description
## ------------
##
## Refitting all of the vital rate functions using a Bayesian framework
## Taking the best fitting models found in ipm_vital_rate_parameter_analysis.R
## and refitting them in a Bayesian framework for easier and more accurate 
## propagation of uncertainty.
##
## Author: Mark Wilber

source("../model_analysis/IPM_functions_for_R.R")
library(rethinking)
library(rstan)

### Loading and formatting temperature data ###

# Load in temperature data
temp_data = read.csv("../../data/formatted/converted_temp_data.csv")

# Log-transform all Bd loads but zeros so that they can techinically range from
# -inf to inf.  This will be important when we are dealing with eviction
# problems later on

log_trans_data = function(data){

  # Function to log-transform all data Bd-data but not the zeros

  nz_ind = data$size != 0
  data$size[nz_ind] = log(data$size[nz_ind])
  nznext_ind = (data$sizeNext != 0) & !(is.na(data$sizeNext))
  data$sizeNext[nznext_ind] = log(data$sizeNext[nznext_ind])

  return(data)

}

temp_data = log_trans_data(temp_data)

# The time between successive observations on a frog in 3 days.  A few time
# points were longer to let's exclude those

days = 3
temp_data = temp_data[temp_data$time == 3, ]

# A set of data without the temperature 26
temp_data_no_26 = temp_data[temp_data$temp != 26, ]

##############################################################################

# Set up a data.frame to hold the parameters of the vital functions
params=list() # Holds linear parameters

drop_zeros = function(data){

  # Function to remove zeros from size and sizeNext

  data_nozeros = data[data$size != 0, ]

  # Drop individuals that lose infection as well as we are considering this a separate process
  ind = data_nozeros$sizeNext != 0
  ind[is.na(ind)] = TRUE
  data_nozeros = data_nozeros[ind, ]
  return(data_nozeros)

}

temp_data_nozeros = drop_zeros(temp_data)
temp_data_no_26_nozeros = drop_zeros(temp_data_no_26)

### Fitting the survival function, s(x) ###

surv_data = temp_data_no_26_nozeros[temp_data_no_26_nozeros$temp == 20, ]
surv_data_stan = list(N=nrow(surv_data), surv=surv_data$surv, 
                        bd_load=as.vector(scale(surv_data$size)))

surv_mod = stan("stan_files/survival_function.stan", data=surv_data_stan, iter=5000,
                                chains=3)

# Saving the parameters and samples

surv_samples = extract(surv_mod)
params$surv_int = mean(surv_samples[['b0']])
params$surv_slope = mean(surv_samples[['b1']])
params$surv_samples = as.data.frame(do.call(cbind, surv_samples[1:2])[
                              sample(1:length(surv_samples[['b0']]), 1000), ])

colnames(params$surv_samples) = c("surv_int", "surv_slope")

# For scaling 
params$surv_mean = mean(temp_data_no_26_nozeros[temp_data_no_26_nozeros$temp == 20, ]$size)
params$surv_sd = sd(temp_data_no_26_nozeros[temp_data_no_26_nozeros$temp == 20, ]$size)

### Fitting the growth function, G(x', x) ###

temp_data_no_26_nozeros_nona = temp_data_no_26_nozeros[!is.na(temp_data_no_26_nozeros$sizeNext), ]

growth_data = temp_data_no_26_nozeros_nona
growth_data_stan = list(N=nrow(growth_data), temp=growth_data$temp, 
                            sizeNow=growth_data$size, 
                            sizeNext=growth_data$sizeNext)

# Fitting the best model identified in ipm_vital_rate_parameter_analysis
growth_mod = stan("stan_files/growth_function.stan", data=growth_data_stan,
                        iter=5000, chains=3)

growth_samples = extract(growth_mod)
params$growth_int = mean(growth_samples[['b0']])
params$growth_size = mean(growth_samples[['b1']])
params$growth_temp = mean(growth_samples[['b2']])
params$growth_sigma2 = mean(growth_samples[['sigma']])
params$growth_sigma2_exp = mean(growth_samples[['delta']])
params$growth_samples = as.data.frame(do.call(cbind, growth_samples[1:5])[
                            sample(1:length(growth_samples[['b0']]), 1000), ])
colnames(params$growth_samples) = c('growth_int', 'growth_size', 
                                    'growth_temp', 'growth_sigma2', 
                                    'growth_sigma2_exp')


### Fitting the initial infection function, G_0(x') ###

alive = temp_data[temp_data$surv == 1, ]
alive_no_26 = temp_data_no_26[temp_data_no_26$surv == 1, ]

transition_fxn = function(data){

  # Find the frogs who went from 0 load to non-zero load
  transition = data[data$size == 0 & data$sizeNext != 0, ]

  ind = (transition$sizeNext > 5)
  print(paste("Dropping", sum(ind), "points"))

  dropped = transition[ind,]

  # Drop values greater than 6 as they are likely PCR error
  transition = transition[!ind, ]

  return(list(transition=transition, dropped=dropped))

}

transition = transition_fxn(alive)$transition
transition_no_26 = transition_fxn(alive_no_26)$transition

clump_data = transition_no_26
clump_data_stan = list(N=nrow(clump_data), sizeNext=clump_data$sizeNext,
                        temp=clump_data$temp)

clump_mod = stan("stan_files/initial_infection_function.stan", 
                                data=clump_data_stan, iter=5000,
                                chains=3) 

# Save the parameters from the clump model
clump_samples = extract(clump_mod)
params$clump_int = mean(clump_samples[['b0']])
params$clump_temp = mean(clump_samples[['b1']])
params$clump_sigma2 = mean(clump_samples[['sigma']])
params$clump_sigma2_exp = mean(clump_samples[['delta']])
params$clump_samples = as.data.frame(do.call(cbind, clump_samples[1:4])[
                            sample(1:length(clump_samples[['b0']]), 1000), ])
colnames(params$clump_samples) = c("clump_int", "clump_temp", "clump_sigma2",
                                    "clump_sigma2_exp")


### Fitting the loss of infection function, l(x) ###

alive_fxn = function(data){
  # Function for preparing data for analysis

  data['loss'] = 0

  # Don't include the zero class
  alive_trun = data[data$size != 0, ]

  alive_trun$loss[alive_trun$sizeNext == 0] = 1

  # Drop cases where loss occurs over 6 ZEs...probably PCR error
  ind = (alive_trun$size > 5 & alive_trun$loss == 1)
  dropped = alive_trun[ind, ]

  print(paste("Dropping", sum(ind), "points"))

  alive_trun = alive_trun[!ind, ]

  return(list(alive_trun=alive_trun, dropped=dropped))

}

alive_trun = alive_fxn(alive)$alive_trun
alive_trun_no_26 = alive_fxn(alive_no_26)$alive_trun


loss_data = alive_trun_no_26
loss_data_stan = list(N=nrow(loss_data), loss=loss_data$loss, 
                        sizeNow=loss_data$size, temp=loss_data$temp)

loss_mod = stan("stan_files/loss_function.stan", data=loss_data_stan,
                          iter=5000, chains=3)

# Saving the parameters for the loss of infection function
loss_samples = extract(loss_mod)
params$loss_int = mean(loss_samples[['b0']])
params$loss_size = mean(loss_samples[['b1']])
params$loss_temp = mean(loss_samples[['b2']])
params$loss_samples = as.data.frame(do.call(cbind, loss_samples[1:3])[
                            sample(1:length(loss_samples[['b0']]), 1000), ])
colnames(params$loss_samples) = c("loss_int", "loss_size", "loss_temp")

### Fitting the transmission functions, phi ### 

# The testing of various models is done in the file `transmission_analysis.R`
# and `pomp_transmission.R`. 
# Here we are just refitting and saving the best models from these analyses

dat = read.csv("../../data/formatted/formatted_dd_data.csv", header=T)
short_dat = dat[dat$Day < 35, ] # Before the crash

# Only look at transitions with times between 4 and 8 days
short_dat = short_dat[(short_dat$time > 3) & (short_dat$time < 9), ]
short_dat = short_dat[complete.cases(short_dat), ] # Drop NAs
short_dat$time_ratio = short_dat$time #/ min(short_dat$time)

# Subset to get the transmission data. 0 -> infected/not infected
trans_data = short_dat[short_dat$size == 0, ]
trans_data$pa = as.integer(trans_data$sizeNext > 0)

## Model 1: Fit density-dependent transmission model with constant time intercept ##

dd_data_stan = list(N=nrow(trans_data),
               infected=trans_data$pa,
               dens_or_prop=trans_data$num_infected*trans_data$time_ratio,
               time=trans_data$time_ratio)

dd_mod = stan("stan_files/fd_dd_transmission_function.stan", data=dd_data_stan,
                      iter=5000, chains=3)

dd_samples = extract(dd_mod)
params$dd_intercept = mean(dd_samples[['b0']])
params$dd_slope = mean(dd_samples[['b1']])
params$dd_samples = as.data.frame(do.call(cbind, dd_samples[1:2])[
                            sample(1:length(dd_samples[['b0']]), 1000), ])
colnames(params$dd_samples) = c("dd_intercept", "dd_slope")

## Model 2: Fitting frequency-dependent transmission model with constant time intercept ##

prop_inf = trans_data$num_infected / trans_data$density
trans_data$prop_inf = prop_inf

fd_data_stan = list(N=nrow(trans_data),
               infected=trans_data$pa,
               dens_or_prop=trans_data$prop_inf*trans_data$time_ratio,
               time=trans_data$time_ratio)

fd_mod = stan("stan_files/fd_dd_transmission_function.stan", data=fd_data_stan,
                      iter=5000, chains=3)

fd_samples = extract(fd_mod)
params$fd_intercept = mean(fd_samples[['b0']])
params$fd_slope = mean(fd_samples[['b1']])
params$fd_samples = as.data.frame(do.call(cbind, fd_samples[1:2])[
                            sample(1:length(fd_samples[['b0']]), 1000), ])
colnames(params$fd_samples) = c("fd_intercept", "fd_slope")


# Model 3: Fitting the constant infection model
const_mod = stan("stan_files/prob_inf_constant.stan", data=fd_data_stan,
                      iter=5000, chains=3)

## Fitting the transmission models with dynamic zoospore pool

# Source the pomp models
source("pomp_transmission.R") # File contains the model fitting
list2env(models_in_manuscript, envir=environment()) # Put results in current environment

# Model 4: Frequency-dependent dynamic zoospore pool
zfd_samples = extract(pomp_mod_fd)
params$zfd_intercept = mean(zfd_samples[['b_zpool']])
params$zfd_slope = mean(zfd_samples[['b_density']])
params$zfd_samples = as.data.frame(do.call(cbind, zfd_samples[1:2])[
                            sample(1:length(zfd_samples[['b_zpool']]), 1000), ])
colnames(params$zfd_samples) = c("zfd_intercept", "zfd_slope")

# Model 5: Density-dependent dynamic zoospore pool
zdd_samples = extract(pomp_mod_dd)
params$zdd_intercept = mean(zdd_samples[['b_zpool']])
params$zdd_slope = mean(zdd_samples[['b_density']])
params$zdd_samples = as.data.frame(do.call(cbind, zdd_samples[1:2])[
                            sample(1:length(zdd_samples[['b_zpool']]), 1000), ])
colnames(params$zdd_samples) = c("zdd_intercept", "zdd_slope")

# Model 6: Only infection from the dynamic zoospore pool
zconst_samples = extract(pomp_mod_zp)
params$zconst_intercept = mean(zconst_samples[['b_zpool']])
params$zconst_samples = as.data.frame(zconst_samples[[1]][
                            sample(1:length(zconst_samples[['b_zpool']]), 1000)])
colnames(params$zconst_samples) = c("zconst_intercept")

## Additional parameters ##

# Survival of uninfected individuals from Briggs et al. 2005 converted to
# a 3 day time scale.
params$class_zero_surv = 0.9^(1 / 122) 

# Default probability of infection parameters from Wilber et al. 2016 that 
# depend on temperature but not host density.  These are overwritten during the
# simulation.
params$prob_inf_int = -1.660243
params$prob_inf_temp = 0.1018972

# Save the parameters and corresponding samples from the posterior distributions
saveRDS(params, "../../results/bayesian_parameter_estimates.rds")





