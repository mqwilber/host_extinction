## Description
## -----------
## Contains some helpful functions for analyzing an IPM model in R. See below
## for the various functions

## Author: Mark Wilber

##############################################################################

library(MASS)
library(ggplot2)
library(data.table)
library(parallel)

local = TRUE # TRUE if you are not running this on the cluster. FALSE if you are

if(local){
  library(reshape2)
  library(RColorBrewer)
  library(fields)
  library(popbio)
  library(grid)
  library(gridExtra)
}


set_discretized_values = function(min_size, max_size, bins){
  # Calculates the necessary parameters to use the midpoint rule to evaluate
  # the IPM model

  # Parameters
  # ----------
  # min_size : The lower bound of the integral
  # max_size : The upper bound of the integral
  # bins : The number of bins in the discretized matrix

  # Returns
  # -------
  # list
  # min_size, max_size, bins, bnd (edges of discretized kernel), y (midpoints),
  # h (width of cells)


  # Set the edges of the discretized kernel
  bnd = min_size+c(0:bins)*(max_size-min_size) / bins

  # Set the midpoints of the discretizing kernel. Using midpoint rule for evaluation
  y = 0.5 * (bnd[1:bins] + bnd[2:(bins + 1)])

  # Width of cells
  h = y[2] - y[1]

  return(list(min_size=min_size, max_size=max_size, bins=bins, bnd=bnd, y=y,
                h=h))

}


get_the_kernel = function(g_xpx, s_x, r_x, bins, y, params, h, 
                            elasg=1, elass=1, elasr=1){
  # Function to get the IPM kernel from the vital rate functions.  
  # Doesn't include the 0 class

  # Parameters
  # ----------

  # g_xpx : the growth function
  # s_x : the survival function
  # r_x : the loss function
  # bins : number of bins in the discretize model
  # y : midpoints (length(y) = bins)
  # params : list of relevant parameters to pass into vital functions
  # h : width of discretized cell. Used for midpoint evaluation
  # elas*: An elasticity argument to pass to each of the vital rate functions
  #         If 1, the default fitted value is returned for the vital rate
  #         function. If different than 1 (i.e. 1.0001) the function is
  #         perturbed by this amount. 
  #         g = growth function, 
  #         s = survival function
  #         r = recovery/loss function 

  # Returns
  # -------
  # : list
  #   P (kernel), G (growth Matrix), R (loss vector), S (survival vector)

  # Make growth matrix
  G = h * outer(y, y, g_xpx, params=params, elas=elasg)

  # Make survival matrix
  S = s_x(y, params=params, elas=elass)

  # Make loss matrix
  R = r_x(y, params=params, elas=elasr)

  # Make the full kernel
  #P = G
  # Faster without a for loop...
  #for(i in 1:bins) P[,i]=G[,i]* S[i] * (1 - R[i])

  P = G %*% diag(S * (1 - R))

  return(list(P=P, S=S, R=R, G=G))

}

get_full_P = function(P, row1, col1, num_add=1, add_top=TRUE){
  # Function takes a P matrix and sticks row1 and col1 onto it.
  #
  # Good for adding extra stages in the IPM framework

  # Parameters
  # ----------
  # P : the transition matrix
  # row1 : vector specifying the probability of any stage transition in P
  #        to the new stage. The length of this vector should be nrow(P) + 1
  # col1 : vector specifying the probability of transitioning from new stage
  #        to old stage. The length of this vector should be ncol(P) + 1
  # num_add : The number of rows and columns you are adding to the matrix
  # add_top : By default sticks them on the top of the matrix. FALSE adds them
  #           to the bottom.

  # Returns
  # -------
  # : matrix
  #   The transition matrix with a 0 class

  if(add_top){ # Add extra columns and rows to the top

    full_P = matrix(NA, nrow=nrow(P) + num_add, ncol=ncol(P) + num_add)
    full_P[1:num_add, ] = row1
    full_P[ , 1:num_add] = col1
    full_P[(num_add + 1):nrow(full_P), (num_add + 1):ncol(full_P)] = P

  } else{ # Add extra rows and columns to the bottom

    tadd = num_add - 1
    full_P = matrix(NA, nrow=nrow(P) + num_add, ncol=ncol(P) + num_add)
    full_P[(nrow(full_P) - tadd):nrow(full_P), ] = row1
    full_P[ ,(ncol(full_P) - tadd):ncol(full_P)] = col1
    full_P[1:(nrow(full_P) - num_add), 1:(ncol(full_P) - num_add)] = P

  }

  return(full_P)

}

absorption_times = function(P){

  # Calculate the mean time to absoprtion (death) given a transition matrix P
  # Calculate the variance as well

  # Parameters
  # ----------
  # P : transition matrix
  #
  # Returns
  # -------
  # list: mean=mean absorption time, var=var absorption time

  I = diag(dim(P)[1])

  # Calculate the fundamental matrix
  N = solve((I - P))

  # Exepected time to absorption starting in transient state x
  mean_absorb_time = colSums(N)
  var_absorb_time = colSums((2*(N %*% N) - N)) - mean_absorb_time^2

  return(list(mean=mean_absorb_time, var=var_absorb_time))

}

get_elasticity_and_sensitivity = function(P, h, y, plot_it=T, save="temp"){
  # Calculate the elasticities and sensitivities of the population growth rate
  # to transition probabilities in the matrix P

  # Parameters
  # ----------
  # P : transition matrix
  # h : cell width
  # y : midpoints of cells, more plotting

  # Returns
  # -------
  # list: elas= elasticity matrix, sens= sensitivity matrix


  lam=Re(eigen(P)$values[1])
  w.eigen = Re(eigen(P)$vectors[,1])
  stable.dist = w.eigen / sum(w.eigen)
  v.eigen = Re(eigen(t(P))$vectors[,1])
  repro.val = v.eigen / v.eigen[1]

  v.dot.w = sum(stable.dist*repro.val) * h

  sens = sensitivity(P)# outer(repro.val,stable.dist) / v.dot.w
  elas = elasticity(P)# matrix(as.vector(sens)*as.vector(P) / lam, nrow=length(y))

  if(plot_it){

    myPalette = colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

    long_elas = format_for_heat(t(elas), y, y, c('size', 'sizeNext', 'value'))
    long_elas$id = "Elasticity"

    # Just plotting elasticity
    zplot = ggplot(long_elas, aes(x=size, y=sizeNext, fill=value)) + geom_tile()
    zplot = zplot + scale_fill_gradientn(colours = myPalette(100), name="Elasticity")
    zplot = zplot + theme_bw() + xlab("Log zoospore load at t") + ylab("Log zoospore load at t + 1")
    zplot = zplot + facet_wrap(~id)

    ggsave(paste(save, ".pdf", sep=""), height=5, width=7)

  }

  return(list(elas=elas, sens=sens))

}

get_stable_dist = function(P){

  # Get the stable distribution of a P matrix based on Caswell 2001

  w.eigen = Re(eigen(P)$vectors[,1])
  stable_dist = w.eigen / sum(w.eigen)
  return(stable_dist)

}

get_vm_ratio = function(P, y){

  # Get the variance to mean ratio of the stable distribution given the full
  # transition matrix P. y is the step size

  stable_dist = get_stable_dist(P)

  cond_stable_dist = stable_dist[2:length(stable_dist)] / (1 - stable_dist[1])
  cond_mean = sum(cond_stable_dist * y)
  cond_var = sum(cond_stable_dist * y^2) - (cond_mean)^2

  return(cond_var / cond_mean)

}

##############################################################################

g0_x = function(x, params, elas=1){
  # Define the initial infection function
  # elas defines function level elasticity
  return(dnorm(x, mean=params$clump_mean, sd=params$clump_sd) * elas)

}

s_x = function(x, params, elas=1) {
  # Define survival function
  # elas defines function level elasticity

  # Survival function is scaled for better estimation
  x_scale = (x - params$surv_mean) / params$surv_sd

  u = exp(params$surv_int + params$surv_size * x_scale)
  return((u / (1 + u)) * elas)

}


r_x = function(x, params, elas=1) {
  # Define the loss function, r(x) (also defined as l(x))
  # elas defines function level elasticity

  u = exp(params$loss_int + params$loss_size*x + params$loss_temp)
  probs = u / (1 + u)

  return(probs * elas)

}


k_xpx = function(xp, x, params){

  # Define the full kernel (for eviction calculations)

  return(g_xpx(xp, x, params) * s_x(x, params) * (1 - r_x(x, params)))

}


g_xpx = function(xp, x, params, elas=1){
  # Tranistion from x to xp given some parameters (growth function)
  # elas defines function level elasticity

  if(params$linear){ # Linear growth function

    return(dnorm(xp, mean=params$growth_int + params$growth_size * x +
       params$growth_temp,
       sd=sqrt(params$growth_sigma2*exp(2 * params$growth_sigma2_exp * x))) 
                          * elas)

  } else{ # Non-linear growth function

    x[x <= params$growth_deriv0] = params$growth_deriv0
    return(dnorm(xp, mean=params$growth_int + params$growth_size * x +
       params$growth_temp + params$growth_size_sq * x**2,
       sd=sqrt(params$growth_sigma2*exp(2 * params$growth_sigma2_exp * x))) * elas)
  }
}

get_zoospore_surv_prob = function(temp, elas=1, zsurv=NULL){
  # Uses a smoothing spline to predict the zoospore survival probability for
  # any temp between 4C and 28C using data from Woodhams et al. 2008, Ecology
  #
  # Parameters
  # ----------
  # temp: float, temperature between 4 and 28
  # elas : float, function-level elasticity. Used to perturb the function.
  #         Include by setting a little bit bigger or 
  #         less than 1. i.e. 1.0001 or .99999
  # zsurv : float, user-defined default survival probability
  #
  # Returns
  # -------
  # : zoospore survival probability

  # Temperature and death rate values from Table 1 in Woodhams et al. 2008
  temp_vals = c(4.0, 14.5, 23.0, 28.0)
  drates = c(0.0027, 0.0035, 0.0160, 0.0410)

  fit = smooth.spline(temp_vals, drates)
  pred_vals = predict(fit, temp)$y

  # Convert to right time scale. Rates are on an hourly time scale and we want
  # a 3 day survival probability.
  surv_prob = exp(-pred_vals * 24 * 3)

  #return(0.2227553)

  if(is.null(zsurv)){
    return(surv_prob * elas)
  } else{ # Set as constant probability if zsurv is not NULL
    return(zsurv * elas)
  }

}

set_temp_params = function(desired_temp, params, linear, elasz=1, zsurv=NULL){
  # Specific parameters for each temperature, given a desired temp
  #
  # Parameters
  # ----------
  # desired_temp: int, the desired temperature at which to calculate params
  # params: list, a list-like object with the base params
  # linear: bool, TRUE is you are using a linear growth function. FALSE otherwise
  # elasz: float, the elasticity for the zoospore survival probability.
  #         Include by setting a little bit bigger or 
  #         less than 1. i.e. 1.0001 or .99999
  # zsurv : float, user-defined default survival probability
  #
  # Returns
  # --------
  # : Updated parameters

  tparams = data.frame(holder=NA)

  # Parameters that are the same for everybody
  tparams$surv_int = params$surv_int
  tparams$surv_size = params$surv_slope
  tparams$surv_mean = params$surv_mean
  tparams$surv_sd = params$surv_sd
  tparams$zspore_surv = get_zoospore_surv_prob(desired_temp, elas=elasz, 
                                                zsurv=zsurv)

  if(linear){ # Assuming a linear effect of size

    # 5 parameters for growth function
    tparams$growth_int = params$growth_int
    tparams$growth_temp = params$growth_temp * desired_temp
    tparams$growth_size = params$growth_size
    tparams$growth_sigma2 = params$growth_sigma2
    tparams$growth_sigma2_exp = params$growth_sigma2_exp

    # Three parameters for the loss function
    tparams$loss_int = params$loss_int
    tparams$loss_size = params$loss_size
    tparams$loss_temp = params$loss_temp * desired_temp

    # Probability of becoming infected
    u = exp(params$prob_inf_int + params$prob_inf_temp * desired_temp)
    tparams$prob_inf = u / (1 + u)

    # Clumped distribution
    tparams$clump_mean = params$clump_int + params$clump_temp * desired_temp
    tparams$clump_sd = sqrt(params$clump_sigma2*exp(2 * params$clump_sigma2_exp * desired_temp))

    # Class zero survival
    tparams$class_zero_surv = params$class_zero_surv
    tparams$linear = TRUE

  } else{ # Assuming a non-linear effect of size

      # Parameters for growth function
    tparams$growth_int = params$growth_int
    tparams$growth_temp = params$growth_temp * desired_temp
    tparams$growth_size = params$growth_size + params$growth_size_temp_int * desired_temp
    tparams$growth_size_sq = params$growth_size_sq
    tparams$growth_sigma2 = params$growth_sigma2
    tparams$growth_sigma2_exp = params$growth_sigma2_exp

    # Find where the derivative is 0
    tparams$growth_deriv0 = (-1 * tparams$growth_size) / (2 * tparams$growth_size_sq)

    # Three parameters for the loss function
    tparams$loss_int = params$loss_int
    tparams$loss_size = params$loss_size
    tparams$loss_temp = params$loss_temp * desired_temp

    # Probability of becoming infected
    u = exp(params$prob_inf_int + params$prob_inf_temp * desired_temp)
    tparams$prob_inf = u / (1 + u)

    # Clumped distribution
    tparams$clump_mean = params$clump_int + params$clump_temp * desired_temp
    tparams$clump_sd = sqrt(params$clump_sigma2*exp(2 * params$clump_sigma2_exp * desired_temp))

    # Class zero survival
    tparams$class_zero_surv = params$class_zero_surv
    tparams$linear = FALSE

  }

  return(tparams)
}


## Various transmission functions based on both model fits and
## what we think is happening in the system

get_infection_probability = function(pother, pfull, method, elas=1){
  # Generic function to return an infection probability given a specified type
  # of transmission function given by the option method
  #
  # Parameters
  # ----------
  # pother : parameter vector
  #         Must contain information regarding density of adults, zoospores,
  #         or proportion infected depending on the transmission function
  #         used.
  # pfull : parameter list
  #         Contains the parameters for the vital rate functions.
  # method : string, either
  #          "density_dependent_individuals"
  #          "frequency_dependent_individuals" 
  #          "zoospore_pool_infection_dd"
  #          "zoospore_pool_infection_fd"
  #          "zoospore_pool_infection_const"
  # elas : float, the elasticity. Include by setting a little 
  #         less than 1. i.e. .99999. 
  #
  # Returns
  # : float, the probability of infection in a time step

  # Chose transmission method
  if(method == "density_dependent_individuals"){

    prob_inf = density_dependent_individuals(pother, pfull)

  } else if(method == "frequency_dependent_individuals"){

    prob_inf = frequency_dependent_individuals(pother, pfull)

  } else if(method == "zoospore_pool_infection_dd") {

    prob_inf = zoospore_pool_infection(pother, pfull, "dd")

  } else if(method == "zoospore_pool_infection_fd"){

    prob_inf = zoospore_pool_infection(pother, pfull, "fd")

  } else if(method == "zoospore_pool_infection_const"){

    prob_inf = zoospore_pool_infection(pother, pfull, "const")

  }else{
    stop(paste(method, "is not a valid transmission function"))
  }

  return(prob_inf * elas)

}

density_dependent_individuals = function(pother, pfull){
  # Density dependent transmission defined by
  #     prob_inf = 1 - exp(-(b1 + b2 * number of infected adults))
  # Constant infection from the environment
  #
  # Returns
  # -------
  # Infection probability for a given time step that depends on the number
  # of infected adults in the previous time step

  beta1 = pfull$dd_intercept*(pother$sim$step_length / pother$ipm$time_convert)
  beta2 = pfull$dd_slope

  # Get the number of infected adults 
  num_inf = pother$ipm$num_inf

  # Multiply by a time step
  num_inf_norm = num_inf * (pother$sim$step_length / pother$ipm$time_convert)
  zs_present = as.integer(pother$ipm$zoo_inf > 0) # Only if zoospores are present

  infection_prob = 1 - exp(-(beta1*zs_present + beta2*num_inf_norm))

  return(infection_prob)

}

frequency_dependent_individuals = function(pother, pfull){
  # Frequency dependent transmission defined by
  #     prob_inf = 1 - exp(-(b1 + b2 * proportion infected adults))
  # Constant infection from the environment
  #
  # Returns
  # -------
  # Infection probability for a given time step that depends on the proportion
  # of infected adults in the previous time step

  beta1 = pfull$fd_intercept*(pother$sim$step_length / pother$ipm$time_convert)
  beta2 = pfull$fd_slope

  # Proportion of infected adults
  prop_inf = pother$ipm$prop_inf 

  # multiply by a time step
  prop_inf_norm = prop_inf * (pother$sim$step_length / pother$ipm$time_convert)
  zs_present = as.integer(pother$ipm$zoo_inf > 0) # Only if zoospores are present

  infection_prob = 1 - exp(-(beta1*zs_present + beta2*prop_inf_norm))

  return(infection_prob)

}

zoospore_pool_infection = function(pother, pfull, method){
  # Function specifying the probability of infection given the total
  # number zoospores in the population
  #       prob_inf = 1 - exp(-(b1*ZE + b2*density/frequency of infecteds))
  #
  # Parameters
  # ----------
  # pother, pfull: list/data.frames holding relevant parameters
  # method: 
  #      "dd" for density-dependent
  #      "fd" for frequency-dependent
  #      "const" for only infection from zoospore pool

  if(method == "dd"){ # Density-dependent transmission

    beta1 = pfull$zdd_intercept
    beta2 = pfull$zdd_slope

    num_inf = pother$ipm$num_inf # Set to density density dependent

  } else if(method == "fd"){

    beta1 = pfull$zfd_intercept
    beta2 = pfull$zfd_slope

    num_inf = pother$ipm$prop_inf # Set to frequency-dependent

  } else if(method == "const"){ # Set to only zoospore pool dependent

    beta1 = pfull$zconst_intercept
    beta2 = 0
    num_inf = 0
  }

  time_fact = (pother$sim$step_length / pother$ipm$time_convert)
  num_inf_norm = num_inf * time_fact

  # Multiply by a time step
  total_z_norm = log(pother$ipm$zoo_inf + 1) * time_fact

  infection_prob = 1 - exp(-(beta1*total_z_norm + beta2*num_inf_norm))

  return(infection_prob)
}


##################### Multiseason simulation functions ########################

get_temp = function(val, start_temp){
  # Uses a cosine function to get a temperature value
  # 
  # Parameters
  # ----------
  # val : float, time step in 122 day cycle
  # start_temp: The starting temperature of the cosine function
  # 
  # Returns
  # -------
  # : temperature at time point
  return(12 + -8 * cos((val* 2*pi / 122) + acos((start_temp - 12) / -8)))
}

logit_fxn = function(x, a, b){
  # Logistic function
  return(exp(a + b*x) / (1 + exp(a + b*x)))
}

get_adult_density = function(all_results, index, num_tads){
  # Function gets the total adult densities in all IPMs in all_results.
  # Used in the simulation to store adults density at each time step.
  # 
  # Parameters
  # ----------
  # all_results : list of matrices
  # index : index of column to extract from each matrix
  # num_tads : number of stages to exclude

  adult_sum = function(x, i, trun) sum(x[-c(1:trun, nrow(x)), i])

  return(sum(unlist(lapply(all_results, adult_sum, index, num_tads))))
}

get_zoospore_density = function(all_results, index, num_exclude){
  # Function gets the total zoospore densities in all IPMs in all_results
  # Used in the simulation to store zoospore density at each time step.
  # Parameters
  # ----------
  # all_results : list of matrices
  # index : index of column to extract from each matrix
  # num_exclude : 1:number of stages to exclude
  # 
  zoo_sum = function(x, i, trun) sum(x[-c(1:trun), i])

  return(sum(unlist(lapply(all_results, zoo_sum, index, num_exclude))))
}

assign_initial_pop_values = function(init_val, mean_val, sd_val, resist_params){
  # Assign initial adult abundances to various trait classes assuming that
  # the distribution of the trait is normally distributed with some
  # mean and variance
  #
  # Parameters
  # ----------
  # init_val : Initial number of adults
  # mean_val : mean heterogeneity value
  # sd_val : standard deviation in heterogeneity value
  # resist_params : list with resistance parameters

  if(length(resist_params$y) != 1){ # May only have a single resistant class

    probs = resist_params$h * dnorm(resist_params$y, mean=mean_val, 
                          sd=sd_val)
    norm_probs = probs / sum(probs)
    init_vect = round(init_val * norm_probs)

  } else {
    init_vect = init_val
  }

  return(init_vect)
}

get_resistance_params = function(mean_val, sd_val, resist_bins, quant=c(0.05, 0.95)){
  # Calculates the binned resistance parameters. Deals with the fact that 
  # variance could be 0 and only included 1 resistance class if that is the
  # case. 
  #
  # Parameters
  # ----------
  # mean_val : float, mean of resistance parameter
  # sd_val : float, standard deviation of resistance parameter
  # resist_bins : number of bins used to discretize the IPM
  # quant : tuple of lower percent and upper percent. 
  #         Assumes the trait is normally distributed and sets the lower and
  #         upper bounds based on a the bounds in a normal distribution 
  #
  # Returns
  # -------
  # : list-like object that contains the midpoints of the discretized IPM

  if(sd_val != 0){

    lower = qnorm(quant[1], mean=mean_val, sd=sd_val)
    upper = qnorm(quant[2], mean=mean_val, sd=sd_val)
    resist_params = set_discretized_values(lower, upper, resist_bins)

  } else{
    resist_params = list()
    resist_params$y = c(mean_val) # Change in growth slope
    resist_params$h = 1 
  }

  return(resist_params)

}


build_demography_matrix = function(pother, slope, y, stages, adults, 
                                        dropZ=FALSE){
  # Building a demography matrix for between year transitions. Includes
  # 
  #
  # Parameters
  # -----------
  # pother : list of parameters from yaml file
  # slope : Slope of growth function, resistance parameter
  # y : mid point of infection loads for infected adults. Necessary to compute
  #     decline in reproduction with infection
  # stages: total number of stages in the demography matrix
  # adults: number of adults in the population for density dependent recruitment
  #         Set to zero if you don't want any density dependence
  # dropZ: bool, TRUE to drop zoospore class, FALSE otherwise
  #
  # Returns
  # -------
  # : A demography matrix used to project populations stage structure between
  #   years

  num_tads = pother$ssmod$tad_stages
  demo_P = array(0, dim=c(stages, stages))
  lower_ident = diag(pother$sim$bins + 2)

  # Set the null transitions
  demo_P[(num_tads + 1):nrow(demo_P),
         (num_tads + 1):ncol(demo_P)] = lower_ident

  # Set the tadpole transitions. 
  # tad_diag defines rows t2, t3, adults and cols t1, t2, t3
  tad_diag = diag(num_tads)
  diag(tad_diag) = pother$ssmod$tad_transitions * pother$ssmod$tad_not_adult

  # Tads can transition to adults. See Briggs et al. 2005
  tad_diag[num_tads, ] = ((1 - pother$ssmod$tad_not_adult) * # Transition to adult
                          pother$ssmod$tad_transitions * # Survive
                          pother$ssmod$tad_meta_surv * # Survive metamorphosis
                          exp(-adults / pother$ssmod$K)) # Ricker recruitment
  demo_P[2:(num_tads + 1), 1:num_tads] = tad_diag


  # Add fecundity. Any infection less than log(1) = 0 set to 0
  y_val = c(0, y) # include uninfected adults as 0
  y_val[y_val <= 0] = 0

  # Exponential decline in reproduction with infection and 
  # Reduce the reproduction probability given the slope of the resistance
  # function. THIS IS WHERE THE RESISTANCE TRADEOFF COMES INTO PLAY 
  repro_probs = with(pother$ssmod,
          uninfected_adult_repro_prob * exp(-1*repro_decline*y_val))

          #* logit_fxn(slope, resistance_trade_off[1], resistance_trade_off[2]))

  # Mean reproduction
  repro_means = pother$ssmod$repro * repro_probs

  # Minus one to avoid the zoospore column. Plus one to start at adults
  demo_P[1, (num_tads + 1):(ncol(demo_P) - 1)] = repro_means

  if(!dropZ){
    return(demo_P) # Don't drop the zoospores
  } else{
    return(demo_P[-nrow(demo_P), -ncol(demo_P)])
  }

}

demography_step = function(n_t, demo_matrix, pother){
  # Perform a between year step in demography matrix
  #
  # Parameters
  # ----------
  # n_t : population vector at time t
  # demo_matrix : demography matrix calculated from build demography_matrix
  # pother : parameters from yaml file
  #
  # Returns
  # -------
  # : The vector at the next time step give n_t and the projection matrix 
  #  demo_matrix

  # Add the mean reproductive output values
  demo_matrix[1, ] = demo_matrix[1, ]*pother$ssmod$repro # rpois(1, pother$ssmod$repro)

  n_t = as.matrix(n_t, nrow=length(n_t), ncol=1)

  demo_step = demo_matrix %*% n_t
  return(demo_step)

}


within_season_simulation = function(n_0, params, pother, pfull, bins, bnd, y, h,
                      kernels, transmission_type, repro_now, 
                      adult_density, slope,
                      stochastic=FALSE, elasphi=1, elasg0=1){
  # Calculates a one step projection given a vector n_0
  #
  # Parameters
  # ----------
  # n_0 : initial population vector
  # params: parameters for the IPM fit from data
  # pother: Parameters stored in yaml file
  # pfull : Parameters from IPM fit from data before conversion
  # bins, bnd, y, h : discretization parameters
  # kernels : density independent IPM kernels
  # transmission_type : a character vector/string specifying the type of
  #                     transmission function to use
  # repro_now: Boolean determining whether reproduction occurs in this time step
  # adult_density: Density of adults at current time step
  # slope : Slope of the Bd growth function
  # stochastic: If TRUE, implements a stochastic version of the IPM
  # elasphi : float, elasticity perturbation of infection function
  # elasg0 : float, elasticity perturbation of initial infection function

  # Returns
  # -------
  # list of updated pop vector and infection probability

  infection_prob = get_infection_probability(pother, pfull, 
                                    transmission_type, elas=elasphi)

  # The zero -> zero transition
  zero_trans = params$class_zero_surv * (1 - infection_prob)

  # The zero to nonzero transitions: Transitioning from zero to infected
  prob_clump = h * g0_x(y, params, elas=elasg0) #dnorm(y, mean=params$clump_mean, sd=params$clump_sd)
  ztnz = params$class_zero_surv * infection_prob * prob_clump
  col1 = c(zero_trans, ztnz)

  # Compute the infected to uninfected row. This is just the probability of
  # surviving with load $x$ ($s(x)$) and then losing an infection $r(x)$.  We
  # just need to multiply $s(x) r(x)$

  # The loss probabilities: Transitioning from non-zero to zero
  nztz = kernels$S * kernels$R
  row1 = c(zero_trans, nztz)

  # Add the infection probabilities onto the P matrix
  full_P = get_full_P(kernels$P, row1, col1)

  # Add the zoospore class onto the p matrix
  zrow = c(0, exp(y)*pother$ipm$zspore_prod, params$zspore_surv)
  zcol = c(rep(0, length(y) + 1), params$zspore_surv)

  full_P = get_full_P(full_P, zrow, zcol, add_top=FALSE)

  # Add the three tadpole classes onto the matrix
  # Assuming tadpoles are immediately infected and have constant load

  full_tad = diag(nrow(full_P) + pother$ssmod$tad_stages)
  # Set tadpole survival probability for each 3 day time step.
  diag(full_tad)[1:pother$ssmod$tad_stages] = pother$ipm$tad_survival
  tad_rows = full_tad[1:pother$ssmod$tad_stages, ]
  tad_cols = full_tad[, 1:pother$ssmod$tad_stages]
  tad_cols[nrow(tad_cols), ] = pother$ipm$mean_tad_load

  # Add on tadpole columns
  full_P = get_full_P(full_P, tad_rows, tad_cols,
                  num_add=pother$ssmod$tad_stages, add_top=TRUE)

  # Determine if reproduction needs to occur at this step
  if(repro_now){
    demo_matrix = build_demography_matrix(pother, slope, y, 
                                    nrow(full_P), adult_density, dropZ=FALSE)
    full_P = demo_matrix %*% full_P # ORDER MATTERS! 
  }

  if(!stochastic){
    n_t_new = nonstochastic_step(full_P, n_0, params, pother) #full_P %*% n_0
  } else{
    n_t_new = stochastic_step(full_P, n_0, params, pother)
  }

  return(list(nt=n_t_new, inf_prob=infection_prob, full_P=full_P))

}


multiseason_simulation = function(pfull, pother, standing_var, resist_bins,
                  hetero_param="growth_size",
                  elas_vals=list(elasphi=1,elasg0=1,elasg=1,elass=1,elasr=1,elasz=1), 
                  sims=50, output=TRUE, zsurv=NULL){
    # A multiseason simulation of the seasonal host-parasite IPM
    #
    # Parameters
    # ----------
    # pfull : list-like object with temperature-dependent IPM parameters
    # pother : list-like objects with IPM and simulation parameters specified 
    #            in model_parameters.yml
    # sims : int, number of simulations to run. Set to one if the simulation
    #              is not stochastic
    # standing_var : Float, the amount of standing variation in resistance.
    #                If 0, resist_bins is set to 1. This is what you want 
    #                if you don't want heterogeneity
    # resist_bins : int, the number of bins used to approximate the continuous
    #                resistance trait.  Approximating the continuous with the
    #                mid-point rule. Set to 1 if you don't want heterogeneity.
    # elas_vals : list of floats, elasticity perturbation of the vital rate 
    #          functions
    #         elasphi : infection/transmission function 
    #         elasg0 : initial infection function
    #         elasg : growth function
    #         elass : survival function
    #         elasr : recovery function
    #         elasz : zoospore survival probability
    # hetero_param : str, the parameter in pfull that shows heterogeneity between
    #                hosts.  By default this is growth_size.
    # zsurv : NULL or float, If null uses empirically estimated zoospore decay
    #           probability. If float, expects a survival probability between
    #           0 and 1 that specifies the zoospore survival prob over a given
    #           time step (3 days in this model)
    # Returns
    # --------
    # : list
    #    The length of the list is the same length as sims.  Each item in 
    #    the list has a length of resist bins where each input is a separate
    #    IPM run with a different resistance value.  If standing var is 0
    #    there can only be one resist bin. 

    # If there is no variance there is only one resist bin
    if(standing_var == 0){
        resist_bins = 1
    }

    # Unpack elasticity vals. TODO: Use list2env()
    elasphi = elas_vals$elasphi; elasg0 = elas_vals$elasg0;
    elasg = elas_vals$elasg; elass = elas_vals$elass; elasr = elas_vals$elasr
    elasz = elas_vals$elasz

    # Set the basic IPM matrix parameters
    steps = with(pother$sim, (days_per_year * years) / step_length)
    matrix_params = with(pother$sim,
                            set_discretized_values(min_size, max_size, bins))
    bnd = matrix_params$bnd
    y = matrix_params$y
    h = matrix_params$h
    bins = pother$sim$bins

    all_temps = array(NA, dim=steps)
    all_temps[1] = get_temp(0, pother$sim$start_temp)

    # Adding 2 to the dimensions to account for uninfected adults and zoospores
    results = with(pother$ssmod,
                array(NA, dim=c(bins + tad_stages + juv_stages + 2, steps + 1)))

    # Initialize vector from yaml file
    #   1. First items are tadpole stages and adult
    #   2. Second - bins items are infected classes
    #   3. The last item initial number of zoospores in the pool
    n_0 = with(pother$ssmod, c(init_tad_vals, init_adult_val, rep(0, bins),
                                        init_zspore_val))
    results[, 1] = n_0
    repro_correct_temp = pother$sim$repro_first

    pfull$ranef_var = standing_var # Set standing variation

    # Set resistance bins
    resist_params = get_resistance_params(pfull[[hetero_param]], 
                                          sqrt(standing_var), resist_bins)

    init_vect = assign_initial_pop_values(pother$ssmod$init_adult_val, 
                          pfull[[hetero_param]],
                          sqrt(standing_var),
                          resist_params)

    init_full = sapply(c(pother$ssmod$init_tad_vals, pother$ssmod$init_adult_val),
                              assign_initial_pop_values, 
                              pfull[[hetero_param]], sqrt(standing_var),
                              resist_params)

    if(resist_bins == 1){
      init_full = t(as.matrix(init_full))
    }

    sim_results = list()

    ########################## BEGIN SIMULATION ################################

    for(sim in 1:sims){ # Loop through number of simulations

        if(output){
          print(paste("Beginning simulation", sim))
        }

        # Reproduction requires this to be true and temperature to be right
        # Ensures that there is only one reproduction per time cycle
        repro = rep(pother$sim$repro_first, resist_bins)

        all_results = replicate(resist_bins, results, simplify=FALSE)

        # Assign the right initial abundances to the given classes.
        for(j in 1:length(all_results)){
          all_results[[j]][1:(pother$ssmod$tad_stages + pother$ssmod$juv_stages + 1), 1] = init_full[j, ]
        }

        for(i in 2:(steps + 1)){ # Looping through time steps

            # Calculate the temperature at the current time
            all_temps[i] = get_temp(i-1, pother$sim$start_temp)

            # Get density across all result at time i - 1
            adult_density = get_adult_density(all_results, i-1, 
                                            pother$ssmod$tad_stages)
            infected_density = get_adult_density(all_results, i-1, 
                                            pother$ssmod$tad_stages + 1)
            zoospore_density = get_zoospore_density(all_results, i-1, 
                                            pother$ssmod$tad_stages + 1 + bins)
            
            pother$ipm$num_inf = infected_density
            pother$ipm$prop_inf = ifelse(adult_density != 0, 
                        infected_density / adult_density, 0) # Can't divide by zero...population is extinct
            pother$ipm$zoo_inf = zoospore_density

            for(j in 1:resist_bins){ # Need to run a separate IPM for each hetero parameter

                # Update the slope parameter
                pfull[[hetero_param]] = resist_params$y[j]
                slope = pfull[[hetero_param]]

                # Calculate the temperature dependent kernel
                params = set_temp_params(all_temps[i], pfull, 
                                          pother$ipm$linear_temp, elasz=elasz,
                                          zsurv=zsurv)

                # TODO: Potential speed up. Think about when you need to
                # update these kernels
                kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h,
                                        elasg=elasg, elass=elass, elasr=elasr)

                # Perform demographic change if time is right
                if(round(all_temps[i], 2) == pother$sim$temp_at_repro){

                    if(repro[j]){
                        repro_now = TRUE
                        repro[j] = !repro[j]
                    }
                    else{
                        repro_now = FALSE
                        repro[j] = !repro[j] # prevents double reproduction at 
                                        # temps that happen twice a year
                    }
                } else{ 
                    repro_now = FALSE 
                }

                # Perform the IPM simulation
                n_t = as.matrix(all_results[[j]][, i - 1], 
                                    nrow=nrow(all_results[[j]]), ncol=1)

                one_step_out = within_season_simulation(n_t, params, pother, 
                                    pfull, bins, bnd, y, h, kernels, 
                                    pother$ipm$transmission_type, repro_now,
                                    adult_density, slope,
                                    stochastic=pother$sim$stochastic, 
                                    elasphi=elasphi, elasg0=elasg0)

                all_results[[j]][, i] =  one_step_out$nt

            } # End resistance slope loop

        } # End season loop

        sim_results[[sim]] = all_results

    } # End simulation loop

    return(sim_results)

}


extinction_prob_daily = function(pfull, pother, standing_var, resist_bins,
                hetero_param ="growth_size", 
                elas_vals=list(elasphi=1,elasg0=1,elasg=1,elass=1,elasr=1,elasz=1)){
  # Calculates extinction probabilities on a three day time scale using an
  # adaptive dynamics approach given in Schreiber and Ross 2016. Assumes
  # density-independent recruitment and transmission in order to use a 
  # branching process assumption.
  #
  # Parameters
  # ----------
  # pfull : list-like object with temperature-dependent IPM parameters
  # pother : list-like objects with IPM and simulation parameters specified 
  #            in model_parameters.yml
  # standing_var : Float, the amount of standing variation in resistance.
  #                If 0, resist_bins is set to 1.
  # resist_bins : int, the number of bins used to approximate the continuous
  #                resistance trait.  Approximating the continuous with the
  #                mid-point rule
  # hetero_param : str, the parameter in pfull that shows heterogeneity between
  #                hosts.  By default this is growth_size.
  # elas_vals : list of floats, elasticity perturbation of the vital rate 
  #          functions
  #         elasphi : infection/transmission function 
  #         elasg0 : initial infection function
  #         elasg : growth function
  #         elass : survival function
  #         elasr : recovery function
  #         elasz : zoospore decay 
  #
  # Returns
  # -------
  # : extinction probability at every three-day time step

  # If there is no variability there is only one resist bin
  if(standing_var == 0){
    resist_bins = 1
  }

  steps = with(pother$sim, (days_per_year * years) / step_length)
  vals = seq(1, steps, by=122)
  #vals = c(10)

  # Set the basic IPM matrix parameters
  matrix_params = with(pother$sim,
                          set_discretized_values(min_size, max_size, bins))
  bnd = matrix_params$bnd
  y = matrix_params$y
  h = matrix_params$h
  bins = pother$sim$bins

  ############### Initialize data vectors for simulation ####################

  #pfull$growth_ranef_var = standing_var
  pfull$ranef_var = standing_var
  resist_params = get_resistance_params(pfull[[hetero_param]], 
                                        sqrt(standing_var), resist_bins)
  # init_vect = assign_initial_pop_values(pother$ssmod$init_adult_val,
  #                         pfull[[hetero_param]],
  #                         sqrt(standing_var),
  #                         resist_params)

  init_full = sapply(c(pother$ssmod$init_tad_vals, pother$ssmod$init_adult_val), 
                          assign_initial_pop_values, 
                          pfull[[hetero_param]], sqrt(standing_var),
                          resist_params)

  if(resist_bins == 1){
    init_full = t(as.matrix(init_full))
  }

  # extinction_probs = list()
  qouts = list()
 
  for(r in 1:length(resist_params$y)){

    pfull[[hetero_param]] = resist_params$y[r]

    # Not including zoospores so only adding one
    q_res = with(pother$ssmod,
                array(NA, dim=c(bins + tad_stages + juv_stages + 1, length(vals) + 1)))
    q_res[, 1] = 0

    slope = pfull[[hetero_param]]
    demo_matrix = build_demography_matrix(pother, slope, y, nrow(q_res) + 1, 0, 
                                            dropZ=TRUE)

    for(i in 1:length(vals)){

      time = vals[i]

      q_up = array(NA, dim=c(nrow(q_res), time)) # Hold results for the recursive for loop
      start = array(0, dim=c(nrow(q_up), 1))

      repro_correct_temp = !pother$sim$repro_first

      for(j in time:1){ # Reverse loop through time
          
          # Compute the transition matrix for the given time step
          # Should be j - 1 because the zeroth temp is the start_temp

          res = build_time_dependent_pgf(j - 1, pfull, pother, bins, bnd, y, h,
                          repro_correct_temp, demo_matrix, elas_vals)

          full_P = res$ft
          repro_correct_temp = res$repro

          if(j == time){
              up_ext = compute_extinction_probs(full_P, 1, start)
          } else{
              up_ext = compute_extinction_probs(full_P, 1, q_up[, j + 1])
          }

          q_up[, j] = up_ext[, 2]

      }

      q_res[, i + 1] = q_up[, 1]

    }

    #extinction_probs[[r]] = q_res[4, ]

    qouts[[r]] = apply(q_res[1:4, ]^init_full[r, ], 2, prod)


  }

  pop_extinct = apply(t(simplify2array(qouts)), 2, prod)
  
  # extinct_matrix = do.call(rbind, extinction_probs)

  # init_matrix = matrix(rep(init_vect, pother$sim$years + 1), nrow=length(init_vect), 
  #                                 ncol=pother$sim$years + 1)

  # pop_extinct = apply(extinct_matrix^init_matrix, 2, prod)

  years = with(pother$sim, (c(0, vals) * step_length) / days_per_year)

  return(list(probs=pop_extinct, time=years))

}


extinction_prob_annual = function(pfull, pother, standing_var, resist_bins,
                hetero_param="growth_size",
                elas_vals=list(elasphi=1,elasg0=1,elasg=1,elass=1,elasr=1,elasz=1)){
  # Calculate extinction probabilities using a multi-type branching process
  # over an entire year. Assumes density-independent recruitment and 
  # transmission in order to use a branching process approximation. 
  #
  # Parameters
  # ----------
  # pfull : list-like object with temperature-dependent IPM parameters
  # pother : list-like objects with IPM and simulation parameters specified 
  #            in model_parameters.yml
  # standing_var : Float, the amount of standing variation in resistance.
  #                If 0, resist_bins is set to 1.
  # resist_bins : int, the number of bins used to approximate the continuous
  #                resistance trait.  Approximating the continuous with the
  #                mid-point rule
  # hetero_param : str, the parameter in pfull that shows heterogeneity between
  #                hosts.  By default this is growth_size.
  # elas_vals : list of floats, elasticity perturbation of the vital rate 
  #          functions
  #         elasphi : infection/transmission function 
  #         elasg0 : initial infection function
  #         elasg : growth function
  #         elass : survival function
  #         elasr : recovery function
  #         elasz : zoospore decay 
  #
  # Results
  # -------
  # : the yearly extinction probabilities calculated from a branching process

  # Setting densities and proportions for branching process to work
  pother$ipm$num_inf = 2 # Set to 2 to match fully stochastic simulation
  pother$ipm$prop_inf = 1
  pother$ipm$zoo_inf = 10000

  # If there is not variability there is only one resistance bin
  if(standing_var == 0){
        resist_bins = 1
    } 

  # Unpack elasticity vals
  elasphi = elas_vals$elasphi; elasg0 = elas_vals$elasg0; 
  elasg = elas_vals$elasg; elass = elas_vals$elass; elasr = elas_vals$elasr
  elasz = elas_vals$elasz

  # Set the basic IPM matrix parameters
  matrix_params = with(pother$sim,
                          set_discretized_values(min_size, max_size, bins))
  bnd = matrix_params$bnd
  y = matrix_params$y
  h = matrix_params$h
  bins = pother$sim$bins

  ##########################################################################

  resist_params = get_resistance_params(pfull[[hetero_param]], 
                                        sqrt(standing_var), resist_bins)
  # init_vect = assign_initial_pop_values(pother$ssmod$init_adult_val, 
  #                         pfull[[hetero_param]],
  #                         sqrt(standing_var),
  #                         resist_params)

  # Get adult and tadpole initial values
  init_full = sapply(c(pother$ssmod$init_tad_vals, pother$ssmod$init_adult_val), 
                          assign_initial_pop_values, 
                          pfull[[hetero_param]], sqrt(standing_var),
                          resist_params)
  if(resist_bins == 1){
    init_full = t(as.matrix(init_full))
  }

  starting_temp = pother$sim$start_temp
  all_temps = get_temp(1:121, starting_temp)
  #extinction_probs = list()
  qouts = list()

  # Loop through different levels of resistance
  for(s in 1:length(resist_params$y)){

      # Set the growth slope
      pfull[[hetero_param]] = resist_params$y[s]

      params = set_temp_params(starting_temp, pfull, pother$ipm$linear_temp, 
                                  elasz=elasz)
      pother$ipm$const_temp = params$prob_inf
      params$prob_inf = get_infection_probability(pother, pfull,
                      pother$ipm$transmission_type, elas=elasphi)

      kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, 
                              elasg=elasg, elass=elass, elasr=elasr)
      full_P = build_transition_extinction(params, pother, bins, bnd, y, h, 
                                          kernels, elasg0=elasg0)

      # Demography matrix and fertility
      demo_matrix = build_demography_matrix(pother, pfull[[hetero_param]], y, 
                                              nrow(full_P) + 1, 0, dropZ=TRUE)

      # Accounts for multiple temperatures within a season
      repro_now = pother$sim$repro_first

      # Calculate the seasonal matrix.  This is using the method from Chapter
      # 13 of the Caswell and combining each time step into a single, yearly
      # transition matrix. Assuming density independence as necessary for the
      # branching process
      for(i in 1:length(all_temps)){

          params = set_temp_params(all_temps[i], pfull, pother$ipm$linear_temp,
                                    elasz=elasz)
          pother$ipm$const_temp = params$prob_inf # Setting this if constant temperature is used
          params$prob_inf = get_infection_probability(pother, pfull,
                                    pother$ipm$transmission_type, elas=elasphi)

          kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, 
                                      elasg=elasg, elass=elass, elasr=elasr)
          tfull_P = build_transition_extinction(params, pother, bins, bnd, y, 
                                    h, kernels, elasg0=elasg0)
          full_P = tfull_P %*% full_P

          # Decide if the temperature and the time is right for reproduction
          if(round(all_temps[i], 2) == pother$sim$temp_at_repro){
              if(repro_now){
                  full_P = demo_matrix %*% full_P
                  repro_now = !repro_now
              } else{
                  repro_now = !repro_now
              }
          }

      }

      time = pother$sim$years # Time over which to explore extinction probability

      # Compute the probability of population extinction using a branching process
      q_start = array(0, dim=c(nrow(full_P), 1))
      qout = compute_extinction_probs(full_P, time, q_start)

      # Extract the probability for a single uninfected adult
      #extinction_probs[[s]] = qout[4, ]

      # Compute extinction prob for initial structure and each extinction class
      qouts[[s]] = apply(qout[1:4, ]^init_full[s, ], 2, prod)


  }


  # Then compute the probability of extinct given the initial population
  # structure stored in init_full
  pop_extinct = apply(t(simplify2array(qouts)), 2, prod)


  # # Calculate probability of extinction with just adults
  # extinct_matrix = do.call(rbind, extinction_probs)

  # init_matrix = matrix(rep(init_vect, time + 1), nrow=length(init_vect), 
  #                               ncol=time + 1)

  # pop_extinct = apply(extinct_matrix^init_matrix, 2, prod)

  # Calculate 
  return(pop_extinct)
  #return(list(pop_ext=pop_extinct, full_probs=full_probs))

}

get_prob_z_inf = function(params, pother){
  # Calculates the probability of gaining infection from zoospore pool. 
  # 
  # You can do it by just gaining an infection from the zoospore pool or
  # gaining it from both the zoospore pool and interactions with 
  # infected individuals.
  # 
  # Parameters
  # ----------
  # params, pother : list of model parameters
  #
  # Results
  # -------
  # : prob infection from either just th zoospore pool or both zoospore pool
  #   and contact with an infected host.


  tpother = pother
  tpother$ipm$num_inf = 0
  tpother$ipm$prop_inf = 0
  prob_inf_z = get_infection_probability(tpother, pfull, pother$ipm$transmission_type)

  # Get prob infection from only DD or FD
  tpother = pother
  tpother$ipm$zoo_inf = 0
  prob_inf_d = get_infection_probability(tpother, pfull, pother$ipm$transmission_type)

  # First term: Prob inf just from zoospore pool, Second term: prob_inf from both
  return(prob_inf_z * (1 - prob_inf_d) + prob_inf_z*prob_inf_d)

}

nonstochastic_step = function(full_P, n_t, params, pother){
  # The non-stochastic step that accounts for the reduction of zoospores in the
  # zoospore pool.
  #
  # Parameters
  # ----------
  # full_P : matrix, the full transition matrix
  # n_t : vector, the vector of observed abundances
  # params : list, The temperature-dependent parameters
  # pother : list, additional parameters.
  #
  # Returns
  # -------
  # : population vector at time t + 1

  # Probability of gaining infection from zoospore pool. You can do it by
  # just gaining an infection from the zoospore pool or gaining it from
  # both the zoospore pool and interactions with infected individuals.
  phi_z = get_prob_z_inf(params, pother)
  S_t = n_t[pother$ssmod$tad_stages + 1] # Extract uninfected adults
  mu_clump = params$clump_mean # Mean number of zoospores gained per clump
  Z_t = n_t[length(n_t)] # Number of zoospores in pool
  reduce = min(c(phi_z * S_t * exp(mu_clump + params$clump_sd^2 / 2), 
                                            Z_t * params$zspore_surv))

  # Update matrix
  n_t_new = full_P %*% n_t

  # Reduce zoospore pool
  n_t_new[length(n_t_new)] = n_t_new[length(n_t_new)] - reduce

  return(n_t_new)

}

stochastic_step = function(full_P, n_t, params, pother){
  # Perform a stochastic step based on the given transition matrix
  # Splits the transition matrix into augmented T matrix and (F)ertility
  # Assumes Poisson reproduction for F and multinomial transitions for
  # augmented T.  This is specific for the frog-Bd IPM.
  # 
  # Parameters
  # ----------
  # full_P : matrix, The transition matrix for a given time step
  # n_t : The population vector at time t - 1
  # params : list, The temperature-dependent parameters
  # pother : list, additional parameters.
  #
  # Returns
  # -------
  # : Updated n_t incorporating demographic stochasticity

  # Step 1: Split into T_aug and F
  T = full_P
  T[1, -1] = 0 # Remove the Fertility matrix, but not tadpoles

  T[nrow(T), -ncol(full_P)] = 0 # Remove the zoospore fertility, but leave survival
  F = array(0, dim=dim(full_P))
  F[1, ] = full_P[1, ] # Create fertility matrix
  F[1, 1] = 0 # Drop tadpole transition in the fertility matrix
  F[nrow(T), ] = full_P[nrow(full_P), ]
  F[nrow(T), ncol(T)] = 0 # Drop zoospore transition in the fertility matrix

  T_aug = rbind(T, 1 - colSums(T)) # Create an augmented T
  T_aug[T_aug < 0] = 0 # This just catches some rounding errors that R makes
  T_aug = t(t(T_aug) / colSums(T_aug)) # Make sure all columns are summing to one

  tmat = array(NA, dim=dim(full_P))

  n_t = c(n_t)
  for(i in 1:length(n_t)){

    # Draw transition from a multinomial
    #t_pred = rmultinom(1, n_t[i], T_aug[, i])

    # Multinomial can overflow so just do this deterministically when
    # zoospore load gets too big
    t_pred = tryCatch({

            tpred = rmultinom(1, n_t[i], T_aug[, i])

          }, error = function(err){

            # Just make the transition deterministic if overflow 
            dead = round(T_aug[nrow(T_aug), i] * n_t[i])
            alive = n_t[i] - dead
            tpred = matrix(c(rep(0, (nrow(T_aug) - 2)), alive, dead), 
                            nrow=nrow(T_aug), ncol=1)
            return(tpred)

          })

    # Reproduction 
    f_pred = rpois(1, n_t[i]*F[1, i])
    z_pred = rpois(1, n_t[i]*F[nrow(T), i])

    z_left = t_pred[nrow(T_aug) - 1, ]

    # Add tadpoles and zoospores onto what is already there
    t_pred[1, ]  = t_pred[1, ] + f_pred
    t_pred[nrow(T_aug) - 1, ] = t_pred[nrow(T_aug) - 1, ] + z_pred
    tmat[, i] = t_pred[1:(nrow(T_aug) - 1), ]
  }

  n_t_new = as.matrix(rowSums(tmat), nrow=nrow(full_P), ncol=1)

  # Reduce zoospore pool based on gain from from susceptibles
  phi_z = get_prob_z_inf(params, pother)
  S_t = n_t[pother$ssmod$tad_stages + 1] # Extract uninfected adults
  S_t_inf = rbinom(1, S_t, phi_z)

  gained_z = round(sum(exp(rnorm(S_t_inf, mean=params$clump_mean, 
                    sd=params$clump_sd))))

  reduce = min(c(gained_z, z_left))

  # Reduce the zoospore pool based on zoospores that were acquired from
  # from zpool
  n_t_new[nrow(n_t_new), 1] = n_t_new[nrow(n_t_new), 1] - reduce

  return(n_t_new)
}

build_time_dependent_pgf = function(time, pfull, pother, bins, bnd, y, h,
                repro_correct_temp, demo_matrix, elas_vals){
  # Builds the time/temperature dependent kernel given some time increment
  # The resulting kernel can be converted to a PGF for extinction time 
  # calculations
  # 
  # Parameters
  # ----------
  # time: A time step integer from 0 to something large
  # pfull: A list of parameters from the empirically fitted IPM
  # pother: A list of parameters from the YAML file
  # bins: Number of bins for the IPM
  # bnd: Boundaries for the discretized IPM
  # y: Midpoints for the discretized IPM
  # repro_correct_temp: bool, if TRUE reproduction occurs the first time a
  #                    temp_at_repro occurs. Otherwise, the second.
  # demo_matrix: The demography transition matrix.
  # elas_vals: list of elasticity values
  #
  # Returns
  # -------
  # : list with the full transition matrix and the updated reproductive 
  #           boolean

  # Unpack elasticity vals
  elasphi = elas_vals$elasphi; elasg0 = elas_vals$elasg0; 
  elasg = elas_vals$elasg; elass = elas_vals$elass; elasr = elas_vals$elasr
  elasz = elas_vals$elasz
  
  temp = get_temp(time, pother$sim$start_temp)
  params = set_temp_params(temp, pfull, pother$ipm$linear_temp, elasz=elasz)

  # Setting approximate values depending on what transmission function is used
  pother$ipm$num_inf = 2; pother$ipm$prop_inf = 1; pother$ipm$zoo_inf = 10000
  pother$ipm$const_temp = params$prob_inf

  params$prob_inf = get_infection_probability(pother, pfull, 
                              pother$ipm$transmission_type, elas=elasphi)

  kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, 
                                        elasg=elasg, elass=elass, elasr=elasr)
  full_transition = build_transition_extinction(params, pother, bins, bnd, y, 
                                    h, kernels, elasg0=elasg0)

  if(round(temp, 2) == pother$sim$temp_at_repro){
      # This if/else prevents double reproduction occurring
      if(repro_correct_temp){ # Correct instance of temperature
          full_transition = demo_matrix %*% full_transition
          repro_correct_temp = !repro_correct_temp
      } else{ 
          repro_correct_temp = !repro_correct_temp
      }
  } 

  return(list(ft=full_transition, repro=repro_correct_temp))
    
}


build_transition_extinction = function(params, pother, bins, bnd, y, h, kernels,
                                            elasg0=1){
  # Builds the density-independent transition matrix used in the branching
  # process approximation
  #
  # Parameters
  # ----------
  # params: parameters for the IPM fit from data
  # pother: Parameters stored in yaml file
  # bins, bnd, y, h : discretization parameters
  # kernels : density independent IPM kernels

  # Returns
  # -------
  # list of updated pop vector and infection probability


  # Calculate the time variant portion of the kernel. This will be updated
  # each time step based on the distribution of zoospores.  Initiall assuming
  # density dependent transmission. The total number of zoospores across all
  # frogs.
  infection_prob = params$prob_inf

  zero_trans = params$class_zero_surv * (1 - infection_prob)

  # The zero to nonzero transitions: Transitioning from zero to infected
  prob_clump = h * g0_x(y, params, elas=elasg0) # dnorm(y, mean=params$clump_mean, sd=params$clump_sd)
  ztnz = params$class_zero_surv * infection_prob * prob_clump
  col1 = c(zero_trans, ztnz)

  # Compute the infected to uninfected row. This is just the probability of
  # surviving with load $x$ ($s(x)$) and then losing an infection $r(x)$.  We
  # just need to multiply $s(x) r(x)$

  # The loss probabilities: Transitioning from non-zero to zero
  nztz = kernels$S * kernels$R
  row1 = c(zero_trans, nztz)

  # Add the infection probabilities onto the P matrix
  full_P = get_full_P(kernels$P, row1, col1)

  # Add the three tadpole classes onto the matrix
  # Assuming tadpoles are immediately infected and have constant load
  full_tad = diag(nrow(full_P) + pother$ssmod$tad_stages)
  # Set tadpole survival probability for each 3 day time step.
  diag(full_tad)[1:pother$ssmod$tad_stages] = pother$ipm$tad_survival
  tad_rows = full_tad[1:pother$ssmod$tad_stages, ]
  tad_cols = full_tad[, 1:pother$ssmod$tad_stages]
  #tad_cols[nrow(tad_cols), ] = pother$ipm$mean_tad_load

  # Add on tadpole columns
  full_P = get_full_P(full_P, tad_rows, tad_cols,
                    num_add=pother$ssmod$tad_stages, add_top=TRUE)

  return(full_P)

}

###########################################################################

build_full_P = function(min_size, max_size, bins, params) {
  # Function that builds the full transition matrix for the host parasite
  # IPM
  #
  # Parameters
  # ---------- 
  # min_size: float, lower bound of the matrix
  # max_size: float, upper bound of the matrix
  # bins: int, number of bins in the matrix
  # params: list, list of parameters to use to build P

  # Returns
  # -------
  # full_P : full transition matrix

  # Get the relevant matrix parameters: midpoints, cell edges, and cell width
  matrix_params = set_discretized_values(min_size, max_size, bins)
  bnd = matrix_params$bnd
  y = matrix_params$y
  h = matrix_params$h

  # Calculate and plot eviction.  Eviction is most important when e is large and. Using the script provided by Williams et al. for checked for eviction

  # Get the kernel of the IPM without the zero class
  kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, plot_it=F)
  P = kernels$P

  # The zero -> zero transition
  zero_trans = params$class_zero_surv * (1 - params$prob_inf)

  # The zero to nonzero transitions: Transitioning to zero to infected
  prob_clump = h * dnorm(y, mean=params$clump_mean, sd=params$clump_sd)
  ztnz = params$class_zero_surv * params$prob_inf * prob_clump
  col1 = c(zero_trans, ztnz)

  # The loss probabilities: Transitioning from non-zero to zero
  nztz = kernels$S * kernels$R
  row1 = c(zero_trans, nztz)

  # Now we just need to stick these extra columns onto the **P** matrix. The
  # most likely thing that is going to happen is you are going to stay in the
  # uninfected class.

  full_P = get_full_P(P, row1, col1, min_size, y, plot_it=F)

  return(full_P)

}

###########################################################################

monte_carlo_params = function(params){
  # For the various vital functions in the IPM, uses the vcov matrix for each
  # set of params to draw new params (assuming multivariate normality under
  # asymptotic likelihood theory) and return a params matrix with the updated
  # params

  # Parameters
  # ----------
  # params : list, properly formatted list of parameters from ipm_vital*.R

  # Returns
  # -------
  # rand_params : list, list where growth, loss, surv, clump, and prob_inf
  # parameters are random draws

  rand_params = params
  vital_names = c('growth', 'loss', 'surv', 'clump', 'prob_inf')

  for(vital_name in vital_names){

    # Draw params for each vital function
    mean_vect = paste(vital_name, "coef", sep="_")
    vcov_mat = paste(vital_name, "vcov", sep="_")
    #print(paste(mean_vect, vcov_mat))
    vital_draw = mvrnorm(n=1, mu=params[[mean_vect]], Sigma=params[[vcov_mat]])
    vital_draw = as.list(vital_draw)

    # Store params if they exist
    rand_params[[paste(vital_name, "int", sep="_")]] = vital_draw[["(Intercept)"]]
    rand_params[[paste(vital_name, "temp", sep="_")]] = vital_draw[["temp"]]
    rand_params[[paste(vital_name, "size", sep="_")]] = vital_draw[["size"]]
    rand_params[[paste(vital_name, "size_sq", sep="_")]] = vital_draw[["I(size^2)"]]
    rand_params[[paste(vital_name, "size_temp_int", sep="_")]] = vital_draw[["size:temp"]]

  }

  return(rand_params)

}


bayesian_params = function(params, vital_names=c('growth', 'loss', 'surv', 
                                                 'clump', 'zfd', 'zdd', 'zconst')){
  # Draw a sample from the posterior distribution for each parameter.
  #
  # This function uses the Bayesian estimates of each vital rate functions
  # to appropriately account for the correlation between parameters within a
  # vital rate function. This assumes that each vital rate function given
  # in `vital_names` has a slot named `*_samples`.
  #
  # Parameters
  # ----------
  # params : list containing posterior samples and parameter for
  #          vital rate functions 
  # vital_names : array-like, contains prefixes for the names of the vital
  #               rate functions from which a posterior sample should be
  #               drawn. Defaults to all vital rate functions.
  # 
  # Returns
  # -------
  # : updated parameter list with parameters drawn from the posterior of each
  #   parameter

  rand_params = params

  for(vital_name in vital_names){ # Loop through

    samples = params[[paste(vital_name, "samples", sep="_")]]
    rsamp = samples[sample(1:nrow(samples), 1), ]

    if(length(colnames(samples)) != 1){

      for(pname in colnames(samples)){
        rand_params[[pname]] = as.numeric(rsamp[pname])
      }

    } else{
      rand_params[[colnames(samples)]] = rsamp
    }
  }

  return(rand_params)

}

perturb_params = function(params, param_names, sigma){
  # Perturb parameters by a log-normal amount with median 1 and log variance
  # sigma^2. The sign of the parameter will not be changed. This is the
  # the method proposed by Sobie (2009)
  #
  # Parameters
  # ----------
  # params : list-like object containing parameters
  # param_names : vector of parameter names
  # sigma : sdlog of lognormal (rlnorm)
  #
  # Returns
  # -------
  # : updated list of parameters with the random perturbation

  updated_params = params

  for(p in param_names){

    # This won't change the sign of the parameter
    updated_params[[p]] = params[[p]]*rlnorm(1, sdlog=sigma)

  }

  return(updated_params)
}

################################

compute_quartiles = function(est_matrix){
  # Computes the median, lower quartile, and upper quartile of est_matrix
  # column wise

  median = apply(est_matrix, 2, median)
  upper = apply(est_matrix, 2, uq)
  lower = apply(est_matrix, 2, lq)

  return(list(lower=lower, median=median, upper=upper))

}

compute_extinction_probs = function(full_P, time, q){
  # Uses the theory presented in Schreiber and Ross 2016 and Caswell 2001
  # to build generating functions. Assumes poisson reproduction
  #
  # Parameters
  # ----------
  # full_P : Transition matrix
  # time : Number of time steps at which to calculate the extinction
  #        probability
  # q : starting probabilities of extinction

  # Using code from Caswell page 494 to compute extinction probabilities
  # Involves calculating the PGF and then predicting the extinction probability

  q = as.matrix(q)
  T = full_P

  # Remove the Fertility matrix, but don't remove the first item which gives
  # the probability of item 1 transitioning to itself (survival probability)
  T[1, -1] = 0 
  F = array(0, dim=dim(full_P))
  F[1, ] = full_P[1, ] # Create fertility matrix
  F[1, 1] = 0 # Don't include the tadpole transition

  T_aug = rbind(T, 1 - colSums(T)) # Create an augmented T
  T_aug[T_aug < 0] = 0 # This just catches some rounding errors that R makes
  T_aug = t(t(T_aug) / colSums(T_aug)) # Make sure all columns are summing to one

  qprime = array(0, dim=c(nrow(F), 1))
  k = nrow(F)
  qout = array(NA, dim=c(nrow(F), time + 1))
  qout[, 1] = q


  for(l in 2:(time + 1)){ # Iterate the generating function through time

      for(j in 1:k){

          gt = t(as.matrix(T_aug[, j])) %*% as.matrix(c(q, 1)) # Multitype PGF 
          gb = exp(F[1, j] * (q[1, ] - 1)) # Poisson birth PGF, 1 - F[1, j] + F[1, j] * q[1, ] # Bernoulli birth
          qprime[j, ] = gt * gb

      }

      q = qprime 
      qout[, l] = q

  }

  return(qout)

}


simulation_extinction_times = function(sim_results){
  # For the results of a multiseason simulation, get whether or not
  # a population is extinction at any given time point.
  #
  # Parameters
  # ----------
  # sim_results : list, output from multiseason_simulation
  #
  # Returns
  # -------
  # : array, contains booleans specifying whether a population went extinct 
  # or not at a given time point for a given simulation.

  sims = length(sim_results)
  rbins = length(sim_results[[1]]) # resistance bins
  ibins = nrow(sim_results[[1]][[1]])
  steps = ncol(sim_results[[1]][[1]])
  extinct_results = array(NA, dim=c(sims, steps))


  for(s in 1:sims){

    sim_pops = array(NA, dim = c(rbins, steps))

    for(r in 1:rbins){
      pop = colSums(sim_results[[s]][[r]][1:(ibins - 1), ])
      sim_pops[r, ] = pop == 0 # Check extinct
    }

    extinct_results[s, ] = apply(sim_pops, 2, all)
 
  }

  return(extinct_results)

}

simulation_trajectories = function(sim_results){
  # Extracts the simulation trajectories from sim_results
  # Returns both the total population trajectories, collapsing on all resist
  # bins for each simulation (pop_trajs) and the simulations for each
  # resistance group (resist_groups)
  #
  # ASSUMING THREE TADPOLES
  #
  # Returns
  # -------
  # : list, pop_traj contains full populations trajectories by simulation
  #

  sims = length(sim_results)
  rbins = length(sim_results[[1]]) # resistance bins
  ibins = nrow(sim_results[[1]][[1]])
  steps = ncol(sim_results[[1]][[1]])
  pop_trajs = array(NA, dim=c(sims, steps))
  resist_groups = list()
  pop_trajs_tads = array(NA, dim=c(sims, steps))
  resist_groups_tads = list()
  zoo_trajs = array(NA, dim=c(sims, steps))

  for(s in 1:sims){

    resist_pops = array(NA, dim=c(rbins, steps))
    resist_pops_tads = array(NA, dim=c(rbins, steps))
    resist_pops_zoos = array(NA, dim=c(rbins, steps))

    for(r in 1:rbins){

      # Adult trajectories
      pop = colSums(sim_results[[s]][[r]][4:(ibins - 1), ])
      resist_pops[r, ] = pop

      # Tadpole trajectories
      pop_tads = colSums(sim_results[[s]][[r]][1:3, ])
      resist_pops_tads[r, ] = pop_tads

      # Zoospore trajectories
      pop_zoos = sim_results[[s]][[r]][ibins, ]
      resist_pops_zoos[r, ] = pop_zoos

    }

    pop_trajs[s, ] = apply(resist_pops, 2, sum)
    resist_groups[[s]] = resist_pops

    pop_trajs_tads[s, ] = apply(resist_pops_tads, 2, sum)
    resist_groups_tads[[s]] = resist_pops_tads

    zoo_trajs[s, ] = apply(resist_pops_zoos, 2, sum)

  }

  return(list(zoo_trajs=zoo_trajs, pop_trajs_tads=pop_trajs_tads, 
            pop_trajs=pop_trajs, resist_groups=resist_groups))

}

simulation_load_dynamics = function(sim_results, y){
  # Assuming three tadpole classes
  # 
  # Parameters
  # ----------
  # sim_results : results from multiseason simulation
  # y : vector with mesh points for IPM approximation
  # 
  # Returns
  # -------
  # : matrix load of infected through time. ncol=time steps, nrow=simulations

  sims = length(sim_results)
  rbins = length(sim_results[[1]]) # resistance bins
  ibins = nrow(sim_results[[1]][[1]])
  steps = ncol(sim_results[[1]][[1]])

  mean_results = var_results = prev_results = array(NA, dim=c(sims, steps))

  tot_pop = simulation_trajectories(sim_results)$pop_trajs

  for(s in 1:sims){ # Loop through simulation

    resist_means = resist_vars = r_weights = array(NA, dim=c(rbins, steps))
    uninfecteds = array(NA, dim=c(rbins, steps))

    for(r in 1:rbins){ # Loop through resistance classes

      # Get proportion of adults in each infection class
      adults = sim_results[[s]][[r]][4:(ibins - 1), ]
      norm_adults = t(t(adults) / colSums(adults))


      r_weights[r, ] = colSums(adults) / tot_pop[s, ] # Weights for mixture

      # Get normalized infected adults
      norm_infected = t(t(norm_adults[2:nrow(norm_adults), ]) / 
                                              (1 - norm_adults[1, ]))

      # If there are no infected the mean should be NA
      # norm_infected[, colSums(norm_infected) == 0] = NA

      mean_by_time = colSums(norm_infected * y)
      var_by_time = colSums(norm_infected * y^2) - mean_by_time^2
      resist_means[r, ] = mean_by_time
      resist_vars[r, ] = var_by_time
      uninfecteds[r, ] = adults[1, ]

    }

    # Compute mean and variance from mixture of resistant types. Sometimes
    # all types might be nan so we need to account for this.  We do that
    # with the ugly code below

    not_all_na = apply(resist_means, 2, function(x) !all(is.na(x)))

    mean_results[s, not_all_na] = colSums(resist_means[, not_all_na, drop=F] * 
                        r_weights[, not_all_na, drop=F], na.rm=T)
    rm_no_na = as.matrix(resist_means[, not_all_na, drop=F], )

    var_results[s, not_all_na] = colSums((resist_means[, not_all_na, drop=F] + resist_vars[, not_all_na, drop=F])
                      * r_weights[, not_all_na, drop=F], na.rm=T) - mean_results[s, not_all_na, drop=F]
    prev_results[s, ] = (1 - colSums(uninfecteds) / tot_pop[s, ])

  }

  return(list(means=mean_results, vars=var_results, prevs=prev_results))

}

trapz = function (x, y){
    # Function trapz from `pracma` package. Computes the area under the curve
    # defined by x and y.

    if (missing(y)) {
        if (length(x) == 0) 
            return(0)
        y <- x
        x <- seq(along = x)
    }
    if (length(x) == 0 && length(y) == 0) 
        return(0)
    if (!(is.numeric(x) || is.complex(x)) || !(is.numeric(y) || 
        is.complex(y))) 
        stop("Arguments 'x' and 'y' must be real or complex vectors.")
    m <- length(x)
    if (length(y) != m) 
        stop("Arguments 'x', 'y' must be vectors of the same length.")
    if (m <= 1) 
        return(0)
    xp <- c(x, x[m:1])
    yp <- c(numeric(m), y[m:1])
    n <- 2 * m
    p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
    p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
    return(0.5 * (p1 - p2))

}

############### Plotting functions ############################

plot_simulation = function(sim_results, pfull, pother){
  # Plotting simulation results

  matrix_vals = with(pother$sim, 
                  set_discretized_values(min_size, max_size, bins))

  steps = with(pother$sim, (days_per_year * years) / step_length)

  time = with(pother$sim, (c(0, 1:steps) * step_length) / days_per_year)

  pop_trajs = simulation_trajectories(sim_results)

  # Format adults for plotting
  dt_pjs = data.table(t(pop_trajs$pop_trajs))
  dt_pjs$time = time
  dt_melt = melt(dt_pjs, id.vars="time")
  colnames(dt_melt) = c("time","sim", "pop_size")

  adult_plot = ggplot(data=dt_melt, aes(x=time, y=log(pop_size + 1))) + geom_line(aes(color=sim), alpha=0.4) + 
                    theme_bw() + guides(color=FALSE) +
                    theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
                    xlab("Time (years)") + ylab("ln(Adults + 1)")

  # Format tadpoles for plotting
  dt_pjs_tad = data.table(t(pop_trajs$pop_trajs_tads))
  dt_pjs_tad$time = time
  dt_melt_tad = melt(dt_pjs_tad, id.vars="time")
  colnames(dt_melt_tad) = c("time","sim", "tad_pop_size")

  tadpole_plot = ggplot(data=dt_melt_tad, aes(x=time, y=log(tad_pop_size + 1))) + 
                    geom_line(aes(color=sim), alpha=0.4) + 
                    theme_bw() + guides(color=FALSE) +
                    theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
                    xlab("Time (years)") + ylab("ln(Tadpoles + 1)")

  # Format mean loads for plotting
  load_res = simulation_load_dynamics(sim_results, matrix_vals$y)

  if(all(is.na(load_res$means))){
    load_res$means[is.na(load_res$means)] = 0

  }

  dt_load = data.table(t(load_res$means))

  dt_load$time = time
  dt_melt_load= melt(dt_load, id.vars="time")
  colnames(dt_melt_load) = c("time","sim", "bd_load")

  load_plot = ggplot(data=dt_melt_load, aes(x=time, y=bd_load)) + 
                    geom_line(aes(color=sim), alpha=0.4) + 
                    theme_bw() + guides(color=FALSE) +
                    theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
                    xlab("Time (years)") + ylab("Mean ln Bd load")

  # Format zoospores for plotting
  dt_pjs_zoo = data.table(t(pop_trajs$zoo_trajs))
  dt_pjs_zoo$time = time
  dt_melt_zoo = melt(dt_pjs_zoo, id.vars="time")
  colnames(dt_melt_zoo) = c("time","sim", "zoo_pop_size")

  zoo_plot = ggplot(data=dt_melt_zoo, aes(x=time, y=log(zoo_pop_size + 1))) + 
                    geom_line(aes(color=sim), alpha=0.4) + 
                    theme_bw() + guides(color=FALSE) + 
                    theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) + 
                    xlab("Time (years)") + ylab("Ln(zoospores + 1)")


  grid_plot = grid.arrange(adult_plot, tadpole_plot, load_plot, zoo_plot, nrow=2, ncol=2)
  return(grid_plot)

}
