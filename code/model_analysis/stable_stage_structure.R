## Description
## -----------
##
## Calculate the stable stage structure for a density-independent, disease
## free frog population. We use these results to set our initial conditions
## for the R. muscosa-Bd model
##
## Author: Mark Wilber

library(yaml)
pother = yaml.load_file("model_parameters.yml")
pfull = readRDS("../../results/bayesian_parameter_estimates.rds")


# Build the density-independent projection matrix

pmat = matrix(NA, nrow=4, ncol=4)
steps = pother$sim$days_per_year / pother$sim$step_length
surv_probs = pother$ipm$tad_survival^steps # yearly tadpole survival prob
nometa_probs = pother$ssmod$tad_not_adult # Probability of tadpole metamorphosing
meta_surv_probs = pother$ssmod$tad_meta_surv # Probability of surviving meta
adult_surv = pfull$class_zero_surv^steps


pmat[1, ] = c(0, 0, 0, (pother$ssmod$repro) * pother$ssmod$uninfected_adult_repro_prob) # T1 row
pmat[2, ] = c(surv_probs[1] * nometa_probs[1], 0, 0, 0) # T2 row
pmat[3, ] = c(0, surv_probs[2] * nometa_probs[2], 0, 0) # T3 row
pmat[4, ] = c(surv_probs[1] * (1 - nometa_probs[1]) * meta_surv_probs[1], 
              surv_probs[2] * (1 - nometa_probs[2]) * meta_surv_probs[2],
              surv_probs[3] * (1 - nometa_probs[3]) * meta_surv_probs[3],
              adult_surv) # Adults

# Get stable stage distribution
w_eigen = Re(eigen(pmat)$vectors[, 1])
stab_dist = w_eigen / sum(w_eigen)

# Compute number of individuals in each class given initial adults
rel_init_pop = stab_dist / stab_dist[4]
absolute_init_pop = round(rel_init_pop * pother$ssmod$init_adult_val, 0)