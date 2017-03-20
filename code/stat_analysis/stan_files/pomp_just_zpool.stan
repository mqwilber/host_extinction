data{
    int<lower=0> N; // Number of time points
    int<lower=0> T; // Number of tanks
    int<lower=0> M; // Max number of individual frogs
    int<lower=0> num_frogs[T]; // Number of frogs in each tank

    real frog_zes[T, N]; // Holds frog zes by tank by time
    real density[T, N]; // Holds densities by tank by time
    real tad_zes[T, N]; // Holds tad zes by tank by time
    int prob_inf[T, M, N]; // Holds the 0 or 1 for infection probabilities
    real time[T, N]; // Holds time between swabs

    int init_tads; // Initial zoospore load on tadpoles
    int D; // Total number of data points that are contributing to the likelihood
    real<lower=0> sigma; // Fix variance of zoospore pool as data. Comment this out and uncomment sigma on line 23 and line 33 for sigma as a parameter.

} parameters{
    
    real<lower=0> b_zpool; // Effect of the zoospore pool on transmission
    real<lower=0> zpool[T, N]; // Holds the true, latent zpool values
    real<lower=0> zpool_init; // Initial zoospore values
    real<lower=0> death_rate; // The survival prob of zoospores
    real<lower=0> mu_z; // Contribution of adults and tadpoles to zpool
    // real<lower=0> sigma; // Variance of unobserved, latent zoospore pool

} model{
    
    real zpool_rate[T, N];

    // Priors
    b_zpool ~ cauchy(0, 1);
    zpool_init ~ lognormal(log(init_tads) - 0.5, 1);
    death_rate ~ normal(0.3, 0.03); // Prior taken from Woodhams et al. 2008
    mu_z ~ exponential(1);
    // sigma ~ cauchy(0, 3);


    for(t in 1:T){ // Loop over tanks

        for(i in 1:N){ // Loop over time points

            // By definition, there are no transitions at the first time point
            // All frogs are uninfected.
            if(i == 1) 
                zpool_rate[t, i] <- zpool_init * exp(-1*death_rate*time[t, i]);
            else{
                for(j in 1:M){ // Loop over individuals

                    // Because Stan doesnt handle ragged arrays need to ignore
                    // the jth row when it is greater than num_frogs
                    
                    if(j <= num_frogs[t]){ 

                        // The observation before needs to be 0 and the current
                        // observation needs to be not equal to -1
                        if(prob_inf[t, j, i - 1] == 0 && prob_inf[t, j, i] != -1){
                            prob_inf[t, j, i] ~ bernoulli(1 - exp(-1*time[t, i]*(b_zpool*log(zpool[t, i - 1] + 1)))); 
                        }
                    }
                }
                zpool_rate[t, i] <- zpool[t, i - 1]*exp(-1*death_rate*time[t, i]) + time[t, i]*(mu_z*tad_zes[t, i - 1] + mu_z*frog_zes[t, i - 1]);

            }

            zpool[t, i] ~ lognormal(log(zpool_rate[t, i]) - sigma^2 / 2, sigma);
        }
    }

} generated quantities {
    
    // Extracting the log-likelihoods for model comparison purposes
    vector[D] log_lik;
    int counter;
    real dev;
    counter <- 0;
    dev <- 0;

    for(t in 1:T){ // Loop over tanks

        for(i in 1:N){ // Loop over time points

            if(i != 1){
                for(j in 1:M){ // Loop over individuals

                    // Because Stan doesnt handle ragged arrays need to ignore
                    // the jth row when it is greater than num_frogs
                    
                    if(j <= num_frogs[t]){ 

                        // The observation before needs to be 0 and the current
                        // observation needs to be not equal to -1
                        if(prob_inf[t, j, i - 1] == 0 && prob_inf[t, j, i] != -1){
                            counter <- counter + 1;
                            log_lik[counter] <- bernoulli_log(prob_inf[t, j, i], 1 - exp(-1*time[t, i]*(b_zpool*log(zpool[t, i - 1] + 1)))); 
                            dev <- dev - (2 * bernoulli_log(prob_inf[t, j, i], 1 - exp(-1*time[t, i]*(b_zpool*log(zpool[t, i - 1] + 1)))));
                        }
                    }
                }

            }

        }
    }
}