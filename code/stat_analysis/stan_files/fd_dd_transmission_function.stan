data{
    int<lower=0> N; // Number of data points
    int infected[N]; // P/A infected

    real dens_or_prop[N]; // Density or proportion of frogs in tank
    real time[N]; // Time step

} parameters{

    // P/A predictors
    real<lower=0> b0; // Effect of time/zoospore pool
    real<lower=0> b1; // Effect of density or proportion infected

} model{
    
    // Priors
    b0 ~ cauchy(0, 1);
    b1 ~ cauchy(0, 1);

    for(i in 1:N){
        infected[i] ~ bernoulli(1 - exp(-1*(b0*time[i] + b1*dens_or_prop[i])));
    }
} generated quantities {
    
    vector[N] log_lik;

    for(i in 1:N){
        log_lik[i] <- bernoulli_log(infected[i], 1 - exp(-1*(b0*time[i] + b1*dens_or_prop[i])));
    }
    
}
