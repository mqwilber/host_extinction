data{
    int<lower=0> N; // Number of data points
    int infected[N]; // P/A infected

    real time[N]; // Time

} parameters{

    // P/A predictors
    real<lower=0> b1; // Intercept, constant transmission

} model{
    
    // Priors
    b1 ~ cauchy(0, 1);

    for(i in 1:N){
        infected[i] ~ bernoulli(1 - exp(-1*(b1*time[i])));
    }
} generated quantities {

    vector[N] log_lik;

    for(i in 1:N){
        log_lik[i] <- bernoulli_log(infected[i], 1 - exp(-1*(b1*time[i])));
    }
}