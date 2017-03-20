data{
    int<lower=0> N; // Number of data points
    int loss[N]; // Lose infection or not
    real sizeNow[N]; // Bd load at time t
    real temp[N]; // Temperature

} parameters{

    // Survival predictors
    real b0; // Intercept
    real b1; // Effect of size of probability of death
    real b2; // Effect of temperature of probability of loss

} model{
    
    // Comment out for uninformative priors
    // b0 ~ normal(0, 5);
    // b1 ~ normal(0, 5);
    // b2 ~ normal(0, 5);

    for(i in 1:N){
        loss[i] ~ bernoulli(inv_logit(b0 + b1*sizeNow[i] + b2*temp[i]));
    }
} 