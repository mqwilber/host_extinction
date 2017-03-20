data{
    int<lower=1> N; // Number of data points
    real sizeNext[N]; // Bd load at time t + 1
    real sizeNow[N]; // Bd load at time t
    int temp[N]; // Temperature
}
parameters{
    real b0; // Intercept effect 
    real b1; // Effect of Bd load at time t
    real b2; // Effect of temperature
    real<lower=0> sigma; // Temperature
    real delta; // Parameter for unequal variance
} 
model{
    vector[N] mu;
    vector[N] power_sd; 

    // Priors
    sigma ~ cauchy(0 ,1);
    b2 ~ normal(0, 5);
    b1 ~ normal(0, 5 );
    b0 ~ normal(0, 5);
    delta ~ normal(0, 5);

    // Define the model with unequal variance
    for ( i in 1:N ) {
        power_sd[i] <- sqrt(sigma * exp(2*delta*sizeNow[i]));
        mu[i] <- b0 + b1*sizeNow[i] + b2*temp[i];
    }
    sizeNext ~ normal( mu , power_sd );
}