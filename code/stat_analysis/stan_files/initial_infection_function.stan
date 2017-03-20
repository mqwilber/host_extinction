data{
    int<lower=1> N; // Length of data
    real sizeNext[N]; // Bd load at time t + 1
    real temp[N]; // Temperature
}
parameters{
    real b0; // Intercept
    real b1; // Effect of temperature
    real<lower=0> sigma; // Variance
    real delta; // Parameter for unequal variance
} 
model{
    vector[N] mu;
    vector[N] power_sd; 

    // Priors
    sigma ~ cauchy( 0 , 1 );
    b0 ~ normal( 0 , 5 );
    b1 ~ normal( 0 , 5 );
    delta ~ normal(0, 5);

    // Initial infection model with unequal variance
    for ( i in 1:N ) {
        power_sd[i] <- sqrt(sigma * exp(2*delta*temp[i]));
        mu[i] <- b0 + b1*temp[i];
    }
    sizeNext ~ normal( mu , power_sd );
}