data{
    int<lower=0> N; // Number of data points
    int surv[N]; // Survived or died
    real bd_load[N]; // Log zoospore load

} parameters{

    // Survival predictors
    real b0; // Intercept
    real b1; // Effect of size of probability of death

} model{
    
    // Comment out for uninformative priors
    b0 ~ normal(0, 5);
    b1 ~ normal(0, 5);

    for(i in 1:N){
        surv[i] ~ bernoulli(inv_logit(b0 + b1*bd_load[i]));
    }
} 
//generated quantities {
//
//    vector[N] log_lik;
//
//    for(i in 1:N){
//        log_lik[i] <- bernoulli_log(surv[i], inv_logit(b0 + b1*bd_load[i]));
//    }
// }