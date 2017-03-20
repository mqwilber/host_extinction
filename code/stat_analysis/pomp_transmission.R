# Description
# -----------
#
# This script runs different models for the latent zoospore pool.
# The script formats the data that needs to be passed into Stan and then runs
# the Stan models.  Formatting the data is a bit of a pain because Stan does
# not handle ragged arrays.  That is why there are so many extra matrices!
#
# Author: Mark Wilber

library(data.table)
library(rstan)
library(rethinking)
library(ggplot2)
library(GGally)
library(gridExtra)

# source("/Users/mqwilber/Dropbox/Documents/Stats/Stan/stan_fxns.R")

# Load in the archival data
arch_data = read.csv("../../data/archival/DensExpt.csv")
arch_data[arch_data$swab_id == "MDF16_13_Jul31", "ZE"] = 0 # This is a missing value that we are setting to 0
arch_data = arch_data[arch_data$swab_id != "MDF8_1_AUG7", ] # This was a repeat swab on the same individual so we are excluding it

# Unexplained crash in zoospore load after day 35 so we exclude. Explained
# in manuscript and Appendix.
arch_data = arch_data[arch_data$Day < 35 & !is.na(arch_data$Day), ] 

# Subset on adults
adults = arch_data[arch_data$Frog_or_Tadpole == F, ]
tads = arch_data[arch_data$Frog_or_Tadpole == T, ]

# Going to also want to calculate the number of infecteds at each time point
# Also the total number of zoospores at each time point 
adults = as.data.table(adults)
tads = as.data.table(tads)

# Keys for survey events...helping to group
surv_list = list(c(0, -3), 3:5, 9:11, 16:20, 23:25, 30:32)

surv_events = list()
count = 1
for(vals in surv_list){
    for(v in vals){
        surv_events[as.character(v)] = count
    }
    count = count + 1
}

labels = array(NA, dim=nrow(adults))
for(i in 1:nrow(adults)){
    labels[i] = surv_events[[as.character(adults$Day[i])]]
}

labels_tad = array(NA, dim=nrow(tads))
for(i in 1:nrow(tads)){
    labels_tad[i] = surv_events[[as.character(tads$Day[i])]]
}

adults$survey_event = labels
tads$survey_event = labels_tad

# Total adult zoospores in tank at time t
adult_ze = adults[, list(tot_ze = sum(ZE, na.rm=T), num_inf=sum(ZE > 0, na.rm=T)), 
                                    by=c("survey_event", "Tank_num")]
tad_ze = tads[, list(tot_ze_tad = sum(ZE, na.rm=T)), 
                                    by=c("survey_event", "Tank_num")]

# Merge tables...oh that is pleasant. Left Outer Join.
setkey(adult_ze, survey_event, Tank_num)
setkey(tad_ze, survey_event, Tank_num)
setkey(adults, survey_event, Tank_num)
merged_adults = merge(merge(adults, adult_ze, all.x=T), tad_ze, all.x=T)

##############################################################################
# Formatting the data latent zoospore pool model in Stan ##################### 
##############################################################################

merged_adults = merged_adults[order(Indiv_id, Day), ]
merged_adults$time = c(0, diff(merged_adults$Day))
merged_adults[time < 0, "time"] = 0 # All negatives are really no time change

# Creating empty lists to hold the various predictors that are going to passed
# into Stan as data.
num_indivs = array(NA, dim=length(unique(merged_adults$Tank_num)))
frog_ze_arrays = list()
time_arrays = list()
tad_ze_arrays = list()
density_arrays = list()
frequency_arrays = list()
prob_inf_arrays = list()
swab_id_arrays = list()


# Some helpful function for dealing with the lack of ragged arrays in Stan
growth_array = function(x, size){
    num_grow = size - length(x)
    return(c(x, array(NA, dim=num_grow)))
}

fill_array = function(tarray, value){
    comp_array = array(value, dim=c(16 - nrow(tarray), 6))
    comp_array = rbind(tarray, comp_array)
    return(comp_array)
}


for(tank in sort(unique(merged_adults$Tank_num))){ # Loop through Tanks

    tank_dat = merged_adults[(Tank_num == tank), ]
    swab_events = tank_dat[, list(swabs=length(ZE)), by=c("Indiv_num")]

    setkey(swab_events, Indiv_num)
    setkey(tank_dat, Indiv_num)

    tank_swabs = merge(tank_dat, swab_events, all.x=T)
    tank_swabs = tank_swabs[(swabs >= 6) & (!is.na(ZE)), ]  

    ids = unique(tank_swabs$Indiv_num)
    num_indivs[tank] = length(ids)

    # If there are any frogs in the tank
    if(length(ids) > 0){

        frog_ze_arrays[[tank]] = lapply(ids, function(x) tank_swabs[Indiv_num == x, tot_ze])[[1]]
        tad_ze_arrays[[tank]] = lapply(ids, function(x) tank_swabs[Indiv_num == x, tot_ze_tad])[[1]]
        density_arrays[[tank]] = lapply(ids, function(x) tank_swabs[Indiv_num == x, num_inf])[[1]]
        frequency_arrays[[tank]] = lapply(ids, function(x) tank_swabs[Indiv_num == x, num_inf] / tank_swabs[Indiv_num == x, InitialFrogs])[[1]]
        prob_inf_arrays[[tank]] = lapply(ids, function(x) as.integer(tank_swabs[Indiv_num == x, ZE] > 0))
        swab_id_arrays[[tank]] = lapply(ids, function(x) as.character(tank_swabs[Indiv_num == x, swab_id]))
        time_arrays[[tank]] = lapply(ids, function(x) tank_swabs[Indiv_num == x, time])[[1]]

        # Extend observations that aren't long enough and set to -1
        inf_array = do.call(rbind, lapply(prob_inf_arrays[[tank]], growth_array, 6))
        inf_array[is.na(inf_array)] = -1
        prob_inf_arrays[[tank]] = inf_array

        tswab_array = do.call(rbind, lapply(swab_id_arrays[[tank]], growth_array, 6))
        swab_id_arrays[[tank]] = tswab_array

        # ttime_array = do.call(rbind, lapply(time_arrays[[tank]], growth_array, 6))
        # ttime_array[is.na(ttime_array)] = -1
        # time_arrays[[tank]] = ttime_array

        # Extend arrays with fewer than 16 individuals
        if(nrow(inf_array) < 16){

            prob_inf_arrays[[tank]] = fill_array(inf_array, -1)
            swab_id_arrays[[tank]] = fill_array(tswab_array, NA) 
            # time_arrays[[tank]] = fill_array(ttime_array, -1) 
        }
    }   
} # End tank loop

# Total number of data points that are contributing to the log likelihood
counter = 0
all_ids = list()
for(t in 1:length(prob_inf_arrays)){
    for(i in 1:6){
        if(i != 1){
            for(j in 1:16){
                if(j <= num_indivs[t]){
                    if(prob_inf_arrays[[t]][j, i - 1] == 0){
                        counter = counter + 1
                        all_ids[[counter]] = swab_id_arrays[[t]][j, i]
                    }
                }

            }
        }
    }
}

all_ids = unlist(all_ids)


##############################################################################
# Fitting the latent zoospore pool model ##################################### 
##############################################################################

dic = function(mod){
    # Calculate the DIC of a Stan model

    # Extract the full deviance
    dev_full = extract(mod)$dev
    dev_bar = mean(dev_full)
    pd_alt = 2*var(dev_full) # Alternative effective parameters given by Gelman
    return(dev_bar + pd_alt)
}

colVars <- function(a) {
    # Get columns variances.
    n <- dim(a)[[1]];
    c <- dim(a)[[2]];
    return(.colMeans(((a - matrix(.colMeans(a, n, c), nrow=n, ncol=c,
                        byrow=TRUE))^2), n, c) * n / (n - 1))
}

waic <- function(stanfit){
    # Calculate WAIC from a Stan model. From Vehtari and Gelman 2014

    log_lik <- extract(stanfit, "log_lik")$log_lik

    dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))

    S <- nrow(log_lik)
    n <- ncol(log_lik)

    lpd <- log(colMeans(exp(log_lik)))

    p_waic <- colVars(log_lik)
    elpd_waic <- lpd - p_waic
    waic <- -2*elpd_waic

    loo_weights_raw <- 1/exp(log_lik-max(log_lik))
    loo_weights_normalized <- loo_weights_raw / matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
    loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))

    elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/colMeans(loo_weights_regularized))
    p_loo <- lpd - elpd_loo

    pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
    total <- colSums(pointwise)
    se <- sqrt(n*colVars(pointwise))

    return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
        p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
        pointwise=pointwise, total=total, se=se))
}

##############################################################################

sigmas = c(1) # Default sigma value chosen based on sensitivity analysis

save_plots = FALSE # Set to TRUE to save the DIC plots

# Uncomment and set save_plots = TRUE for DIC and WAIC plots
# Sensitivity analysis of latent zoospore model to sigma
# sigmas = c(0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)

waic_results = list()
dic_results = list()
b_zpool_results = list()
b_density_results = list()

for(i in 1:length(sigmas)){ # Loops through different values of sigma

    sigma =  sigmas[i]

    pomp_dd = list(N=6,
                M=16,
                T=16,
                num_frogs=num_indivs,
                frog_zes=frog_ze_arrays,
                tad_zes=tad_ze_arrays,
                prob_inf=prob_inf_arrays,
                density=density_arrays,
                time=time_arrays,
                init_tads=2000,
                D=counter,
                sigma=sigma)

    pomp_mod_dd = stan("stan_files/pomp_zpool_with_decay.stan", data=pomp_dd, 
                                   warmup=2000, iter=6000, chains=3, cores=3)

    pomp_mod_zp = stan("stan_files/pomp_just_zpool.stan", data=pomp_dd, 
                                    warmup=2000, iter=6000, chains=3, cores=3)

    pomp_fd = list(N=6,
                M=16,
                T=16,
                num_frogs=num_indivs,
                frog_zes=frog_ze_arrays,
                tad_zes=tad_ze_arrays,
                prob_inf=prob_inf_arrays,
                time=time_arrays,
                density=frequency_arrays,
                init_tads=2000,
                D=counter,
                sigma=sigma)

    pomp_mod_fd = stan("stan_files/pomp_zpool_with_decay.stan", data=pomp_fd, 
                        warmup=2000, iter=6000, chains=3, cores=3)

    # Store the models where sigma = 1 in a list. Used in other analyses
    if(sigma == 1){
        models_in_manuscript = list(pomp_mod_dd=pomp_mod_dd, 
                                    pomp_mod_zp=pomp_mod_zp, 
                                    pomp_mod_fd=pomp_mod_fd)
    }

    # Instead of host density, let's use zoospore density in the density-dependent
    # componentns

    # pomp_dd_load = list(N=6,
    #             M=16,
    #             T=16,
    #             num_frogs=num_indivs,
    #             frog_zes=frog_ze_arrays,
    #             tad_zes=tad_ze_arrays,
    #             prob_inf=prob_inf_arrays,
    #             density=lapply(frog_ze_arrays, function(x) log(x + 1)),
    #             time=time_arrays,
    #             init_tads=2000,
    #             D=counter,
    #             sigma=sigma)

    # pomp_mod_dd_load = stan("stan_files/pomp_zpool_with_decay.stan", data=pomp_dd_load, 
    #                                warmup=2000, iter=6000, chains=3, cores=3)

    ## Do some model comparison

    # WAIC comparisons
    waic_results[[i]] = sapply(list(pomp_mod_dd, pomp_mod_fd, pomp_mod_zp), 
                                    function(mod) waic(mod)$total[1])

    # DIC comparisons
    dic_results[[i]] = sapply(list(pomp_mod_dd, pomp_mod_fd, pomp_mod_zp), dic)

    # Extract coefficients
    b_density_results[[i]] = sapply(list(pomp_mod_dd, pomp_mod_fd, pomp_mod_zp),
                                    function(mod) quantile(extract(mod)$b_density, 
                                                c(0.025, 0.5, 0.975)))

    b_zpool_results[[i]] = sapply(list(pomp_mod_dd, pomp_mod_fd, pomp_mod_zp),
                                 function(mod) quantile(extract(mod)$b_zpool, 
                                                c(0.025, 0.5, 0.975)))

} # End sigma loop


if(save_plots){

    ## Plot the results from the sensitivity analysis

    # WAIC results
    waic_array = do.call(rbind, waic_results)
    colnames(waic_array) = c("density dependent", "frequency dependent", "zoospore pool")
    waic_plot = melt(waic_array)
    colnames(waic_plot) = c("sigma", "model", "WAIC")
    waic_plot$sigma = rep(sigmas, 3)

    waic_gg = ggplot(waic_plot, aes(x=sigma, y=WAIC, color=model)) + 
                            geom_line() + geom_point() +
                            theme_bw() +
                            xlab(expression(sigma)) +
                            annotate("text", x=0.1, y=350, label="A.")

    # DIC results
    dic_array = do.call(rbind, dic_results)
    colnames(dic_array) = c("density dependent", "frequency dependent", "zoospore pool")
    dic_plot = melt(dic_array)
    colnames(dic_plot) = c("sigma", "model", "DIC")
    dic_plot$sigma = rep(sigmas, 3)

    dic_gg = ggplot(dic_plot, aes(x=sigma, y=DIC, color=model)) + 
                            geom_line() + geom_point() +
                            theme_bw() +
                            xlab(expression(sigma)) +
                            annotate("text", x=0.1, y=1500, label="B.")

    # Combine these plots
    pdf("../../results/waic_dic_compare.pdf", width=8, height=6)
    grid.arrange(waic_gg, dic_gg)
    dev.off()


    # Coefficent plots

    medians = t(sapply(b_zpool_results, function(x) x[2, ]))
    colnames(medians) = c("density dependent", "frequency dependent", "zoospore pool")
    medians = melt(medians)
    colnames(medians) = c("sigma", "model", "median")
    medians$sigma = rep(sigmas, 3)

    uppers = melt(t(sapply(b_zpool_results, function(x) x[3, ])))
    lowers = melt(t(sapply(b_zpool_results, function(x) x[1, ])))

    medians$lower = lowers$value
    medians$upper = uppers$value

    zpool_gg = ggplot(medians, aes(x=sigma, y=log(median), color=model)) + 
                        geom_line() + geom_point() +
                        geom_errorbar(aes(ymin=log(lower), ymax=log(upper), color=model), 
                                    alpha=0.8) +
                        theme_bw() + 
                        ylab(expression(paste("ln(", beta[0], ": zoospore pool coefficient)", sep=""))) +
                        xlab(expression(sigma)) +
                        annotate("text", x=0, y=-2, label="A.")


    medians = t(sapply(b_density_results, function(x) x[2, ]))
    colnames(medians) = c("density dependent", "frequency dependent", "zoospore pool")
    medians = melt(medians)
    colnames(medians) = c("sigma", "model", "median")
    medians$sigma = rep(sigmas, 3)

    uppers = melt(t(sapply(b_density_results, function(x) x[3, ])))
    lowers = melt(t(sapply(b_density_results, function(x) x[1, ])))

    medians$lower = lowers$value
    medians$upper = uppers$value

    density_gg = ggplot(medians, aes(x=sigma, y=median, color=model)) + 
                        geom_line() + geom_point() +
                        geom_errorbar(aes(ymin=lower, ymax=upper, color=model), 
                                    alpha=0.8) +
                        theme_bw() + ylab(expression(paste(beta[1], ": density/frequency coefficient", sep=""))) +
                        xlab(expression(sigma)) +
                        annotate("text", x=0, y=0.875, label="B.")


    pdf("../../results/coef_plots.pdf", width=8, height=6)
    grid.arrange(zpool_gg, density_gg)
    dev.off()

}


## Some diagnostic plots of the estimated, unobserved zoospore pool

# zpool = extract(pomp_mod_dd)$zpool
# median_zpool = apply(zpool, c(2, 3), median)
# lower_zpool = apply(zpool, c(2, 3), quantile, c(0.025))
# upper_zpool = apply(zpool, c(2, 3), quantile, c(0.975))

# melted_median = melt(median_zpool)
# melted_lower = melt(lower_zpool)
# melted_upper = melt(upper_zpool)

# colnames(melted_median) = c("tank", "time", "zspore")
# colnames(melted_lower) = c("tank", "time", "lower")
# colnames(melted_upper) = c("tank", "time", "upper")
# melted_median = cbind(melted_median, melted_lower$lower, melted_upper$upper)
# colnames(melted_median) = c("tank", "time", "median", "lower", "upper")


# ggplot(data=melted_median, aes(x=time, y=log(median + 1))) + geom_line() +
#         geom_ribbon(aes(ymin=log(lower + 1), ymax=log(upper + 1)), alpha=0.25) +
#         facet_wrap(~as.factor(tank))