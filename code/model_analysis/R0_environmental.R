## Description
## ------------
## Script for calculating R0 with and without environmental transmission from 
## the R. muscosa-Bd IPM. Also makes the necessary plots
## 
## Author: Mark Wilber

# Load the parameters
source("IPM_functions_for_R.R")
library(ggplot2)
library(plyr)
library(yaml)
library(RColorBrewer)
library(GGally)
library(gridExtra)

calc_R0 = function(params, pfull, y, bins, s_0, h, S_init){
    # Function calculates R0 with and without an environmental reservoir
    # assuming density dependent transmission
    #
    # Parameters
    # -----------
    # y : vector of midpoints for load classes
    # params : list of temperature-dependent params as calculated by set_temp_params
    # pfull : list of all parameters from "../../results/bayesian_parameter_estimates.rds"
    #
    # Returns
    # -------
    # : c(R0 for host + environmental transmission, R0 just for host transmission)

    # First compute transition probabilities
    G = h * outer(y, y, g_xpx, params=params)
    S = s_x(y, params=params)
    R = r_x(y, params=params)

    U_host = G %*% diag(S * (1 - R)) # The transition matrix
    U_full = get_full_P(U_host, c(log(exp(y) + 1), params$zspore_surv), 
                           c(rep(0, bins), params$zspore_surv), 
                           add_top=FALSE)

    hostcol = c(s_0 * h * g0_x(y, params) * pfull$zdd_slope * S_init, 0)
    zcol = c(s_0 * h * g0_x(y, params) * pfull$zdd_intercept * S_init, 0)
    host_matrix = sapply(1:bins, function(x) hostcol)

    M_host = host_matrix[1:bins, ]
    M_full = cbind(host_matrix, zcol)

    # Calculate R0 with and without environmental zoospore pool
    R_mat_full = M_full %*% solve(diag(bins + 1) - U_full)
    R_mat_host = M_host %*% solve(diag(bins) - U_host)

    R0_full = max(Re(eigen(R_mat_full)$values))
    R0_host = max(Re(eigen(R_mat_host)$values))

    return(c(R0_full, R0_host))

}


# Read and set params
pfull = readRDS("../../results/bayesian_parameter_estimates.rds")
pother = yaml.load_file("model_parameters.yml")

############################
# Parameters for IPM model #
############################

linear = TRUE
min_size = -5
max_size = 18
bins = 50

matrix_vals = set_discretized_values(min_size, max_size, bins)
y = matrix_vals$y # Centers of discretized matrix
h = matrix_vals$h # Delta between bin edges

S_init_vals = seq(1, 16, len=100) # Initial number of susceptible hosts
s_0 = 1 # Survival probability is 1
temps = seq(4, 20, len=100) # Initial temperatures are 1

# Arrays for storing results
R0s_full = array(NA, dim=c(length(temps), length(S_init_vals))) # With host and environ
R0s_host = array(NA, dim=c(length(temps), length(S_init_vals))) # With just host

################################################################
## Step 2: Calculate R0 with and without an environmental reservoir     ##
################################################################

for(j in 1:length(S_init_vals)){

    for(i in 1:length(temps)){

        desired_temp = temps[i]
        params = set_temp_params(desired_temp, pfull, linear)
        #params$zspore_surv = 0

        R0_vals = calc_R0(params, pfull, y, bins, s_0, h, S_init_vals[j]) 
        R0s_full[i, j] = R0_vals[1] 
        R0s_host[i, j] = R0_vals[2] 

    }
}


plot_list = list()
data_list = list(R0s_host, R0s_full)
trans_type = c("W/out zoospore pool", "With zoospore pool")

for(m in 1:length(data_list)){

    R0s_exact = data_list[[m]]
    rownames(R0s_exact) = temps
    colnames(R0s_exact) = S_init_vals

    R0_data = melt(R0s_exact)
    colnames(R0_data) = c("temp", "S_init", "R0")
    R0_data$prob = 1 - 1 / R0_data$R0
    R0_data$prob[R0_data$prob < 0] = 0
    R0_data$trans_type = trans_type[m]

    plot_list[[m]] = R0_data

}

plot_data = do.call(rbind, plot_list)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
inv_plot = ggplot(plot_data, aes(x=temp, y=S_init, fill=prob)) + geom_tile() +
             scale_fill_gradientn(colours = myPalette(100), name="Invasion probability") +
             theme_bw() + facet_wrap(~ trans_type) +
             xlab("Temperature") + ylab(expression("Initial susceptible frog density, " ~ m^-3)) +
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
             geom_segment(aes(x=4,xend=20,y=4,yend=4), color="black") +
             geom_segment(aes(x=15,xend=15,y=1,yend=16), color="black", linetype="dashed")

ggsave("../../results/R0_heat_map.pdf", inv_plot, width=8, height=5)


################################################
## Step 3: Calculate a slice with uncertainty fixed density ##
################################################

set.seed(2)
draws = 500
S_init = 4
temps = seq(4, 20, len=20)

R0s_full_uncertain = array(NA, dim=c(length(temps), draws))
R0s_host_uncertain = array(NA, dim=c(length(temps), draws))

for(j in 1:draws){

    rand_params = bayesian_params(pfull)

    print(paste("Working on simulation", j))

    for(i in 1:length(temps)){

        desired_temp = temps[i]
        params = set_temp_params(desired_temp, rand_params, linear)
        #params$zspore_surv = 0

        R0_vals = calc_R0(params, rand_params, y, bins, s_0, h, S_init) 
        R0s_full_uncertain[i, j] = R0_vals[1] 
        R0s_host_uncertain[i, j] = R0_vals[2] 

    }
}

R0s_full_data = as.data.frame(t(apply(R0s_full_uncertain, 1, quantile, c(0.025, 0.5, 0.975))))
R0s_host_data = as.data.frame(t(apply(R0s_host_uncertain, 1, quantile, c(0.025, 0.5, 0.975))))

# Rename the columns
colnames(R0s_full_data) = c("lower", "median", "upper")
colnames(R0s_host_data) = c("lower", "median", "upper")

# Plot environmental R0
tplot = ggplot(data=R0s_full_data, aes(x=temps, y=log(median))) +
                geom_line(color="red3", size=1) +
                geom_ribbon(aes(ymin=log(lower), ymax=log(upper)), alpha=0.1) +
                geom_hline(yintercept=0, linetype="dashed") +
                theme_bw() + 
                ylab(expression(paste("ln(", R[0], ")", sep=""))) + 
                xlab(expression(paste("Temperature,", degree*C))) +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text=element_text(size=13),
                      axis.title=element_text(size=15))

# Plot host R0
tplot = tplot + geom_line(data=R0s_host_data, color="blue3", size=1, aes(x=temps, y=log(median))) + 
        geom_ribbon(data=R0s_host_data, aes(ymin=log(lower), ymax=log(upper)), alpha=0.1) + 
        annotate("text", x=9, y=2,  label="With zoopore pool") +
        annotate("text", x=15, y=1,  label="W/out zoopore pool")

ggsave("../../results/R0_density_slice.pdf", tplot, width=8, height=6)

################################################
## Step 4: Calculate a slice with uncertainty fixed temp ##
################################################

set.seed(2)
draws = 500
S_init_vals = seq(1, 16, len=50) # Initial number of susceptible hosts
temp = 15

R0s_full_uncertain = array(NA, dim=c(length(S_init_vals), draws))
R0s_host_uncertain = array(NA, dim=c(length(S_init_vals), draws))

for(j in 1:draws){

    rand_params = bayesian_params(pfull)

    print(paste("Working on simulation", j))

    for(i in 1:length(S_init_vals)){

        desired_temp = temp
        params = set_temp_params(desired_temp, rand_params, linear)
        #params$zspore_surv = 0

        R0_vals = calc_R0(params, rand_params, y, bins, s_0, h, S_init_vals[i]) 
        R0s_full_uncertain[i, j] = R0_vals[1] 
        R0s_host_uncertain[i, j] = R0_vals[2] 

    }
}

R0s_full_data = as.data.frame(t(apply(R0s_full_uncertain, 1, quantile, c(0.025, 0.5, 0.975))))
R0s_host_data = as.data.frame(t(apply(R0s_host_uncertain, 1, quantile, c(0.025, 0.5, 0.975))))

# Rename the columns
colnames(R0s_full_data) = c("lower", "median", "upper")
colnames(R0s_host_data) = c("lower", "median", "upper")

# Plot environmental R0
tplot_dens = ggplot(data=R0s_full_data, aes(x=S_init_vals, y=log(median))) +
                geom_line(color="red3", size=1) +
                geom_ribbon(aes(ymin=log(lower), ymax=log(upper)), alpha=0.1) +
                geom_hline(yintercept=0, linetype="dashed") +
                theme_bw() + 
                ylab(expression(paste("ln(", R[0], ")", sep=""))) + 
                xlab(expression("Initial susceptible frog density," ~ m^-3)) +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text=element_text(size=13),
                      axis.title=element_text(size=15))

# Plot host R0
tplot_dens = tplot_dens + geom_line(data=R0s_host_data, color="blue3", size=1, aes(x=S_init_vals, y=log(median))) + 
        geom_ribbon(data=R0s_host_data, aes(ymin=log(lower), ymax=log(upper)), alpha=0.1) + 
        annotate("text", x=8, y=3.6,  label="With zoopore pool") +
        annotate("text", x=12, y=2.1,  label="W/out zoopore pool")

ggsave("../../results/R0_temperature_slice.pdf", tplot_dens, width=8, height=6)

