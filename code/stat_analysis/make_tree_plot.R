# Description
# -----------
#
# Making a plot of pruned regression trees for the sensitivity analyses that 
# explore the global sensitivity of extinction to transmission, resistance,
# and tolerance.  Make pruned regression trees for three different
# transmission functions.
#
# Author: Mark Wilber

library(tree)
library(plyr)
library(rpart)
library(rpart.plot)
set.seed(1)

# Different transmission functions on which to run regression tree analysis
trans_fxns = c("zoospore_pool_infection_const",
            "zoospore_pool_infection_fd",
            "zoospore_pool_infection_dd")

# Labels for nodes
labels = c("s(x), intercept",
          "s(x), size",
          "G(x', x), intercept",
          "G(x', x), temperature",
          "G(x', x), size",
          "G(x', x), variance",
          "G(x', x), variance effect",
          "G0(x'), intercept",
          "G0(x'), temperature", 
          "G0(x'), variance",
          "G0(x'), variance effect",
          "l(x), intercept",
          "l(x), temperature",
          "l(x), size", 
          "Transmission, zoospore pool", "Transmission, density/frequency")

sigmas = c(0.3) #c(0.1, 0.2, 0.3, 0.4, 0.5)
for(sigma in sigmas){

  for(trans in trans_fxns){

      dat = read.csv(
          paste("/Users/mqwilber/Repos/density_dependent_ipm/results/sens_results_", 
          trans, "_", sigma, ".csv", sep=""))

      # Uncomment and comment above for bayesian sensivity analysis data
      # dat = read.csv(
      #     paste("/Users/mqwilber/Repos/density_dependent_ipm/results/sens_results_", 
      #     trans, "_bayesian.csv", sep=""))
      dat = dat[, -1]
      dat$y = as.factor(dat$y)
      dat$y = revalue(dat$y, c("1"="extinct", "0"="persist"))

      # Drop extra parameter if necessary
      if(trans == "zoospore_pool_infection_const")
        colnames(dat)[1:(length(labels) - 1)] = labels[-length(labels)]
      else
        colnames(dat)[1:length(labels)] = labels

      # Fit the regression tree using rpart
      treemod = rpart(y ~ ., data=dat, cp=0.01, parms=list(split="gini"))

      # Get the minimum CV Error for all alphas
      treecv = printcp(treemod)
      min_error = which.min(treecv[, "xerror"])
      min_cp = treecv[min_error, "CP"]
      prune_treemod = prune.rpart(treemod, min_cp)

      pdf(paste("../../results/tree_plot_", trans, "_", sigma, ".pdf", sep=""))

      # Uncomment and comment above for saving bayesian sens. results
      # pdf(paste("../../results/tree_plot_", trans, "_bayesian.pdf", sep=""))
      rpart.plot(prune_treemod, clip.right.labs=FALSE, varlen=0, extra=0)

      dev.off()

  }
}