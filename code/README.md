Description of the folders and files contained in the code directory.  All code files also contain a description at the top of the file.

The `code` directory contains two folders:

1. `model_analysis/`
2. `stat_analysis/`


`stat_analysis/`
    
- This folder contains the scripts necessary to fit the various vital rate functions and transmission functions to the laboratory and mesocosm data. The results are saved as `*.rds` files which are used in the scripts iben in `model_analysis/`

`model_analysis/`

- This folder contains the scripts that perform the various simulations and sensitivity analyses on the hybrid model defined in the main text `../docs/ipm_extinction_manusript.pdf`. 


This directory also contains a `run_all.R` script.  This script will replicated all analyses in `..docs/ipm_extinction_*.pdf`.  This `run_all.R` script also shows the order in which the scripts need to be run. For example,
the statistical analysis of the vital rate functions and the transmission functions needs to be run first to that these parameters can be used in the following analyses. 

**WARNING**: With the default number of simulations and cores used in the analysis the `run_all.R` script may take >24 hours to run.  If you want to speed this up in order to play around with how the scripts are working (but, of course, at the cost of not getting a good summary of the stochastic simulations!), make the following changes to reduce the number of simulations.

1. `sensitivity_analysis.R`: Line 50 -> Change to `SIMS = 10`
2. `sensitivity_analysis_bayesian.R`:  Line 50 -> Change to `SIMS = 10`
3. `zdecay_plots.R`: Line 30 -> Change to `sims = 10`
4. `stochastic_multiseason_simulation.R`: Line 40 -> Change to `SIMS = 10`
