---
title: "dadi"
author: "Enrico"
date: "5/16/2022"
output: html_document
---

## Defining models to test

The first step of running dadi is defining the demographic models you want to test.

To do so, you have to write them in the form of python functions. Most info on how to do that can be found at:

[dadi - home](https://dadi.readthedocs.io/en/latest/)

[dadi - model specifying](https://dadi.readthedocs.io/en/latest/user-guide/specifying-a-model/)

For our study, we have decided to test 6 different demographic scenarios for our two populations:
  - model_1_a : iberian sudden bottleneck, eurasian exponential decline
  - model_1_b : iberian exponential decline, eurasian exponential decline
  - model_2_a : iberian double sudden bottleneck, eurasian exponential decline
  - model_2_b : iberian exponential decline followed by sudden bottleneck, eurasian exponential decline
  - model_2_c : iberian sudden bottleneck followed by exponential decline, eurasian exponential decline
  - model_2_d : iberian double exponential decline, eurasian exponential decline

In all of our models, population split is followed by a period of symmetrical gene-flow. After population declines, gene-flow is allowed to be asymmetrical.

The dadi python functions for these models can be found in [my_models.py](./my_models.py)

## Parameter optimization and model likelihood calculation

After defining the models to test, we can proceed to optimize the parameters, so we can find the values that best fit the data to the model. The strategy we adopt for this task is based on the tool [dadi_pipeline](https://github.com/dportik/dadi_pipeline):

Using a customized version of the dadi_pipeline script [dadi_run_models_optimize.py](./dadi_run_models_optimize.py), select a model (first argument) and a pair of populations  and we use the function Optimize_Routine defined in [Optimize_Functions.py](./Optimize_Functions.py) (from dadi_pipeline) to run 4 rounds of optimization for the selected model, using the [Nelder-Mead](https://dadi.readthedocs.io/en/latest/api/dadi/Inference.html#dadi.Inference.optimize_log_fmin) method of optimization.

Initial parameters are first selected at random within the bounds defined in [dadi_run_models_optimize.py](./dadi_run_models_optimize.py), and then each replicate (see below) will use last replicate's optimized parameters as the initial parameters.

The first round will have 10 consecutive replicates, each of which will use a maximum of 3 iterations each, and [perturbing](https://dadi.readthedocs.io/en/latest/api/dadi/Misc.html#dadi.Misc.perturb_params) the initial parameters 3-fold.

The second round will have 20 consecutive replicates, each of which will use a maximum of 5 iterations each, and [perturbing](https://dadi.readthedocs.io/en/latest/api/dadi/Misc.html#dadi.Misc.perturb_params) the initial parameters 2-fold.

The third round will have 30 consecutive replicates, each of which will use a maximum of 10 iterations each, and perturbing the initial parameters ([dadi.Misc.perturb_params](https://dadi.readthedocs.io/en/latest/api/dadi/Misc.html#dadi.Misc.perturb_params)) 2-fold.

The fourth round will have 40 consecutive replicates, each of which will use a maximum of 15 iterations each, and [perturbing](https://dadi.readthedocs.io/en/latest/api/dadi/Misc.html#dadi.Misc.perturb_params) the initial parameters 1-fold.

**Before running dadi_run_models_optimize.py, make sure the proper input data, projections and grid sizes are defined in the script**

The optimization script [dadi_run_models_optimize.py](./dadi_run_models_optimize.py) is run 20 times for each model, to check for convergence of the optimization process.
```{bash}
# activate dadi environment
conda activate dadi
# run 
for i in {1..20};  do   echo "repetition ${i}";   python dadi_run_models_optimize.py <model_name> <pop_pair>; done
```
The results from the 20 runs will all go into a single table, in which the header will be repeated, and dividing each run's results. To extract each run into its own table we can run the script [split_optimized_table.sh](./split_optimized_table.sh) for the specific model and population pair
```{bash}
./split_optimized_table.sh <model_name> <pop_pair>
```
To visualize the results of the optimization process I have wrote an R [script](./plot_optimizations.R) that takes the model name and the population pair and will output 3 plots:
 - the *lolli_plot* will have the log-likelihood of each all the 20 runs' 100 replicates, color coded by round, to show the trajectory of the log-likelihood across each run
 - *lolli_fourth_panels_plot* will again plot the log-likelihood of the 20 runs, but only from the fourth round of optimization, as a line plot, to check if the optimization has converged
 - *lolli_fourth_together_plot* will show the log-likelihood from the fourth round of all of the 20 runs together, as a dot plot, to check if there are any particular runs that behave differently at convergence
```{bash}
Rscript plot_optimizations.R <model_name> <pop_pair>
```
Additionally, we want to see how the 2D-SFS (see dadi [documentation](https://dadi.readthedocs.io/en/latest/user-guide/plotting/)) of our data compares to the one generated by our model, using the optimized parameters. The dadi function [dadi.Plotting.plot_2d_comp_multinom](https://dadi.readthedocs.io/en/latest/api/dadi/Plotting.html#dadi.Plotting.plot_2d_comp_multinom) can be used to extract the set of parameters with highest log-likelihood and produce a plot that will show:
 - the 2D-SFS from the data
 - the 2D-SFS from the model with highest scoring parameters
 - a joint spectrum showing the Anscombe residuals of the data-model comparison
 - an histogram of the Anscombe residuals
To do so, there is an interactive jupyter notebook [script](./plotting_interactively.ipynb) that can be run running jupyter notebook in the dadi folder. Alternatively the python script [dadi_plot_2D_comp_multinom.py](./dadi_plot_2D_comp_multinom.py) will automatically produce the same plot for the selected model and population pair:
```{bash}
python dadi_plot_2D_comp_multinom.py <model_name> <pop_pair>
```