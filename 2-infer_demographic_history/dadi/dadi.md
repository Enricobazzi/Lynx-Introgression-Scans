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

After defining the models to test, we can proceed to optimize the parameters, so we can find the values that best fit the data to the model. The strategy we adopt for this task is based on the tool [dadi_pipeline](https://github.com/dportik/dadi_pipeline) ([Portik et al. 2017](https://onlinelibrary.wiley.com/doi/10.1111/mec.14266)):

Using a customized version of the dadi_pipeline script [dadi_run_models_optimize.py](./dadi_run_models_optimize.py), select a model (first argument) and a pair of populations  and we use the function Optimize_Routine defined in [Optimize_Functions.py](./Optimize_Functions.py) (from dadi_pipeline) to perform 4 consecutive rounds of optimizations for each model. Each round consisted of multiple replicates, using the parameter estimates from the best scoring replicate (highest log-likelihood) to seed searches in the following round. We used the default settings in dadi_pipeline for each round (replicates = 10, 20, 30, 40; maxiter = 3, 5, 10, 15; fold = 3, 2, 2, 1), and we optimized parameters using the Nelder-Mead method ([optimize_log_fmin](https://dadi.readthedocs.io/en/latest/api/dadi/Inference.html#dadi.Inference.optimize_log_fmin)).

To find which [projection](https://dadi.readthedocs.io/en/latest/user-guide/frequently-asked-questions/#1-im-projecting-my-data-down-to-a-smaller-frequency-spectrum-what-sample-sizes-should-i-project-down-to) (downsampling) of our data would yield the highest number of SNPs for a particular population pair, I ran a python [script](./best_projections.py) and used its output to manually edit the optimization [script](./dadi_run_models_optimize.py).

The optimization script [dadi_run_models_optimize.py](./dadi_run_models_optimize.py) is then run 20 times for each model, to check for convergence of the optimization process. For that we use a custom [script](./cesga_dadi_opti.sh) that will run one optimization round, for one model and one population pair:

```{bash}
# for lpa-wel
for model in model_1_a model_1_b model_2_a model_2_b model_2_c
 do
  for i in {1..20}
   do
    echo "sbatch ${model} of lpa-wel number ${i}"
    sbatch -t 04-00:00 -n 1 -c 1 --mem=60GB cesga_dadi_opti.sh ${model} lpa-wel ${i}
  done
done
```

A different table will be generated containing the results from each the 20 runs with the format dadi_<model>_<pop_pair>_<n> indicating the model, population pair and the repetition number.

To create a single table with the results of all of the optimizations ordered by likelihood I ran an R [script](./order_optimization_results.R). This table will be used to interpret and visualize the dadi results as well as part of the input for the simulations, as it will indicate optimized parameter values for the best scoring model.

## Interpreting and displaying results

To visualize the way likelihood changed throughout the optimization process I have wrote an R [script](./plot_optimizations.R) that takes the table of the ordered optimization results (produced [here](./order_optimization_results.R)) and output a few plots that reflect different aspects of the optimization process:
 - the *lolli_plot* will have the log-likelihood of each all the 20 runs' 100 replicates, to show the trajectory of the log-likelihood across each run
 - *lolli_fourth_panels_plot* will again plot the log-likelihood of the 20 runs, but only the fourth round of each optimization to check if the optimization has converged
 - *bar_top10_lolli_plot* will show the top 10 rounds of highest log-likelihood, as a bar plot color coded by run number, to check from how many different optimization rounds the best scoring optimizations come from. 

Additionally, we want to see how the 2D-SFS (see dadi [documentation](https://dadi.readthedocs.io/en/latest/user-guide/plotting/)) of our data compares to the one generated by our model, using the optimized parameters. The dadi function [dadi.Plotting.plot_2d_comp_multinom](https://dadi.readthedocs.io/en/latest/api/dadi/Plotting.html#dadi.Plotting.plot_2d_comp_multinom) can be used to produce a plot that will show:
 - the 2D-SFS of the data
 - the 2D-SFS of a given model with a given set of parameter values
 - a joint spectrum showing the Anscombe residuals of the data-model comparison
 - an histogram of the Anscombe residuals

The python script [plot_functions](./plot_functions.py) contains the functions that load data into dadi from a VCF, extract parameters from the table of ordered optimization results, build a model using a set of parameters and plot the 2D comparisons. Using an interactive jupyter notebook [script](./plotting_interactively.ipynb) the comparisons between the data and the 4 best scoring parameters for each model were plotted.

## Processing the results into real world parameters

To convert the parameters obtained in the optimization runs into real world parameters of time, population size and migration rates I used custom python script [process_dadi_results](./process_dadi_results.py), which contains both the functions to mathematically convert the parameters and a "main" function that takes the table of ordered optimization results and prints a table with the converted results of the best 4 optimization rounds of each model.
