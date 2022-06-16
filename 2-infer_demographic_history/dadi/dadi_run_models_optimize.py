'''
Usage: python dadi_run_models_optimize.py <model_name>
This is a modified version of the 'dadi_Run_Optimizations.py' script from the package dadi_pipeline:

https://github.com/dportik/dadi_pipeline

in which we run optimizations for our set of models that have been defined in my_models.py

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary, as well as the  my_models.py script, which
has all the model definitions.

General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates.

 The starting parameters are initially random, but after each round is complete the parameters
 of the best scoring replicate from that round are used to generate perturbed starting parameters
 for the replicates of the subsequent round.
 
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.
 
#########################################################################
##### DECIDE IF GOING WITH RANDOM INITIAL PARAMETERS OR DEFINE THEM #####
#########################################################################

Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688

Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject.

########################################################
##### DECIDE WHAT TO DO ABOUT LINKAGE BETWEEN SNPs #####
########################################################

Citations:
 If you use these scripts or the main diversification models for your work, please
 cite the following publication:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266
 
 If you use the additional diversification models or the island models set please cite 
 the following publication:
    Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K.,
    Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. Sky, sea,
    and forest islands: diversification in the African leaf-folding frog
    Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae).
    Journal of Biogeography 45: 1781-1794. 
    doi: 10.1111/jbi.13365
'''

import sys
import os
import numpy
import dadi
import pylab
import argparse
from datetime import datetime
import Optimize_Functions
import my_models
import subprocess

#===========================================================================
# Arguments
#===========================================================================

parser = argparse.ArgumentParser()
parser.add_argument("model", type=str,
                    help="name of the model as defined in my_models.py")
parser.add_argument("pops", type=str,
                    help="population abbreviations of pair chosen")
args = parser.parse_args()

chosen_model = args.model
pops = args.pops

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************
# vcf file
vcf_prefix = '/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted.filter5.phased.fixed.'
vcf_suffix = '.miss.rd_fil.intergenic.vcf'
datafilename = [vcf_prefix, pops, vcf_suffix]

datafile = ''.join(datafilename)

# pop file
pop_prefix = '/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression/'
pop_suffix = '_dadi_popfile.txt'
popfilename = [pop_prefix, pops, pop_suffix]

popfile = ''.join(popfilename)

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict_vcf(datafile, popfile)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids = pops.split("-")

#**************
#projection sizes, in ALLELES not individuals (see best_projections.py to know how these are calculated)
if pops == "lpa-wel":
    proj_pop1 = 34
    proj_pop2 = 36

if pops == "lpa-eel":
    proj_pop1 = 34
    proj_pop2 = 34

if pops == "lpa-sel":
    proj_pop1 = 34
    proj_pop2 = 22

proj = [proj_pop1, proj_pop2]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

#================================================================================
# Calling external models from the my_models.py script
#================================================================================
'''
We will use a function from the Optimize_Functions.py script for our optimization routines:
 
 Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded=True, 
                        reps=None, maxiters=None, folds=None, in_params=None, 
                        in_upper=None, in_lower=None, param_labels=None, optimizer="log_fmin")
 
   Mandatory Arguments =
    fs:  spectrum object name
    pts: grid size for extrapolation, list of three values
    outfile:  prefix for output naming
    model_name: a label to help label the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or from a separate python model script, 
          ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).
   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order
     optimizer: a string, to select the optimizer. Choices include: "log" (BFGS method), "log_lbfgsb" (L-BFGS-B method), 
                "log_fmin" (Nelder-Mead method), and "log_powell" (Powell's method).

I have set the optimization routine to be the same for eachcmodel using the optional lists below,
which are included as optional arguments for each model. This particular configuration will run 4 rounds as follows:
Round1 - 10 replicates, maxiter = 3, fold = 3
Round2 - 20 replicates, maxiter = 5, fold = 2
Round3 - 30 replicates, maxiter = 10, fold = 2
Round4 - 40 replicates, maxiter = 15, fold = 1
'''


# create a prefix based on the population names to label the output files
# ex. Pop1_Pop2
prefix = "_".join(pop_ids)

#**************
# make sure to define your extrapolation grid size (based on your projections)
pts = [50, 60, 70]

#**************
# set the number of rounds here
rounds = 4

# define the lists for optional arguments
# you can change these to alter the settings of the optimization routine
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]

if (chosen_model == "model_1_a") or (chosen_model == "model_1_b"):
# Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr, m, m_12, m_21 
	in_lower = [1e-4, 1e-4, 1e-4, 0, 1e-4, 0, 0, 0, 0]
	in_upper = [100, 100, 10, 1, 10, 1, 40, 40, 40]
elif (chosen_model == "model_2_a") or (chosen_model == "model_2_b") or (chosen_model == "model_2_c") or (chosen_model == "model_2_d"):
# Tsplit, Tbot2, Tbot1, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21
	in_lower = [1e-4, 1e-4, 1e-4, 1e-4, 0, 0, 1e-4, 0, 0, 0, 0, 0, 0]
	in_upper = [100, 100, 100, 10, 1, 1, 10, 1, 40, 40, 40, 40, 40]
	

#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True

'''
Direct to correct Optimize_Function.Optimize_Routine call based on chosen_model
'''

# Split into two populations, no migration.
if chosen_model == "model_1_a":
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "model_1_a", my_models.model_1_a, rounds, 9, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_lower=in_lower, in_upper=in_upper,
                                        param_labels = "Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr, m, m_12, m_21")
# Split into two populations, with continuous symmetric migration.
elif chosen_model == "model_1_b":
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "model_1_b", my_models.model_1_b, rounds, 9, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_lower=in_lower, in_upper=in_upper,
                                        param_labels = "Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr, m, m_12, m_21")
# Split into two populations, with continuous asymmetric migration.
elif chosen_model == "model_2_a":
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "model_2_a", my_models.model_2_a, rounds, 13, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_lower=in_lower, in_upper=in_upper,
                                        param_labels = "Tsplit, Tbot2, Tbot1, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21")
# Split with continuous symmetric migration, followed by isolation.
elif chosen_model == "model_2_b":
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "model_2_b", my_models.model_2_b, rounds, 13, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_lower=in_lower, in_upper=in_upper,
                                        param_labels = "Tsplit, Tbot2, Tbot1, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21")
# Split with continuous asymmetric migration, followed by isolation.
elif chosen_model == "model_2_c":
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "model_2_c", my_models.model_2_c, rounds, 13, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_lower=in_lower, in_upper=in_upper,
                                        param_labels = "Tsplit, Tbot2, Tbot1, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21")
# Split with no gene flow, followed by period of continuous symmetrical gene flow.
elif chosen_model == "model_2_d":
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "model_2_d", my_models.model_2_d, rounds, 13, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_lower=in_lower, in_upper=in_upper,
                                        param_labels = "Tsplit, Tbot2, Tbot1, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21")
