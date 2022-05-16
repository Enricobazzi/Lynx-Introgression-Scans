import sys
import os
import numpy
import dadi
from dadi import Numerics, PhiManip, Integration, Misc
from dadi.Spectrum_mod import Spectrum
import pylab
import matplotlib.pyplot as plt
import my_models
import argparse
import pandas

'''
plot the comparison between data and model with optimized parameters

Usage: python dadi_plot_2D_comp_multinom.py <model_name>
'''

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
# Import data to create DATA joint-site frequency spectrum
#===========================================================================

#**************
# vcf file
datafile = 'data/dummy_dani_variants.vcf.gz'
# pop file
popfile = 'data/dummy_popfile.txt'

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict_vcf(datafile, popfile)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids = pops.split("_")

#**************
#projection sizes, in ALLELES not individuals
proj = [38,26]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#===========================================================================
# create MODEL joint-site frequency spectrum from optimized parameters
#===========================================================================
# parameters - optimized for maximum log-likelihood
all_best_lls = []

for i in range(1, 21):
    
    model_file_name = ["tables/", pops, ".", chosen_model, ".optimized.r_", str(i), ".txt"]
    model_file_name = ''.join(model_file_name)
    if os.path.exists(model_file_name):
        model_optimized_allrounds = pandas.read_csv(model_file_name, sep='\t')
        ll_col = model_optimized_allrounds.loc[:,"log-likelihood"]
        best_ll = model_optimized_allrounds.iloc[ll_col.idxmax(),]
        all_best_lls.append(best_ll)

all_best_lls_df = pandas.concat(all_best_lls, axis=1).T
best_ll = all_best_lls_df.sort_values('log-likelihood', ascending=False)[:1]
param_string = best_ll.iloc[0,6].split(",")
params = [float(l) for l in param_string]

print("for", chosen_model, "best repetition was", best_ll.iloc[0,1])
print("with a log-likelihood of", best_ll.iloc[0,2], "and the following parameters:")
if len(params) == 9:
    param_def = ["Tsplit", "Tbot", "iber_a", "iber_pr", "eura_a", "eura_pr", "m", "m_12", "m_21"]
if len(params) == 13:
    param_def = ["Tsplit", "Tbot", "Tbot_a", "iber_a", "iber_pr_a", "iber_pr", "eura_a", "eura_pr", "m", "ma_12", "ma_21", "m_12", "m_21"]
for i in range(0, len(params)):
        print(param_def[i], ":", params[i])
ns = proj
pts = 60


if chosen_model == "model_1_a":
	model = my_models.model_1_a(params=params, ns=ns, pts=pts)

elif chosen_model == "model_1_b":
	model = my_models.model_1_b(params=params, ns=ns, pts=pts)

elif chosen_model == "model_2_a":
	model = my_models.model_2_a(params=params, ns=ns, pts=pts)

elif chosen_model == "model_2_b":
	model = my_models.model_2_b(params=params, ns=ns, pts=pts)

elif chosen_model == "model_2_c":
	model = my_models.model_2_c(params=params, ns=ns, pts=pts)

elif chosen_model == "model_2_d":
	model = my_models.model_2_d(params=params, ns=ns, pts=pts)

model_folded = model.fold()


#===========================================================================
# create MODEL joint-site frequency spectrum from optimized parameters
#===========================================================================
plot_file_name = ["plots/", chosen_model, ".2D_comp_multinom.pdf"]
plot_file_name = ''.join(plot_file_name)

fig = pylab.figure(1)
fig.clear()
dadi.Plotting.plot_2d_comp_multinom(model, fs, pop_ids=pop_ids)
fig.savefig(plot_file_name)


