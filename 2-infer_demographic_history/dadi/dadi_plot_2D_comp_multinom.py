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
args = parser.parse_args()

chosen_model = args.model

#===========================================================================
# Import data to create DATA joint-site frequency spectrum
#===========================================================================

#**************
# vcf file
datafile = 'dummy_dani_variants.vcf.gz'
# pop file
popfile = 'dummy_popfile.txt'

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict_vcf(datafile, popfile)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids = ['sm','ki']

#**************
#projection sizes, in ALLELES not individuals
proj = [19,13]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#===========================================================================
# create MODEL joint-site frequency spectrum from optimized parameters
#===========================================================================
# parameters - optimized for maximum log-likelihood
model_file_name = ["sm_ki.", chosen_model, ".optimized.txt"]
model_file_name = ''.join(model_file_name)
model_optimized_allrounds = pandas.read_csv(model_file_name, sep='\t')
ll_col = model_optimized_allrounds.loc[:,"log-likelihood"]
param_string = model_optimized_allrounds.iloc[ll_col.idxmax(),6].split(",")
params = [float(i) for i in param_string]

ns = proj

pts = 30

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
plot_file_name = [chosen_model, ".2D_comp_multinom.pdf"]
plot_file_name = ''.join(plot_file_name)

fig = pylab.figure(1)
fig.clear()
dadi.Plotting.plot_2d_comp_multinom(model, fs, pop_ids=pop_ids)
fig.savefig(plot_file_name)


