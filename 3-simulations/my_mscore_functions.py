import numpy as np
import dadi
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

"""
Each of these functions corresponds to one of the models defined by my_models.py, and it will
build a ms core function for that model, indicating the demography to be simulated by ms.

The core function together with the optimized parameters specific to the model, can be plugged
in a dadi function called dadi.Misc.ms_command that will generate the full ms command.

*** Note that I have not included migration or recombination in these commands as of yet ***

"""

##########################################################################################

def model_1_a_mscore(params):
    """
    ms core command corresponding to model_1_a
    model_1_a - one bottleneck, abrupt change
	# population split with different sizes
	# abrupt change in Iberian lynx population at time Tbot
	# exponential change in Eurasian lynx population as function of time since Tbot
    """
    Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr = params
    
	# Bottleneck time:
    Tb = Tbot*Tsplit
    
	# contemporary populations:
    iber_c = iber_a*iber_pr
    eura_c = eura_a*eura_pr
    
	# growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb
    
	# command
    command = "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
              "-eg 0 2 %(eura_alpha)f "\
              "-en %(Tb)f 1 %(iber_a)f "\
              "-en %(Tb)f 2 %(eura_a)f "\
              "-ej %(Tsum)f 2 1"

	# transform dadi parameters to ms parameters
    sub_dict = {'iber_a':iber_a, 'eura_a':eura_a, 'iber_c':iber_c, 'eura_c':eura_c,
                'eura_alpha':2*eura_alpha,
                'Tb':Tb/2, 'Tsum':(Tb+Tsplit)/2}
    
    return command % sub_dict
              
##########################################################################################

def model_1_b_mscore(params):
    """
    ms core command corresponding to model_1_b
    model_1_b - one bottleneck, exponential change
     # population split with different sizes
     # exponential change in Iberian lynx population as function of time since Tbot
     # exponential change in Eurasian lynx population as function of time since Tbot
    """
    Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr = params
    
    # Bottleneck time:
    Tb = Tbot*Tsplit
    
    # contemporary populations:
    iber_c = iber_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha = np.log(iber_c/iber_a)/Tb
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb
    
    # command
    command = "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
              "-eg 0 1 %(iber_alpha)f "\
              "-eg 0 2 %(eura_alpha)f "\
              "-en %(Tb)f 1 %(iber_a)f "\
              "-en %(Tb)f 2 %(eura_a)f "\
              "-ej %(Tsum)f 2 1"
    
    # transform dadi parameters to ms parameters
    sub_dict = {'iber_a':iber_a, 'eura_a':eura_a, 'iber_c':iber_c, 'eura_c':eura_c,
                'iber_alpha':2*iber_alpha, 'eura_alpha':2*eura_alpha,
                'Tb':Tb/2, 'Tsum':(Tb+Tsplit)/2}
    
    return command % sub_dict

##########################################################################################

def model_2_a_mscore(params):
    """
    ms core command corresponding to model_2_a
    model_2_a - two bottlenecks, both abrupt change:
     # population split with different sizes
     # abrupt change in Iberian lynx population at time Tbot_a
     # exponential change in Eurasian lynx population as function of time since Tbot_a
     # abrupt change in Iberian lynx population at time Tbot
    """
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr = params
    
    # Bottleneck times
    Tb1 = Tbot_a*Tsplit
    Tb2 = Tbot*Tbot_a*Tsplit
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb1
    
    # command
    command = "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
              "-eg 0 2 %(eura_alpha)f "\
              "-en %(Tb2)f 1 %(iber_i)f "\
              "-en %(Tb1)f 1 %(iber_a)f "\
              "-en %(Tb1)f 2 %(eura_a)f "\
              "-ej %(Tsum)f 2 1"
    
    # transform dadi parameters to ms parameters
    sub_dict = {'iber_a':iber_a, 'eura_a':eura_a, 'iber_i':iber_i, 'iber_c':iber_c, 'eura_c':eura_c,
                'eura_alpha':2*eura_alpha,
                'Tb2':Tb2/2, 'Tb1':(Tb2+Tb1)/2, 'Tsum':(Tb2+Tb1+Tsplit)/2}
    
    return command % sub_dict

##########################################################################################

def model_2_b_mscore(params):
    """
    ms core command corresponding to model_2_b
    model_2_b - two bottlenecks, first exponential change then abrupt
     # population split with different sizes 
     # exponential change in Iberian lynx population as function of time between Tbot_a and Tbot
     # exponential change in Eurasian lynx population as function of time since Tbot_a
     # abrupt change in Iberian lynx population at time Tbot
    """
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr = params
    
    # Bottleneck times
    Tb1 = Tbot_a*Tsplit
    Tb2 = Tbot*Tbot_a*Tsplit
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha = np.log(iber_i/iber_a)/Tb1
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb1
    
    
    # command
    command = "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
              "-eg 0 2 %(eura_alpha)f "\
              "-en %(Tb2)f 1 %(iber_i)f "\
              "-eg %(Tb2)f 1 %(iber_alpha)f "\
              "-en %(Tb1)f 1 %(iber_a)f "\
              "-en %(Tb1)f 2 %(eura_a)f "\
              "-ej %(Tsum)f 2 1"
    
    # transform dadi parameters to ms parameters
    sub_dict = {'iber_a':iber_a, 'eura_a':eura_a, 'iber_i':iber_i, 'iber_c':iber_c, 'eura_c':eura_c,
                'iber_alpha':2*iber_alpha, 'eura_alpha':2*eura_alpha,
                'Tb2':Tb2/2, 'Tb1':(Tb2+Tb1)/2, 'Tsum':(Tb2+Tb1+Tsplit)/2}
    
    return command % sub_dict

def model_2_c_mscore(params):
    """
    ms core command corresponding to model_2_c
    model_2_c - two bottlenecks, first abrupt change then exponential
     # population split with different sizes
     # abrupt change in Iberian lynx population at time Tbot_a
     # exponential change in Eurasian lynx population as function of time since Tbot_a
     # exponential change in Iberian lynx population as function of time since Tbot
    """
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr = params
    
    # Bottleneck times
    Tb1 = Tbot_a*Tsplit
    Tb2 = Tbot*Tbot_a*Tsplit
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha = np.log(iber_c/iber_i)/Tb2
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb1
    
    
    # command
    command = "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
              "-eg 0 1 %(iber_alpha)f "\
              "-eg 0 2 %(eura_alpha)f "\
              "-en %(Tb2)f 1 %(iber_i)f "\
              "-en %(Tb1)f 1 %(iber_a)f "\
              "-en %(Tb1)f 2 %(eura_a)f "\
              "-ej %(Tsum)f 2 1"
    
    # transform dadi parameters to ms parameters
    sub_dict = {'iber_a':iber_a, 'eura_a':eura_a, 'iber_i':iber_i, 'iber_c':iber_c, 'eura_c':eura_c,
                'iber_alpha':2*iber_alpha, 'eura_alpha':2*eura_alpha,
                'Tb2':Tb2/2, 'Tb1':(Tb2+Tb1)/2, 'Tsum':(Tb2+Tb1+Tsplit)/2}
    
    return command % sub_dict
    