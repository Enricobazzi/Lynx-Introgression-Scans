import numpy as np
import dadi
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

"""
Each of these functions corresponds to one of the models defined by my_models.py, and it will
build a ms function for that model, indicating the demography to be simulated by ms.

The following information is required for the function to be generated:

 - Optimized parameters of the model
 - Theta of the optimized model
 - Description of the migration event to be included:
 	- m0 for no migration
 	- m21 for migration from Lynx lynx into Lynx pardinus (from pop2 into pop1)
  	- m12 for migration from Lynx pardinus into Lynx lynx (from pop1 into pop2)

Values will be drawn from uniform distributions of:
 
 - parameters = from 0.5 * optimized parameter to 1.5 * optimized parameter
 - theta = from 0.5 * optimized theta to 1.5 * optimized theta
 - time of migration = from 0 to 0.25 * time of divergence
 - probability of migrantion = from 0 to 1

"""
##########################################################################################

#### FUNCTIONS TO DRAW RANDOM NUMBER FROM UNIFORM DISTRIBUTION AROUND PARAMETER VALUE ####

def draw_uni(param):
    """
    function to draw a random number from a uniform distribution between 0.5 and 1.5 times the value
    """
    low_bound = param * 0.5
    up_bound = param * 1.5
    draw = np.random.uniform(low=low_bound, high=up_bound)
    return draw

#############################################

def draw_uni_cap(param):
    """
    function to draw a random number from a uniform distribution between 0.5 and 1.5 times the value

    including an upper limit to the drawn value of 1:

    if a values > 1 are present in the uniform distribution, values between the value and 1 will be as likely
    to be drawn as values between 0.5 and the value

    sorry you had to read that
    """
    low_bound = param * 0.5
    up_bound = param * 1.5
    draw = np.random.uniform(low=low_bound, high=up_bound)
    if draw > 1:
        draw = 1 - (draw - 1)
    return draw

##########################################################################################

############### FUNCTIONS GENERATE AN MS COMMAND FROM THE INPUT PARAMETERS ###############

def model_1_a_ms(params):
    """
    ms command corresponding to model_1_a
    model_1_a - one bottleneck, abrupt change
	# population split with different sizes
	# abrupt change in Iberian lynx population at time Tbot
	# exponential change in Eurasian lynx population as function of time since Tbot
    """
    Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr, theta = params
    
    # draw parameter values from uniform distribution around optimized value
    Tsplit = draw_uni(Tsplit)
    Tbot = draw_uni_cap(Tbot)
    iber_a = draw_uni(iber_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
        
    # theta rescaling to 10KB/GenomeLength
    rescale_factor = 10000/2400000000
    theta = theta * rescale_factor
    
	# bottleneck time:
    Tb = Tbot*Tsplit
        
	# contemporary populations:
    iber_c = iber_a*iber_pr
    eura_c = eura_a*eura_pr
    
	# growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb
    
    # recombination rho:
    rho = theta * (0.019 / 0.0006)
    nsites = theta * 10
    
    # transform dadi units to ms units
    Tsum = (Tb+Tsplit)/2
    Tb = Tb/2
    eura_alpha = 2 * eura_alpha
        
	# command
    command = "-t %(theta)f -I 2 1 1 "\
              "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
              "-eg 0 2 %(eura_alpha)f "\
              "-en %(Tb)f 1 %(iber_a)f "\
              "-en %(Tb)f 2 %(eura_a)f "\
              "-ej %(Tsum)f 2 1 "\
              "-r %(rho)f %(nsites)f"

	# transform dadi parameters to ms parameters
    sub_dict = {'theta':theta,
                'iber_c':iber_c, 'eura_c':eura_c,
                'eura_alpha':eura_alpha,
                'iber_a':iber_a, 'eura_a':eura_a,
                'Tb':Tb, 'Tsum':Tsum,
                'rho':rho, 'nsites':nsites}
    
    return command % sub_dict

#############################################

def model_1_a_ms_m21(params):
    """
    ms command corresponding to model_1_a with migration from pop2 into pop1
    model_1_a - one bottleneck, abrupt change
	# population split with different sizes
	# abrupt change in Iberian lynx population at time Tbot
	# exponential change in Eurasian lynx population as function of time since Tbot
	# migration event from from Eurasian lynx into Iberian Lynx (from pop2 into pop1)
    """
    Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr, theta = params
    
    # draw parameter values from uniform distribution around optimized value
    Tsplit = draw_uni(Tsplit)
    Tbot = draw_uni_cap(Tbot)
    iber_a = draw_uni(iber_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
        
    # theta rescaling to 10KB/GenomeLength
    rescale_factor = 10000/2400000000
    theta = theta * rescale_factor
    
	# bottleneck time:
    Tb = Tbot*Tsplit
        
	# contemporary populations:
    iber_c = iber_a*iber_pr
    eura_c = eura_a*eura_pr
    
	# growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb
    
    # recombination rho:
    rho = theta * (0.019 / 0.0006)
    nsites = theta * 10
    
    # transform dadi units to ms units:
    Tsum = (Tb+Tsplit)/2
    Tb = Tb/2
    eura_alpha = 2 * eura_alpha
    
    # migration parameters
    Tm = np.random.uniform(0, 0.25 * Tsum)
    Mprob =  np.random.uniform(0,1)
    
    # command if time of migration is < than time of bottleneck:
    if Tm < Tb:
        command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-ev %(Tm)f 2 1 %(Mprob)f "\
                  "-en %(Tb)f 1 %(iber_a)f "\
                  "-en %(Tb)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"

    # command if time of bottleneck <= than time of migration:
    if Tb <= Tm:
	    command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb)f 1 %(iber_a)f "\
                  "-en %(Tb)f 2 %(eura_a)f "\
                  "-ev %(Tm)f 2 1 %(Mprob)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"
                  
	# transform dadi parameters to ms parameters
    sub_dict = {'theta':theta,
                'iber_c':iber_c, 'eura_c':eura_c,
                'eura_alpha':eura_alpha,
                'Tm':Tm, 'Mprob':Mprob,
                'iber_a':iber_a, 'eura_a':eura_a,
                'Tb':Tb, 'Tsum':Tsum,
                'rho':rho, 'nsites':nsites}
    
    return command % sub_dict
    
#############################################

def model_1_a_ms_m12(params):
    """
    ms command corresponding to model_1_a with migration from pop1 into pop2
    model_1_a - one bottleneck, abrupt change
	# population split with different sizes
	# abrupt change in Iberian lynx population at time Tbot
	# exponential change in Eurasian lynx population as function of time since Tbot
	# migration event from from Iberian lynx into Eurasian Lynx (from pop1 into pop2)
    """
    Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr, theta = params
    
    # draw parameter values from uniform distribution around optimized value
    Tsplit = draw_uni(Tsplit)
    Tbot = draw_uni_cap(Tbot)
    iber_a = draw_uni(iber_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
        
    # theta rescaling to 10KB/GenomeLength
    rescale_factor = 10000/2400000000
    theta = theta * rescale_factor
    
	# bottleneck time:
    Tb = Tbot*Tsplit
        
	# contemporary populations:
    iber_c = iber_a*iber_pr
    eura_c = eura_a*eura_pr
    
	# growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb
    
    # recombination rho:
    rho = theta * (0.019 / 0.0006)
    nsites = theta * 10
    
    # transform dadi units to ms units:
    Tsum = (Tb+Tsplit)/2
    Tb = Tb/2
    eura_alpha = 2 * eura_alpha
    
    # migration parameters
    Tm = np.random.uniform(0, 0.25 * Tsum)
    Mprob =  np.random.uniform(0,1)
    
    # command if time of migration is < than time of bottleneck:
    if Tm < Tb:
        command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-ev %(Tm)f 1 2 %(Mprob)f "\
                  "-en %(Tb)f 1 %(iber_a)f "\
                  "-en %(Tb)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"

    # command if time of bottleneck <= than time of migration:
    if Tb <= Tm:
	    command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb)f 1 %(iber_a)f "\
                  "-en %(Tb)f 2 %(eura_a)f "\
                  "-ev %(Tm)f 1 2 %(Mprob)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"
                  
	# transform dadi parameters to ms parameters
    sub_dict = {'theta':theta,
                'iber_c':iber_c, 'eura_c':eura_c,
                'eura_alpha':eura_alpha,
                'Tm':Tm, 'Mprob':Mprob,
                'iber_a':iber_a, 'eura_a':eura_a,
                'Tb':Tb, 'Tsum':Tsum,
                'rho':rho, 'nsites':nsites}
    
    return command % sub_dict

#############################################

def model_2_d_ms(params):
    """
    ms command corresponding to model_2_d
    model_2_d - two bottlenecks both exponential
     # population split with different sizes
     # exponential change in Iberian lynx population as function of time between Tbot_a and Tbot
     # exponential change in Eurasian lynx population as function of time since Tbot_a
     # exponential change in Iberian lynx population as function of time since Tbot
    """
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, theta = params

    # draw parameter values from uniform distribution around optimized value
    Tsplit = draw_uni(Tsplit)
    Tbot = draw_uni_cap(Tbot)
    Tbot_a = draw_uni_cap(Tbot_a)
    iber_a = draw_uni(iber_a)
    iber_pr_a = draw_uni_cap(iber_pr_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
        
    # theta rescaling to 10KB/GenomeLength
    rescale_factor = 10000/2400000000
    theta = theta * rescale_factor
    
    # Bottleneck times
    Tb1 = Tbot_a*Tsplit
    Tb2 = Tbot*Tbot_a*Tsplit
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha_2 = np.log(iber_c/iber_i)/Tb2
    iber_alpha_1 = np.log(iber_i/iber_a)/Tb1
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb1
 
    # recombination rho:
    rho = theta * (0.019 / 0.0006)
    nsites = theta * 10
   
    # transform dadi units to ms units:
    Tsum = (Tb1+Tb2+Tsplit)/2
    Tb1 = (Tb2+Tb1)/2
    Tb2 = Tb2/2
    eura_alpha = 2 * eura_alpha
    iber_alpha_1 = 2 * iber_alpha_1
    iber_alpha_2 = 2 * iber_alpha_2

	# command
    command = "-t %(theta)f -I 2 1 1 "\
              "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
              "-eg 0 1 %(iber_alpha_2)f "\
              "-eg 0 2 %(eura_alpha)f "\
              "-en %(Tb2)f 1 %(iber_i)f "\
              "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
              "-en %(Tb1)f 1 %(iber_a)f "\
              "-en %(Tb1)f 2 %(eura_a)f "\
              "-ej %(Tsum)f 2 1 "\
              "-r %(rho)f %(nsites)f"

	# transform dadi parameters to ms parameters
    sub_dict = {'theta':theta,
                'iber_c':iber_c, 'eura_c':eura_c,
                'iber_alpha_2':iber_alpha_2,
                'iber_alpha_1':iber_alpha_1,
                'eura_alpha':eura_alpha,
                'iber_i':iber_i, 'iber_a':iber_a, 'eura_a':eura_a,
                'Tb2':Tb2, 'Tb1':Tb1, 'Tsum':Tsum,
                'rho':rho, 'nsites':nsites}
    
    return command % sub_dict

#############################################

def model_2_d_ms_m21(params):
    """
    ms command corresponding to model_2_d with migration from pop2 into pop1
    model_2_d - two bottlenecks both exponential
     # population split with different sizes
     # exponential change in Iberian lynx population as function of time between Tbot_a and Tbot
     # exponential change in Eurasian lynx population as function of time since Tbot_a
     # exponential change in Iberian lynx population as function of time since Tbot
     # migration event from from Eurasian lynx into Iberian Lynx (from pop2 into pop1)
    """
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, theta = params

    # draw parameter values from uniform distribution around optimized value
    Tsplit = draw_uni(Tsplit)
    Tbot = draw_uni_cap(Tbot)
    Tbot_a = draw_uni_cap(Tbot_a)
    iber_a = draw_uni(iber_a)
    iber_pr_a = draw_uni_cap(iber_pr_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
        
    # theta rescaling to 10KB/GenomeLength
    rescale_factor = 10000/2400000000
    theta = theta * rescale_factor
    
    # Bottleneck times
    Tb1 = Tbot_a*Tsplit
    Tb2 = Tbot*Tbot_a*Tsplit
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha_2 = np.log(iber_c/iber_i)/Tb2
    iber_alpha_1 = np.log(iber_i/iber_a)/Tb1
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb1
 
    # recombination rho:
    rho = theta * (0.019 / 0.0006)
    nsites = theta * 10
   
    # transform dadi units to ms units:
    Tsum = (Tb1+Tb2+Tsplit)/2
    Tb1 = (Tb2+Tb1)/2
    Tb2 = Tb2/2
    eura_alpha = 2 * eura_alpha
    iber_alpha_1 = 2 * iber_alpha_1
    iber_alpha_2 = 2 * iber_alpha_2

    # migration parameters
    Tm = np.random.uniform(0, 0.25 * Tsum)
    Mprob =  np.random.uniform(0,1)

    # command if time of migration is < than time of second bottleneck:
    if Tm < Tb2:
        command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-ev %(Tm)f 2 1 %(Mprob)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"
                  
    # command if time of second bottleneck is <= than time of migration < time of first bottleneck
    if Tb2 <= Tm and Tm < Tb1:
        command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-ev %(Tm)f 2 1 %(Mprob)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"

    # command if time of first bottleneck is <= than time of migration
    if Tb1 <= Tm:
        command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ev %(Tm)f 2 1 %(Mprob)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"

	# transform dadi parameters to ms parameters
    sub_dict = {'theta':theta,
                'iber_c':iber_c, 'eura_c':eura_c,
                'Tm':Tm, 'Mprob':Mprob,
                'iber_alpha_2':iber_alpha_2,
                'iber_alpha_1':iber_alpha_1,
                'eura_alpha':eura_alpha,
                'iber_i':iber_i, 'iber_a':iber_a, 'eura_a':eura_a,
                'Tb2':Tb2, 'Tb1':Tb1, 'Tsum':Tsum,
                'rho':rho, 'nsites':nsites}
    
    return command % sub_dict

#############################################

def model_2_d_ms_m12(params):
    """
    ms command corresponding to model_2_d with migration from pop1 into pop2
    model_2_d - two bottlenecks both exponential
     # population split with different sizes
     # exponential change in Iberian lynx population as function of time between Tbot_a and Tbot
     # exponential change in Eurasian lynx population as function of time since Tbot_a
     # exponential change in Iberian lynx population as function of time since Tbot
	# migration event from from Iberian lynx into Eurasian Lynx (from pop1 into pop2)
    """
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, theta = params

    # draw parameter values from uniform distribution around optimized value
    Tsplit = draw_uni(Tsplit)
    Tbot = draw_uni_cap(Tbot)
    Tbot_a = draw_uni_cap(Tbot_a)
    iber_a = draw_uni(iber_a)
    iber_pr_a = draw_uni_cap(iber_pr_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
        
    # theta rescaling to 10KB/GenomeLength
    rescale_factor = 10000/2400000000
    theta = theta * rescale_factor
    
    # Bottleneck times
    Tb1 = Tbot_a*Tsplit
    Tb2 = Tbot*Tbot_a*Tsplit
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha_2 = np.log(iber_c/iber_i)/Tb2
    iber_alpha_1 = np.log(iber_i/iber_a)/Tb1
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tb1
 
    # recombination rho:
    rho = theta * (0.019 / 0.0006)
    nsites = theta * 10
   
    # transform dadi units to ms units:
    Tsum = (Tb1+Tb2+Tsplit)/2
    Tb1 = (Tb2+Tb1)/2
    Tb2 = Tb2/2
    eura_alpha = 2 * eura_alpha
    iber_alpha_1 = 2 * iber_alpha_1
    iber_alpha_2 = 2 * iber_alpha_2

    # migration parameters
    Tm = np.random.uniform(0, 0.25 * Tsum)
    Mprob =  np.random.uniform(0,1)

    # command if time of migration is < than time of second bottleneck:
    if Tm < Tb2:
        command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-ev %(Tm)f 1 2 %(Mprob)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"
                  
    # command if time of second bottleneck is <= than time of migration < time of first bottleneck
    if Tb2 <= Tm and Tm < Tb1:
        command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-ev %(Tm)f 1 2 %(Mprob)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"

    # command if time of first bottleneck is <= than time of migration
    if Tb1 <= Tm:
        command = "-t %(theta)f -I 2 1 1 "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ev %(Tm)f 1 2 %(Mprob)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)f"

	# transform dadi parameters to ms parameters
    sub_dict = {'theta':theta,
                'iber_c':iber_c, 'eura_c':eura_c,
                'Tm':Tm, 'Mprob':Mprob,
                'iber_alpha_2':iber_alpha_2,
                'iber_alpha_1':iber_alpha_1,
                'eura_alpha':eura_alpha,
                'iber_i':iber_i, 'iber_a':iber_a, 'eura_a':eura_a,
                'Tb2':Tb2, 'Tb1':Tb1, 'Tsum':Tsum,
                'rho':rho, 'nsites':nsites}
    
    return command % sub_dict

#############################################

              
##########################################################################################

