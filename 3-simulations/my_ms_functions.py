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
 	- m21 for migration from Lynx lynx into Lynx pardinus (from pop2 into pop1)
  - m12 for migration from Lynx pardinus into Lynx lynx (from pop1 into pop2)
 	- m0 for no migration
 - Population 2 (Eurasian lynx) name

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
    to be drawn as values between 0.5 times the value and the value

    sorry you had to read that
    """
    low_bound = param * 0.5
    up_bound = param * 1.5
    draw = np.random.uniform(low=low_bound, high=up_bound)
    if draw > 1:
        draw = np.random.uniform(low=param, high=1)
    return draw

##########################################################################################

############### FUNCTIONS GENERATE AN MS COMMAND FROM THE INPUT PARAMETERS ###############

def model_2_d_ms_m21(params):
    """
    Use input parameters from dadi to build an MS command corresponding to model_2_d
    with one migration from pop2 into pop1
    
    model_2_d consists of population divergence followed by:
     - two consecutive exponential bottlenecks in iberian lynx
     - one exponential bottleneck in eurasian lynx

    """
    
    Tsplit, Tbot2, Tbot1, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, theta, pop2 = params
    
    """
    the input parameters are defined as follows from the dadi function:
    
    Tsplit = time span between divergence and first bottleneck event
    
    Tbot1 = time span between first and second bottleneck events
           
    Tbot2 = time span between second bottleneck event and present
    
    iber_a = population size of iberian lynx at time of divergence (ancient)
    
    iber_pr_a = proportion of iberian lynx population left after first bottleneck
                size after first bottleneck = iber_a * iber_pr_a
                between 0 and 1

    iber_pr = proportion of iberian lynx population left after second bottleneck
                size after second bottleneck = iber_a * iber_pr_a * iber_pr
                between 0 and 1
    
    eura_a = population size of eurasian lynx at time of divergence (ancient)
    
    eura_pr = proportion of eurasian lynx population after bottleneck
                size after bottleneck = eura_a * eura_pr
                between 0 and 1
    
    theta = optimal multiplicative scaling factor between model and data (dadi.Inference.optimal_sfs_scaling)
            value for whole "callable" genome
    
    pop2 = name of the eurasian lynx population
    """
    
    # draw parameter values from uniform distribution around optimized value:
    Tsplit = draw_uni(Tsplit)
    Tbot1 = draw_uni(Tbot1)
    Tbot2 = draw_uni(Tbot2)
    iber_a = draw_uni(iber_a)
    iber_pr_a = draw_uni_cap(iber_pr_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
    if pop2 == "ki" or "wel":
        n_pop2 = 40
    if pop2 == "eel":
        n_pop2 = 38
    if pop2 == "sel":
        n_pop2 = 24
    
    """
    input values are expected to be optimized
    
    the function draw_uni will draw a value from a uniform distribution that goes 
    from 0.5 * value to 1.5 * value
    
    the function draw_uni_cap will do the same but it is applied to values that go
    from 0 to 1 so that values above the value are equally as likely to be drawn as
    values below the value
    
    sample size of the eurasian lynx (population 2) is variable and depends on what population is chosen
    sample size of iberian lynx is always the same so we don't need to define it
    """
    
    # theta rescaling to 100KB/EffectiveGenomeLength - see (1-prepare_dataset.md)
    if pop2 == "wel":
        rescale_factor = 100000/612311182
    if pop2 == "eel":
        rescale_factor = 100000/612178519
    if pop2 == "sel":
        rescale_factor = 100000/612085036

    theta = theta * rescale_factor
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha_2 = np.log(iber_c/iber_i)/Tbot2
    iber_alpha_1 = np.log(iber_i/iber_a)/Tbot1
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tbot1
 
    # recombination rho:
    rho = theta * (0.019 / 0.006)
    nsites = 10000
   
    # transform dadi units to ms units:
    Tsum = (Tbot2+Tbot1+Tsplit)/2
    Tb1 = (Tbot2+Tbot1)/2
    Tb2 = Tbot2/2
    eura_alpha = 2 * eura_alpha
    iber_alpha_1 = 2 * iber_alpha_1
    iber_alpha_2 = 2 * iber_alpha_2

    # migration parameters
    Tm = np.random.uniform(0, 0.25 * Tsum)
    Mprob =  np.random.uniform(0,1)
   
    """
    the parameter values need to be transformed and rescaled in order to be used by MS:
    
    theta = rescaled from effective genome size to the size of the window to be simulated
    
    proportional times and sizes are calculated from original values:
        iber_i = size after first bottleneck (intermediate) = iber_a * iber_pr_a
        iber_c = size after second bottleneck (contemporary) = iber_a*iber_pr_a*iber_pr
        eura_c = size after bottleneck (contemporary) = eura_a * eura_pr
        
    growth rates are calculated as the log of the ratio between starting size and ending size (going back in time)
        divided by the time length of the bottleneck event

    rho is calculated from theta considering a recombination rate of 1.9 cM/Mbp and
        a mutation rate of 6*10e-9 per site per generation
    
    nsite is set to 10000 to have significantly more recombination events than segregating sites in the
        simulated data (~200 expected)
    
    units are transformed from dadi to ms scales by dividing times by a factor of 2
    
    time in MS is also the sum of individual time frames calculated by dadi, going from present to the past
    
    migration parameters are also generated by drawing the time of migration from a uniform distribution 
        between 0 and 0.25 * time of divergence, 
        and the probability of migration from a uniform distribution between 0 and 1
    """
    
    # command if time of migration is < than time of second bottleneck:
    if Tm < Tb2:
        command = "-t %(theta)f -I 2 38 %(n_pop2)i "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-ev %(Tm)f 2 1 %(Mprob)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)i"
                  
    # command if time of second bottleneck is <= than time of migration < time of first bottleneck
    if Tb2 <= Tm and Tm < Tb1:
        command = "-t %(theta)f -I 2 38 %(n_pop2)i "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-ev %(Tm)f 2 1 %(Mprob)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)i"

    # command if time of first bottleneck is <= than time of migration
    if Tb1 <= Tm:
        command = "-t %(theta)f -I 2 38 %(n_pop2)i "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ev %(Tm)f 2 1 %(Mprob)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)i"

    """
    The actual MS command!
    
    Timeline going from present to past:
    
    (0) the migration event at time Tm is positioned always before Tsum (by definition)
        and can fall anywhere relatively to the two bottleneck events, hence the if statements

    (1) iberian is population 1 and eurasian is population 2
    
    (2) the two populations have sizes iber_c and eura_c (iberian and eurasian contemporary)
    
    (3) iberian lynx experiences exponential change:
            from population sizes iber_c to iber_i (iberian intermediate)
            from time 0 to time Tb2 (time of second bottleneck)
            calculated as iber_alpha_2 (growth factor)
    
    (4) eurasian lynx experiences exponential change:
            from population sizes eura_c to eura_a (eurasian ancient)
            from time 0 to time Tb1 (time of first bottleneck)
            calculated as eura_alpha (growth factor)

    (5) iberian population size is set to iber_i at time Tb2 to stop exponential change
    
    (6) iberian lynx experiences another exponential change:
            from population sizes iber_i to iber_a (iberian ancient)
            from time Tb2 to time Tb1 (time of first bottleneck)
            calculated as iber_alpha_1 (growth factor)
    
    (7) iberian and eurasian population sizes are set to iber_a and eura_a
        at time TB1 to stop exponential change
    
    (8) all the individuals from one population are moved into the other at
        time Tsum (divergence event backwards)
    
    """

    # substitute values in command:
    sub_dict = {'theta':theta, 'n_pop2':n_pop2,
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
    Use input parameters from dadi to build an MS command corresponding to model_2_d
    with one migration from pop1 into pop2
    
    model_2_d consists of population divergence followed by:
     - two consecutive exponential bottlenecks in iberian lynx
     - one exponential bottleneck in eurasian lynx

    """
    
    Tsplit, Tbot2, Tbot1, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, theta, pop2 = params
    
    """
    the input parameters are defined as follows from the dadi function:
    
    Tsplit = time span between divergence and first bottleneck event
    
    Tbot1 = time span between first and second bottleneck events
           
    Tbot2 = time span between second bottleneck event and present
    
    iber_a = population size of iberian lynx at time of divergence (ancient)
    
    iber_pr_a = proportion of iberian lynx population left after first bottleneck
                size after first bottleneck = iber_a * iber_pr_a
                between 0 and 1

    iber_pr = proportion of iberian lynx population left after second bottleneck
                size after second bottleneck = iber_a * iber_pr_a * iber_pr
                between 0 and 1
    
    eura_a = population size of eurasian lynx at time of divergence (ancient)
    
    eura_pr = proportion of eurasian lynx population after bottleneck
                size after bottleneck = eura_a * eura_pr
                between 0 and 1
    
    theta = optimal multiplicative scaling factor between model and data (dadi.Inference.optimal_sfs_scaling)
            value for whole "callable" genome
    
    pop2 = name of the eurasian lynx population
    """
    
    # draw parameter values from uniform distribution around optimized value:
    Tsplit = draw_uni(Tsplit)
    Tbot1 = draw_uni(Tbot1)
    Tbot2 = draw_uni(Tbot2)
    iber_a = draw_uni(iber_a)
    iber_pr_a = draw_uni_cap(iber_pr_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
    if pop2 == "ki" or "wel":
        n_pop2 = 40
    if pop2 == "eel":
        n_pop2 = 38
    if pop2 == "sel":
        n_pop2 = 24
    
    """
    input values are expected to be optimized
    
    the function draw_uni will draw a value from a uniform distribution that goes 
    from 0.5 * value to 1.5 * value
    
    the function draw_uni_cap will do the same but it is applied to values that go
    from 0 to 1 so that values above the value are equally as likely to be drawn as
    values below the value
    
    sample size of the eurasian lynx (population 2) is variable and depends on what population is chosen
    sample size of iberian lynx is always the same so we don't need to define it
    """
    
    # theta rescaling to 100KB/EffectiveGenomeLength - see (1-prepare_dataset.md)
    if pop2 == "wel":
        rescale_factor = 100000/612311182
    if pop2 == "eel":
        rescale_factor = 100000/612178519
    if pop2 == "sel":
        rescale_factor = 100000/612085036

    theta = theta * rescale_factor
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha_2 = np.log(iber_c/iber_i)/Tbot2
    iber_alpha_1 = np.log(iber_i/iber_a)/Tbot1
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tbot1
 
    # recombination rho:
    rho = theta * (0.019 / 0.006)
    nsites = 10000
   
    # transform dadi units to ms units:
    Tsum = (Tbot2+Tbot1+Tsplit)/2
    Tb1 = (Tbot2+Tbot1)/2
    Tb2 = Tbot2/2
    eura_alpha = 2 * eura_alpha
    iber_alpha_1 = 2 * iber_alpha_1
    iber_alpha_2 = 2 * iber_alpha_2

    # migration parameters
    Tm = np.random.uniform(0, 0.25 * Tsum)
    Mprob =  np.random.uniform(0,1)
   
    """
    the parameter values need to be transformed and rescaled in order to be used by MS:
    
    theta = rescaled from effective genome size to the size of the window to be simulated
    
    proportional times and sizes are calculated from original values:
        iber_i = size after first bottleneck (intermediate) = iber_a * iber_pr_a
        iber_c = size after second bottleneck (contemporary) = iber_a*iber_pr_a*iber_pr
        eura_c = size after bottleneck (contemporary) = eura_a * eura_pr
        
    growth rates are calculated as the log of the ratio between starting size and ending size (going back in time)
        divided by the time length of the bottleneck event

    rho is calculated from theta considering a recombination rate of 1.9 cM/Mbp and
        a mutation rate of 6*10e-9 per site per generation
    
    nsite is set to 10000 to have significantly more recombination events than segregating sites in the
        simulated data (~200 expected)
    
    units are transformed from dadi to ms scales by dividing times by a factor of 2
    
    time in MS is also the sum of individual time frames calculated by dadi, going from present to the past
    
    migration parameters are also generated by drawing the time of migration from a uniform distribution 
        between 0 and 0.25 * time of divergence, 
        and the probability of migration from a uniform distribution between 0 and 1
    """
    
    # command if time of migration is < than time of second bottleneck:
    if Tm < Tb2:
        command = "-t %(theta)f -I 2 38 %(n_pop2)i "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-ev %(Tm)f 1 2 %(Mprob)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)i"
                  
    # command if time of second bottleneck is <= than time of migration < time of first bottleneck
    if Tb2 <= Tm and Tm < Tb1:
        command = "-t %(theta)f -I 2 38 %(n_pop2)i "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-ev %(Tm)f 1 2 %(Mprob)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)i"

    # command if time of first bottleneck is <= than time of migration
    if Tb1 <= Tm:
        command = "-t %(theta)f -I 2 38 %(n_pop2)i "\
                  "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
                  "-eg 0 1 %(iber_alpha_2)f "\
                  "-eg 0 2 %(eura_alpha)f "\
                  "-en %(Tb2)f 1 %(iber_i)f "\
                  "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
                  "-en %(Tb1)f 1 %(iber_a)f "\
                  "-en %(Tb1)f 2 %(eura_a)f "\
                  "-ev %(Tm)f 1 2 %(Mprob)f "\
                  "-ej %(Tsum)f 2 1 "\
                  "-r %(rho)f %(nsites)i"

    """
    The actual MS command!
    
    Timeline going from present to past:
    
    (0) the migration event at time Tm is positioned always before Tsum (by definition)
        and can fall anywhere relatively to the two bottleneck events, hence the if statements

    (1) iberian is population 1 and eurasian is population 2
    
    (2) the two populations have sizes iber_c and eura_c (iberian and eurasian contemporary)
    
    (3) iberian lynx experiences exponential change:
            from population sizes iber_c to iber_i (iberian intermediate)
            from time 0 to time Tb2 (time of second bottleneck)
            calculated as iber_alpha_2 (growth factor)
    
    (4) eurasian lynx experiences exponential change:
            from population sizes eura_c to eura_a (eurasian ancient)
            from time 0 to time Tb1 (time of first bottleneck)
            calculated as eura_alpha (growth factor)

    (5) iberian population size is set to iber_i at time Tb2 to stop exponential change
    
    (6) iberian lynx experiences another exponential change:
            from population sizes iber_i to iber_a (iberian ancient)
            from time Tb2 to time Tb1 (time of first bottleneck)
            calculated as iber_alpha_1 (growth factor)
    
    (7) iberian and eurasian population sizes are set to iber_a and eura_a
        at time TB1 to stop exponential change
    
    (8) all the individuals from one population are moved into the other at
        time Tsum (divergence event backwards)
    
    """

    # substitute values in command:
    sub_dict = {'theta':theta, 'n_pop2':n_pop2,
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

def model_2_d_ms_m0(params):
    """
    Use input parameters from dadi to build an MS command corresponding to model_2_d
    without any migration events
    
    model_2_d consists of population divergence followed by:
     - two consecutive exponential bottlenecks in iberian lynx
     - one exponential bottleneck in eurasian lynx

    """
    
    Tsplit, Tbot2, Tbot1, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, theta, pop2 = params
    
    """
    the input parameters are defined as follows from the dadi function:
    
    Tsplit = time span between divergence and first bottleneck event
    
    Tbot1 = time span between first and second bottleneck events
           
    Tbot2 = time span between second bottleneck event and present
    
    iber_a = population size of iberian lynx at time of divergence (ancient)
    
    iber_pr_a = proportion of iberian lynx population left after first bottleneck
                size after first bottleneck = iber_a * iber_pr_a
                between 0 and 1

    iber_pr = proportion of iberian lynx population left after second bottleneck
                size after second bottleneck = iber_a * iber_pr_a * iber_pr
                between 0 and 1
    
    eura_a = population size of eurasian lynx at time of divergence (ancient)
    
    eura_pr = proportion of eurasian lynx population after bottleneck
                size after bottleneck = eura_a * eura_pr
                between 0 and 1
    
    theta = optimal multiplicative scaling factor between model and data (dadi.Inference.optimal_sfs_scaling)
            value for whole "callable" genome
    
    pop2 = name of the eurasian lynx population
    """
    
    # draw parameter values from uniform distribution around optimized value:
    Tsplit = draw_uni(Tsplit)
    Tbot1 = draw_uni(Tbot1)
    Tbot2 = draw_uni(Tbot2)
    iber_a = draw_uni(iber_a)
    iber_pr_a = draw_uni_cap(iber_pr_a)
    iber_pr = draw_uni_cap(iber_pr)
    eura_a = draw_uni(eura_a)
    eura_pr = draw_uni_cap(eura_pr)
    theta = draw_uni(theta)
    if pop2 == "ki" or "wel":
        n_pop2 = 40
    if pop2 == "eel":
        n_pop2 = 38
    if pop2 == "sel":
        n_pop2 = 24
    
    """
    input values are expected to be optimized
    
    the function draw_uni will draw a value from a uniform distribution that goes 
    from 0.5 * value to 1.5 * value
    
    the function draw_uni_cap will do the same but it is applied to values that go
    from 0 to 1 so that values above the value are equally as likely to be drawn as
    values below the value
    
    sample size of the eurasian lynx (population 2) is variable and depends on what population is chosen
    sample size of iberian lynx is always the same so we don't need to define it
    """
    
    # theta rescaling to 100KB/EffectiveGenomeLength - see (1-prepare_dataset.md)
    if pop2 == "wel":
        rescale_factor = 100000/612311182
    if pop2 == "eel":
        rescale_factor = 100000/612178519
    if pop2 == "sel":
        rescale_factor = 100000/612085036

    theta = theta * rescale_factor
    
    # intermediate population:
    iber_i = iber_a*iber_pr_a
    # contemporary populations:
    iber_c = iber_a*iber_pr_a*iber_pr
    eura_c = eura_a*eura_pr
    
    # growth rate in iberian:
    iber_alpha_2 = np.log(iber_c/iber_i)/Tbot2
    iber_alpha_1 = np.log(iber_i/iber_a)/Tbot1
    
    # growth rate in eurasian:
    eura_alpha = np.log(eura_c/eura_a)/Tbot1
 
    # recombination rho:
    rho = theta * (0.019 / 0.006)
    nsites = 10000
   
    # transform dadi units to ms units:
    Tsum = (Tbot2+Tbot1+Tsplit)/2
    Tb1 = (Tbot2+Tbot1)/2
    Tb2 = Tbot2/2
    eura_alpha = 2 * eura_alpha
    iber_alpha_1 = 2 * iber_alpha_1
    iber_alpha_2 = 2 * iber_alpha_2

    """
    the parameter values need to be transformed and rescaled in order to be used by MS:
    
    theta = rescaled from effective genome size to the size of the window to be simulated
    
    proportional times and sizes are calculated from original values:
        iber_i = size after first bottleneck (intermediate) = iber_a * iber_pr_a
        iber_c = size after second bottleneck (contemporary) = iber_a*iber_pr_a*iber_pr
        eura_c = size after bottleneck (contemporary) = eura_a * eura_pr
        
    growth rates are calculated as the log of the ratio between starting size and ending size (going back in time)
        divided by the time length of the bottleneck event

    rho is calculated from theta considering a recombination rate of 1.9 cM/Mbp and
        a mutation rate of 6*10e-9 per site per generation
    
    nsite is set to 10000 to have significantly more recombination events than segregating sites in the
        simulated data (~200 expected)
    
    units are transformed from dadi to ms scales by dividing times by a factor of 2
    
    time in MS is also the sum of individual time frames calculated by dadi, going from present to the past
    """
    
    # command:
    command = "-t %(theta)f -I 2 38 %(n_pop2)i "\
              "-n 1 %(iber_c)f -n 2 %(eura_c)f "\
              "-eg 0 1 %(iber_alpha_2)f "\
              "-eg 0 2 %(eura_alpha)f "\
              "-en %(Tb2)f 1 %(iber_i)f "\
              "-eg %(Tb2)f 1 %(iber_alpha_1)f "\
              "-en %(Tb1)f 1 %(iber_a)f "\
              "-en %(Tb1)f 2 %(eura_a)f "\
              "-ej %(Tsum)f 2 1 "\
              "-r %(rho)f %(nsites)i"

    """
    The actual MS command!
    
    Timeline going from present to past:
    
    (1) iberian is population 1 and eurasian is population 2
    
    (2) the two populations have sizes iber_c and eura_c (iberian and eurasian contemporary)
    
    (3) iberian lynx experiences exponential change:
            from population sizes iber_c to iber_i (iberian intermediate)
            from time 0 to time Tb2 (time of second bottleneck)
            calculated as iber_alpha_2 (growth factor)
    
    (4) eurasian lynx experiences exponential change:
            from population sizes eura_c to eura_a (eurasian ancient)
            from time 0 to time Tb1 (time of first bottleneck)
            calculated as eura_alpha (growth factor)

    (5) iberian population size is set to iber_i at time Tb2 to stop exponential change
    
    (6) iberian lynx experiences another exponential change:
            from population sizes iber_i to iber_a (iberian ancient)
            from time Tb2 to time Tb1 (time of first bottleneck)
            calculated as iber_alpha_1 (growth factor)
    
    (7) iberian and eurasian population sizes are set to iber_a and eura_a
        at time TB1 to stop exponential change
    
    (8) all the individuals from one population are moved into the other at
        time Tsum (divergence event backwards)
    
    """

    # substitute values in command:
    sub_dict = {'theta':theta, 'n_pop2':n_pop2,
                'iber_c':iber_c, 'eura_c':eura_c,
                'iber_alpha_2':iber_alpha_2,
                'iber_alpha_1':iber_alpha_1,
                'eura_alpha':eura_alpha,
                'iber_i':iber_i, 'iber_a':iber_a, 'eura_a':eura_a,
                'Tb2':Tb2, 'Tb1':Tb1, 'Tsum':Tsum,
                'rho':rho, 'nsites':nsites}
    
    return command % sub_dict

##########################################################################################

