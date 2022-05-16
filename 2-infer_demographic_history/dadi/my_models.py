import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

"""
Define models

Times of bottlenecks are expressed as proportions of time of divergence
Effective population sizes are expressed as proportions of ancestral population size
"""

"""
model 1_a - one bottleneck, abrupt change
 # population split with different sizes and symmetric migration rates
 # abrupt change in Iberian lynx population at time Tbot
 # exponential change in Eurasian lynx population as function of time since Tbot
 # asymmetric migration rates since Tbot
"""
def model_1_a(params, ns, pts):
    Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr, m, m_12, m_21 = params
    xx = Numerics.default_grid(pts)
 # ancestral population
    phi = PhiManip.phi_1D(xx)
 # population split
    phi = PhiManip.phi_1D_to_2D(xx, phi)
 # integrate split with different sizes and symmetric migration rates
    phi = Integration.two_pops(phi, xx, T=Tsplit, nu1=iber_a, nu2=eura_a, m12=m, m21=m)
 # function defining exponential change in Eurasian lynx as function of time since Tbot
    eurac_func = lambda t: eura_a*((eura_pr*eura_a)/eura_a)**(t/(Tbot*Tsplit))
 # integrate the abrupt population change in iberian at time Tbot and asymmetric migration
    phi = Integration.two_pops(phi, xx, T=Tbot*Tsplit, nu1=iber_pr*iber_a, nu2=eurac_func, m12=m_12, m21=m_21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

"""
model 1_b - one bottleneck, exponential change
 # population split with different sizes and symmetric migration rates
 # exponential change in Iberian lynx population as function of time since Tbot
 # exponential change in Eurasian lynx population as function of time since Tbot
 # asymmetric migration rates since Tbot
"""
def model_1_b(params, ns, pts):
    Tsplit, Tbot, iber_a, iber_pr, eura_a, eura_pr, m, m_12, m_21 = params
    xx = Numerics.default_grid(pts)

# ancestral population
    phi = PhiManip.phi_1D(xx)
# population split
    phi = PhiManip.phi_1D_to_2D(xx, phi)
# integrate split with different sizes and symmetric migration rates
    phi = Integration.two_pops(phi, xx, T=Tsplit, nu1=iber_a, nu2=eura_a, m12=m, m21=m)
# function defining exponential change in Iberian lynx as function of time since Tbot
    iberc_func = lambda t: iber_a*((iber_pr*iber_a)/iber_a)**(t/(Tbot*Tsplit))
# function defining exponential change in Eurasian lynx as function of time since Tbot
    eurac_func = lambda t: eura_a*((eura_pr*eura_a)/eura_a)**(t/(Tbot*Tsplit))
# integrate the exponential change and asymmetric migration
    phi = Integration.two_pops(phi, xx, T=Tbot*Tsplit, nu1=iberc_func, nu2=eurac_func, m12=m_12, m21=m_21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

"""
model 2_a - two bottlenecks, both abrupt change
 # population split with different sizes and symmetric migration rates
 # abrupt change in Iberian lynx population at time Tbot_a
 # exponential change in Eurasian lynx population as function of time since Tbot_a
 # asymmetric migration rates between Tbot_a and Tbot
 # abrupt change in Iberian lynx population at time Tbot
 # asymmetric migration rates since Tbot
"""
def model_2_a(params, ns, pts):
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21 = params
    xx = Numerics.default_grid(pts)

 # ancestral population
    phi = PhiManip.phi_1D(xx)
 # population split
    phi = PhiManip.phi_1D_to_2D(xx, phi)
 # integrate split with different sizes and symmetric migration rates
    phi = Integration.two_pops(phi, xx, T=Tsplit, nu1=iber_a, nu2=eura_a, m12=m, m21=m)
 # function defining exponential change in Eurasian lynx as function of time since Tbot_a
    eurac_func = lambda t: eura_a*((eura_pr*eura_a)/eura_a)**(t/(Tbot_a*Tsplit))
 # integrate first abrupt population change (iber_pr_a) in iberian at time Tbot_a, and 
 # the exponential change in Eurasian lynx as function of time since Tbot_a, and 
 # asymmetric migration rates ma
    phi = Integration.two_pops(phi, xx, T=Tbot_a*Tsplit, nu1=iber_pr_a*iber_a, nu2=eurac_func, m12=ma_12, m21=ma_21)
 # integrate second abrupt population change (iber_pr) in iberian at time Tbot and asymmetric migration rates m
    phi = Integration.two_pops(phi, xx, T=Tbot*Tbot_a*Tsplit, nu1=iber_pr*iber_pr_a*iber_a, nu2=eurac_func, m12=m_12, m21=m_21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

"""
model 2_b - two bottlenecks, first exponential change then abrupt
 # population split with different sizes and symmetric migration rates
 # exponential change in Iberian lynx population as function of time between Tbot_a and Tbot
 # exponential change in Eurasian lynx population as function of time since Tbot_a
 # asymmetric migration rates between Tbot_a and Tbot
 # abrupt change in Iberian lynx population at time Tbot
 # asymmetric migration rates since Tbot
"""
def model_2_b(params, ns, pts):
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21 = params
    xx = Numerics.default_grid(pts)

 # ancestral population
    phi = PhiManip.phi_1D(xx)
 # population split
    phi = PhiManip.phi_1D_to_2D(xx, phi)
 # integrate split with different sizes and symmetric migration rates
    phi = Integration.two_pops(phi, xx, T=Tsplit, nu1=iber_a, nu2=eura_a, m12=m, m21=m)
 # function defining exponential change in Eurasian lynx as function of time since Tbot_a
    eurac_func = lambda t: eura_a*((eura_pr*eura_a)/eura_a)**(t/(Tbot_a*Tsplit))
 # function defining exponential change as function of time between Tbot_a and Tbot
    iberc_func_a = lambda t: iber_a*((iber_pr_a*iber_a)/iber_a)**(t/(Tbot_a*Tsplit))
 # integrate first exponential population change (iber_func_a) in iberian since Tbot_a, and 
 # the exponential change in Eurasian lynx as function of time since Tbot_a, and 
 # asymmetric migration rates ma
    phi = Integration.two_pops(phi, xx, T=Tbot_a*Tsplit, nu1=iberc_func_a, nu2=eurac_func, m12=ma_12, m21=ma_21)
 # integrate second abrupt population change (iber_pr) in iberian at time Tbot and asymmetric migration rates m
    phi = Integration.two_pops(phi, xx, T=Tbot*Tbot_a*Tsplit, nu1=iber_pr*iber_pr_a*iber_a, nu2=eurac_func, m12=m_12, m21=m_21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

"""
model 2_c - two bottlenecks, first abrupt change then exponential
 # population split with different sizes and symmetric migration rates
 # abrupt change in Iberian lynx population at time Tbot_a
 # exponential change in Eurasian lynx population as function of time since Tbot_a
 # asymmetric migration rates between Tbot_a and Tbot
 # exponential change in Iberian lynx population as function of time since Tbot
 # asymmetric migration rates since Tbot
"""
def model_2_c(params, ns, pts):
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21 = params
    xx = Numerics.default_grid(pts)

 # ancestral population
    phi = PhiManip.phi_1D(xx)
 # population split
    phi = PhiManip.phi_1D_to_2D(xx, phi)
 # integrate split with different sizes and symmetric migration rates
    phi = Integration.two_pops(phi, xx, T=Tsplit, nu1=iber_a, nu2=eura_a, m12=m, m21=m)
 # function defining exponential change in Eurasian lynx as function of time since Tbot_a
    eurac_func = lambda t: eura_a*((eura_pr*eura_a)/eura_a)**(t/(Tbot_a*Tsplit))
 # integrate first abrupt population change (iber_pr_a) in iberian at time Tbot_a and asymmetric migration rates ma
    phi = Integration.two_pops(phi, xx, T=Tbot_a*Tsplit, nu1=iber_pr_a*iber_a, nu2=eurac_func, m12=ma_12, m21=ma_21)
 # function defining exponential change as function of time between Tbot and present
    iberc_func = lambda t: iber_pr_a*iber_a*((iber_pr*iber_pr_a*iber_a)/iber_pr_a*iber_a)**(t/(Tbot*Tbot_a*Tsplit))
 # integrate second exponential population change (iber_func) in iberian at time Tbot and asymmetric migration rates m
    phi = Integration.two_pops(phi, xx, T=Tbot*Tbot_a*Tsplit, nu1=iberc_func, nu2=eurac_func, m12=ma_12, m21=ma_21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

"""
model 2_d - two bottlenecks, both exponential change
 # population split with different sizes and symmetric migration rates
 # exponential change in Iberian lynx population as function of time between Tbot_a and Tbot
 # exponential change in Eurasian lynx population as function of time since Tbot_a
 # asymmetric migration rates between Tbot_a and Tbot
 # exponential change in Iberian lynx population as function of time since Tbot
 # asymmetric migration rates since Tbot
"""
def model_2_d(params, ns, pts):
    Tsplit, Tbot, Tbot_a, iber_a, iber_pr_a, iber_pr, eura_a, eura_pr, m, ma_12, ma_21, m_12, m_21 = params
    xx = Numerics.default_grid(pts)

 # ancestral population
    phi = PhiManip.phi_1D(xx)
 # population split
    phi = PhiManip.phi_1D_to_2D(xx, phi)
 # integrate split with different sizes and symmetric migration rates
    phi = Integration.two_pops(phi, xx, T=Tsplit, nu1=iber_a, nu2=eura_a, m12=m, m21=m)
 # function defining exponential change as function of time between Tbot_a and Tbot
    iberc_func_a = lambda t: iber_a*((iber_pr_a*iber_a)/iber_a)**(t/(Tbot_a*Tsplit))
 # function defining exponential change in Eurasian lynx as function of time since Tbot_a
    eurac_func = lambda t: eura_a*((eura_pr*eura_a)/eura_a)**(t/Tbot_a*Tsplit)
 # integrate first exponential population change (iber_func_a) in iberian at time Tbot_a and asymmetric migration rates ma
    phi = Integration.two_pops(phi, xx, T=Tbot_a*Tsplit, nu1=iberc_func_a, nu2=eurac_func, m12=ma_12, m21=ma_21)
 # function defining exponential change as function of time between Tbot and present
    iberc_func = lambda t: iber_pr_a*iber_a*((iber_pr*iber_pr_a*iber_a)/iber_pr_a*iber_a)**(t/(Tbot*Tbot_a*Tsplit))
 # integrate second exponential population change (iber_func) in iberian at time Tbot and asymmetric migration rates m
    phi = Integration.two_pops(phi, xx, T=Tbot*Tbot_a*Tsplit, nu1=iberc_func, nu2=eurac_func, m12=ma_12, m21=ma_21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


