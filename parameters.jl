display_plots = true #Whether to open Plots during the computatin. If false, then they are still saved as files in any case.

#BLACK HOLE PARAMETERS

m0s = mp*1e10*ones(1) #BH mass in gram
ks = 0.0*ones(1) #Burden exponent
qs = 0.5*ones(1) #Burden q (when memory burden sets in)
bursts = trues(1) #Whether T~const. or T~1/M during memory burden
a_bh = 0.0 #Value of initial BH spin parameter
two_stages = true #Whether to include two stages of the evaporation (does not work with fit!)

#PARTICLE DM PARAMETERS

m_ncdms = 1e-1*ones(1) #DM mass in GeV
s = 0.5 #Spin of DM particle
g = 1 #DOF of DM particle

#COMOLOGICAL PARAMETERS

#Specify one of the following three parameters:
β = nothing #Initial PBH density fraction. If set to nothing, f is used.
f = nothing #Initial PBH dark matter fraciton. If set to nothing then Ω_ncdm is used.
Ω_ncdm = 0.12 #If set to nothing, then it is computed.
md = nothing #If set to nothing, then it is computed whether there is matter domination phase.

Ω_dm = 0.1201075 #Total dark matter density

#PSD PARAMETERS

use_fits = falses(1) #Whether to use a PSD file or an analytical expression

n_e_class = 200 #Number of values for q=E/T for CLASS
q_min_class = 0.01 #Minimum value for q=E/T for CLASS
q_max_class = 100 #Maximum value for q=E/T for CLASS
n_e = 1000 #Number of values for q=E/T for plots and fitting
q_min = 1e-4 #Minimum value for q=E/T for plots and fitting
q_max = 1e4 #Maximum value for q=E/T for plots and fitting