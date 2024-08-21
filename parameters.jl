#BLACK HOLE PARAMETERS

m0s = [1e1,1e1,1e1,1e1] #BH mass in gram
ks = [0.5,0.5,0.5,0.5] #Burden exponent
qs = [0.5,0.5,0.5,0.5] #Burden q (when memory burden sets in)
bursts = [false,false,false,false] #Whether T~const. or T~1/M during memory burden
a_bh = 0.0 #Value of initial BH spin parameter
two_stages = true #Whether to include two stages of the evaporation (does not work with fit!)

#PARTICLE DM PARAMETERS

m_ncdms = [0.01,0.1,1,10] #DM mass in GeV
s = 0.5 #Spin of DM particle
g = 1 #DOF of DM particle

#COMOLOGICAL PARAMETERS

#Specify one of the following three parameters:
β = nothing #Initial PBH density fraction. If set to nothing, f is used.
f = nothing #Initial PBH dark matter fraciton. If set to nothing then Ω_ncdm is used.
Ω_ncdm = 0.12 #If set to nothing, then it is computed.

Ω_dm = 0.1201075 #Total dark matter density

#PSD PARAMETERS

use_fits = [false,false,false,false] #Whether to use a PSD file or an analytical expression

n_e_class = 200 #Number of values for q=E/T for CLASS
q_min_class = 0.01 #Minimum value for q=E/T for CLASS
q_max_class = 100 #Maximum value for q=E/T for CLASS
n_e = 1000 #Number of values for q=E/T for CLASS
q_min = 1e-4 #Minimum value for q=E/T for plots and fitting
q_max = 1e3 #Maximum value for q=E/T for plots and fitting