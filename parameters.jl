display_plots = true #Whether to open Plots during the computatin. If false, then they are still saved as files in any case.

#BLACK HOLE PARAMETERS

m0s = 1e1*ones(20) #BH mass in gram
ks = LinRange(0.01,0.5,20) #Burden exponent
qs = 0.5*ones(20) #Burden q (when memory burden sets in)
bursts = falses(20) #Whether T~const. or T~1/M during memory burden
a_bh = 0.0 #Value of initial BH spin parameter
two_stages = true #Whether to include two stages of the evaporation (does not work with fit!)

#PARTICLE DM PARAMETERS
T2_ncdm = 0.69054 * (qs.*m0s/mp).^(1/2) .* sqrt.((1 .- qs.^3)./qs.^3 .+ 12.567.^(ks) .* (qs.*m0s/mp).^(2*ks) * 3)
p2_rms = 3.23280617733076
Ω2 = 0.12 * 2 * qs.^2 ./ (1 .+ qs.^2)

T1_ncdm = 0.69054 * (1 .- qs.^3).^(1/2) .*(m0s/mp).^(1/2)
p1_rms = 5.057163286675063
Ω1 = 0.12 * (1 .- qs.^2) ./ (1 .+ qs.^2)

m_ncdms = 1e-6*1.75*5.3^(4/3).*sqrt.(Ω1/0.12*p1_rms^2 .*T1_ncdm.^2 .+ Ω2/0.12*p2_rms^2 .*T2_ncdm.^2)

#m_ncdms = 1*ones(20) #DM mass in GeV
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

use_fits = falses(20) #Whether to use a PSD file or an analytical expression

n_e_class = 200 #Number of values for q=E/T for CLASS
q_min_class = 0.01 #Minimum value for q=E/T for CLASS
q_max_class = 100 #Maximum value for q=E/T for CLASS
#qs_max_class = [10,100,1000,10000,100000] 
n_e = 1000 #Number of values for q=E/T for plots and fitting
q_min = 1e-4 #Minimum value for q=E/T for plots and fitting
q_max = 1e4 #Maximum value for q=E/T for plots and fitting
