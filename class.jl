###################################################
# STILL MISSING
#
# - Radiation Domination (beta > beta_c)
# - Two stages (how does it work?)
# - Mass functions
###################################################

using Pkg
Pkg.activate("./")
Pkg.instantiate()

using DelimitedFiles
using Integrals
using Interpolations
using Trapz
using SciPy
using LaTeXStrings

#Loading the code...

include("../Code/units.jl") #Units and definitions
include("../Code/data.jl") #Data from BlackHawk, Constraints and others
include("../Code/interpolation.jl") #Routines for interpolation
include("../Code/cosmology.jl") #Cosmological background
include("../Code/technical.jl")
include("../Code/evaporation.jl") #Main module which contains the evaporation physics
include("../Code/constraints.jl") #Functions to obtain constraints on f_PBH


#BLACK HOLE PARAMETERS

m0s = mp*[1e5,1e5] #BH mass in gram
ks = [0,0] #Burden exponent
qs = [0,0] #Burden q (when memory burden sets in)
bursts = [false,false] #Whether T~const. or T~1/M during memory burden
a_bh = 0.0 #Value of initial BH spin parameter

#PARTICLE DM PARAMETERS

m_ncdms = [0.01,0.01] #DM mass in GeV
s = 0.5 #Spin of DM particle
g = 1 #DOF of DM particle

#COMOLOGICAL PARAMETERS

#Specify one of the following three parameters:
β = nothing #Initial PBH density fraction. If set to nothing, f is used.
f = nothing #Initial PBH dark matter fraciton. If set to nothing then Ω_ncdm is used.
Ω_ncdm = 0.12 #If set to nothing, then it is computed.

Ω_dm = 0.1201075 #Total dark matter density

#PSD PARAMETERS

use_fits = [false,true] #Whether to use a PSD file or an analytical expression

n_e_class = 200 #Number of values for q=E/T for CLASS
q_min_class = 0.01 #Minimum value for q=E/T for CLASS
q_max_class = 100 #Maximum value for q=E/T for CLASS
n_e = 1000 #Number of values for q=E/T for CLASS
q_min = 1e-4 #Minimum value for q=E/T for plots and fitting
q_max = 1e3 #Maximum value for q=E/T for plots and fitting


using PyPlot
rcdefaults()
rc("figure",dpi=300,figsize=(4,3))
rc("font",size=12,family="serif",serif="cmr10")
rc("mathtext",fontset="cm")
rc("axes",unicode_minus=false)
rc("xtick",direction="in",top=true)
rc("ytick",direction="in",right=true)
rc("savefig",bbox="tight")
ion()

area = zeros(length(m0s))

for i=1:length(m0s)
    
    m0 = m0s[i]
    q = qs[i]
    k = ks[i]
    m_ncdm = m_ncdms[i]
    burst = bursts[i]
    use_fit = use_fits[i]
    
    println("m=",m0,"g, q=",q,", k=",k,", burst=",burst,", m_dm=",m_ncdm)
    
    T0 = bh_temp(m0)

    e = 10 .^ LinRange(log10(T0*q_min),log10(T0*q_max),n_e)
    e_class = 10 .^ LinRange(log10(T0*q_min_class),log10(T0*q_max_class),n_e_class)

    t_start = bh_formation(m0)
    t_end = bh_lifetime_estimate(m0,a_bh,m_fin=q*m0)

    a_end = t_to_a(t_end)
    #a_end_approx = 3.33e-31 * (1-q^3)^(1/2)*(m0/mp)^(3/2)
    #println("aev_1 ",abs(1-a_end/a_end_approx))
    
    T_ncdm = T0*GeV_to_K/2.7255*a_end
    #T_ncdm_approx = 0.69054 * (1-q^3)^(1/2) *(m0/mp)^(1/2)
    println("Tncdm ",T_ncdm)

    Φ1, Φ2 = spectrum_redshift(e,m0,s,g,a_bh,
        t_start=t_start,t_end=t_end,go=false,k=k,burst=burst,q=q,instant=false,n_steps=1000,split_flux=true)
    Φ1_class, Φ2_class = spectrum_redshift(e_class,m0,s,g,a_bh,
        t_start=t_start,t_end=t_end,go=false,k=k,burst=burst,q=q,instant=false,n_steps=1000,split_flux=true)

    N1 = trapz(e,Φ1)
    #N1_approx = 2.85e-2*m0^2/mp^2*(1-q^2)
    #println("N_1 ",abs(1-N1/N1_approx))
    
    if !isnothing(β)
        Ω1 = N1 * m_ncdm * GeV_to_g * 0.26 * f_βm(β,m0) / m0
        Ω1_approx = 2.85e-2 * m0 * m_ncdm*GeV_to_g / mp^2 * 0.26 * f_βm(β,m0) * (1 - q^2)
        println("Deviation of Ω_ncdm ",round(100*abs(1-Ω1_approx/Ω1),sigdigits=4)," %")
    elseif !isnothing(f)
        Ω1 = N1 * m_ncdm * GeV_to_g * 0.26 * f / m0
        Ω1_approx = 2.85e-2 * m0 * m_ncdm*GeV_to_g / mp^2 * 0.26 * f * (1 - q^2)
        println("Deviation of Ω_ncdm ",round(100*abs(1-Ω1_approx/Ω1),sigdigits=4)," %")
    elseif burst
        Ω1 = Ω_ncdm * (1-q^2)
    else
        Ω1 = Ω_ncdm * (1-q^2)/(1+q^2)
    end
    
    println("Ω_ncdm ",Ω1)

    ff(x,a,b,c,d,e) = spectrum_redshift_approximation.(x,bh_lifetime_estimate(m0,a_bh,m_fin=q*m0),bh_lifetime_estimate(m0,a_bh,k=k,burst=burst,q=q),a,b,c,d,e)
    fit = optimize.curve_fit(ff,e/T0,Φ1*T0^3/(mp*g_to_GeV)^2,[1,1,1,1,1])
    Φ_fit = ff.(e/T0,fit[1]...)
    
    p_mean = trapz(e/T0,e/T0 .* Φ1)/trapz(e/T0,Φ1)
    p_mean_approx = trapz(e/T0,e/T0  .* Φ_fit)/trapz(e/T0,Φ_fit)
    p_mean1_err = round(100*abs(1-p_mean_approx/p_mean),sigdigits=4)
    println("Deviation of mean momentum: ",p_mean1_err," %")

    p_rms = sqrt(trapz(e/T0,(e/T0).^2 .* Φ1)/trapz(e/T0,Φ1))
    p_rms_approx = sqrt(trapz(e/T0,(e/T0).^2  .* Φ_fit)/trapz(e/T0,Φ_fit))
    p_rms1_err = round(100*abs(1-p_rms_approx/p_rms),sigdigits=4)
    println("Deviation of rms momentum: ",p_rms1_err," %")
    flush(stdout)
    
    fig = figure(figsize=(8,3))
    ax = subplot(121)
    p = loglog(e/T0, (T0^3/(mp*g_to_GeV)^2)*Φ1,"--")
    plot(e_class/T0, (T0^3/(mp*g_to_GeV)^2)*Φ1_class,color=p[1].get_color())
    plot(e/T0,Φ_fit,"-.",color="black")
    ylim(1e-20,1e5)
    xlim(q_min,q_max)
    grid()
    xlabel(L"q = E/T_F")
    title(L"%$(i) | First stage at $t_{\rm ev,1}$")
        
    ax = subplot(122)
    plot(e_class/T0, (T0^3/(mp*g_to_GeV)^2)*Φ1_class)
    plot(e/T0,Φ_fit,"-.",color="black")
    ylim(0,maximum((T0^3/(mp*g_to_GeV)^2)*Φ1_class)*1.3)
    xlim(0,20)
    grid()
    xlabel(L"q = E/T_F")
    show(false)
    savefig(string("Plots/",i,".png"))
    
    writedlm(string("Input/psd_",i,".dat"),[e/T0 (T0^3/(mp*g_to_GeV)^2)*Φ1])
    println("Running CLASS...")

    if use_fit
        println("Using fit")
        run(`./run.sh $(i) $(m_ncdm*1e9) $(Ω1) $(Ω_dm-Ω1) $(T_ncdm) $(fit[1][1]) $(fit[1][2]) $(fit[1][3]) $(fit[1][4]) $(fit[1][5])`)
    else
        println("Using psd file")
        run(`./run.sh $(i) $(m_ncdm*1e9) $(Ω1) $(Ω_dm-Ω1) $(T_ncdm)`)
    end
    
    println("Computing area criterium...")
    data_cdm = readdlm("Output/CDM00_pk.dat",skipstart=4)
    data_pbh = readdlm(string("Output/",i,"_pk.dat"),skipstart=4)
    cdm_interp = linear_interpolation(log.(data_cdm[:,1]),log.(data_cdm[:,2].*data_cdm[:,1]))
    pbh_interp = linear_interpolation(log.(data_pbh[:,1]),log.(data_pbh[:,2].*data_pbh[:,1]));
    k_ev = 10 .^ LinRange(log10(0.5),log10(20),1000)
    cdm_int = zeros(length(k_ev))
    pbh_int = zeros(length(k_ev))
        
    for i=1:length(k_ev)
        k_int = 10 .^ LinRange(log10(k_ev[i]),log10(500),1000)
        cdm_int[i] = trapz(k_int,exp.(cdm_interp(log.(k_int))))
        pbh_int[i] = trapz(k_int,exp.(pbh_interp(log.(k_int))))
    end
    area[i] = 1 - trapz(k_ev,pbh_int ./ cdm_int)/trapz(k_ev,cdm_int ./ cdm_int)
    println("Area: ",area[i])
    println()
    println("-----------------------------------------------")
    println()
    #ylim(1e-20,1e5)
    #xlim(1e-6,1e3)
    #grid()
    #xlabel(L"x = E/T_F")
    #title(L"First stage at $t_{\rm q}$")

    #T0 = T0/q

    #t_start = bh_formation(m0)
    #t_end = bh_lifetime_estimate(m0,a,k=k,burst=burst,q=q)
    #println("t_end ",t_end)

    #a_end2 = t_to_a(t_end)
    #if burst
    #    a_end2_approx = 3.336e-31 * 12.567^(k/2) * q^(3/2+k) * (m0/mp)^(3/2+k) / sqrt(1+2/3*k)
    #else
    #    a_end2_approx = 3.336e-31 * 12.567^(k/2) * q^(3/2+k) * (m0/mp)^(3/2+k) * sqrt(3)
    #end
    #println("aev_2 ",abs(1-a_end2/a_end2_approx))
    
    #T_ncdm2 = T0*GeV_to_K/2.7255*a_end2
    #if burst
    #    T_ncdm2_approx = 0.6905 * 12.567^(k/2) * q^(1/2+k) * (m0/mp)^(1/2 + k) / sqrt(1 + 2/3*k)
    #else
    #    T_ncdm2_approx = 0.6905 * 12.567^(k/2) * q^(1/2+k) * (m0/mp)^(1/2 + k) * sqrt(3)
    #end
    #println("Tncdm_2 ",abs(1-T_ncdm2/T_ncdm2_approx))


    #Φ1, Φ2 = spectrum_redshift(e,m0,s,g,a,
    #    t_start=t_start,t_end=t_end,go=false,k=k,burst=burst,q=q,instant=false,n_steps=1000,split_flux=true)


    #ax = subplot(122)

    #loglog(e/T0, (T0^3/(mp*g_to_GeV)^2)*Φ1)

    #shift = (bh_lifetime_estimate(m0,a,k=k,burst=burst,q=q)/bh_lifetime_estimate(m0,a,m_fin=q*m0))^(1/2)
    #plot(e/T0/shift,1/q^3*shift*Φ_fit,"-.",color="black")

    #loglog(e/T0, (T0^3/(mp*g_to_GeV)^2)*Φ2)

    #ff(x,a,b,c,d,e) = spectrum_redshift_approximation2.(x,bh_lifetime_estimate(m0,a,m_fin=q*m0),bh_lifetime_estimate(m0,a,k=k,burst=burst,q=q),a,b,c,d,e)
    #fit2 = optimize.curve_fit(ff,e/T0,Φ2*T0^3/(mp*g_to_GeV)^2,[1,1,1,1,1])
    #Φ_fit = ff.(e/T0,fit2[1]...)

    #plot(e/T0,Φ_fit,"-.",color="black")

    #p_mean2 = trapz(e/T0,e/T0 .* Φ2)/trapz(e/T0,Φ2)
    #p_mean2_approx = trapz(e/T0,e/T0  .* Φ_fit)/trapz(e/T0,Φ_fit)
    #println("p_mean_2 ",abs(1-p_mean2_approx/p_mean2))
    #p_mean2_err = round(abs(1-p_mean2_approx/p_mean2),sigdigits=4)
    #println("p_rms_2 ",p_mean2_err)

    #p_rms2 = sqrt(trapz(e/T0,(e/T0).^2 .* Φ2)/trapz(e/T0,Φ2))
    #p_rms2_approx = sqrt(trapz(e/T0,(e/T0).^2  .* Φ_fit)/trapz(e/T0,Φ_fit))
    #p_rms2_err = round(abs(1-p_rms2_approx/p_rms2),sigdigits=4)
    #println("p_rms_2 ",p_rms2_err)

    #N2 = trapz(e,Φ2)
    #if burst
    #    N2_approx = 2.85e-2*m0^2/mp^2*q^2
    #else
    #    N2_approx = 2.85e-2*2*m0^2/mp^2*q^2
    #end
    #Ntot_approx = 2.85e-2*m0^2/mp^2*(1+q^2)

    #println("N_2 ",abs(1-N2/N2_approx))
    #println("N_tot ",abs(1-(N1 + N2)/Ntot_approx))
    #T_ncdm1r = round(T_ncdm,sigdigits=5)
    #a_end1r = round(a_end,sigdigits=5)
    #T_ncdm2r = round(T_ncdm2,sigdigits=5)
    #a_end2r = round(a_end2,sigdigits=5)
    #fit1r = round.(fit[1],sigdigits=5)
    #fit2r = round.(fit2[1],sigdigits=5)
    #N1r = round(N1,sigdigits=5)
    #N2r = round(N2,sigdigits=5)

                                #println("MBH,k,q,burst,TNCDM1,aev1,N1,a1,b1,c1,d1,e1,pmean_err,prms_err,TNCDM2,aev2,N2,a2,b2,c2,d2,e2,pmean_err,prms_err")
    #println(m0,",",k,",",q,",",burst,",",T_ncdm1r,",",a_end1r,",",N1r,",",
    #fit1r[1],",",fit1r[2],",",fit1r[3],",",fit1r[4],",",fit1r[5],",",
    #p_mean1_err,",",p_rms1_err,",",T_ncdm2r,",",a_end2r,",",N2r,",",
    #fit2r[1],",",fit2r[2],",",fit2r[3],",",fit2r[4],",",fit2r[5],",",
    #p_mean2_err,",",p_rms2_err)

    #println("MBH,k,q,burst,TNCDM1,aev1,N1,TNCDM2,aev2,N2")
    #println(m0,",",k,",",q,",",burst,",",T_ncdm1r,",",a_end1r,",",N1r,",",
    #T_ncdm2r,",",a_end2r,",",N2r)

    #ylim(1e-20,1e5)
    #xlim(1e-6,1e3)
    #xlabel(L"x = qE/T_F")
    #grid()
    #title(L"First + second stage at $t_{\rm ev}$")
end

writedlm("Output/area_crit.csv",[m0s ks qs bursts m_ncdms use_fits area])


fig = figure()
data_cdm = readdlm("Output/CDM00_pk.dat",skipstart=4)
for i=1:length(m0s)
    data_pbh = readdlm(string("Output/",i,"_pk.dat"),skipstart=4)
    plot(data_pbh[:,1],data_pbh[:,2]./data_cdm[:,2],label=i)
end

xscale("log")
xlim(1e-3,1e3)
grid()
title("Transfer function")
legend()
savefig("Plots/transfer.png")
show(true)