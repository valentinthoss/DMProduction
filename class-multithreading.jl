###################################################
# STILL MISSING
#
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

include("Code/units.jl") #Units and definitions
include("Code/data.jl") #Data from BlackHawk, Constraints and others
include("Code/interpolation.jl") #Routines for interpolation
include("Code/cosmology.jl") #Cosmological background
include("Code/technical.jl")
include("Code/evaporation.jl") #Main module which contains the evaporation physics
include("Code/constraints.jl") #Functions to obtain constraints on f_PBH
include("parameters.jl") #Load the parameters for the scan


area = zeros(length(m0s))
    
pbhini = open("PBH.ini")
pbh_ini_string = read(pbhini,String)
close(pbhini)
run(`mkdir -p Plots`)
run(`mkdir -p Output`)

Threads.@threads for i=1:length(m0s)

    parameters = []
    
    m0 = m0s[i]
    q = qs[i]
    k = ks[i]
    m_ncdm = m_ncdms[i]
    burst = bursts[i]
    use_fit = use_fits[i]

    println("m=",m0,"g, q=",q,", k=",k,", burst=",burst,", m_dm=",m_ncdm," GeV")
    
    T0 = bh_temp(m0)

    e = 10 .^ LinRange(log10(T0*q_min),log10(T0*q_max),n_e)
    e_class = 10 .^ LinRange(log10(T0*q_min_class),log10(T0*q_max_class),n_e_class)

    t_start = bh_formation(m0)
    if k > 0
        t_end = bh_lifetime_estimate(m0,a_bh,m_fin=q*m0)
    else
        t_end = bh_lifetime_estimate(m0,a_bh)
    end
    #t_end = bh_lifetime_estimate(m0,a_bh,k=k,burst=burst,q=q)
    
    β_c = t_to_a(t_start)/t_to_a(t_end)
    
    if isnothing(md)
        if !isnothing(β)
            β_comp = β
        elseif !isnothing(f)
            β_comp = β_fm(f,m0)
        else
            if burst
                N_approx = 2.85e-2*m0^2/mp^2
            else
                N_approx = 2.85e-2*m0^2/mp^2 * (1+q^2)
            end
            f_comp = Ω_ncdm * m0 / (N_approx * m_ncdm * GeV_to_g * Ω_dm)
            β_comp = β_fm(f_comp,m0)
        end
        if β_comp > β_c
            md_comp = true
        else
            md_comp = false
        end
    else
        md_comp = md
        f_comp = Ω_ncdm * m0 / (N_approx * m_ncdm * GeV_to_g * Ω_dm)
    end
    
    a_end = t_to_a(t_end)
    if md_comp
        a_end = a_end*(9/16)^(1/4)
    end
    
    T_ncdm = T0*GeV_to_K/2.7255*a_end
    
    Φ1, Φ2 = spectrum_redshift(e,m0,s,g,a_bh,
        t_start=t_start,t_end=t_end,go=false,k=k,burst=burst,q=q,instant=false,n_steps=1000,split_flux=true,md=md_comp)
    Φ1_class, Φ2_class = spectrum_redshift(e_class,m0,s,g,a_bh,
        t_start=t_start,t_end=t_end,go=false,k=k,burst=burst,q=q,instant=false,n_steps=1000,split_flux=true,md=md_comp)
    
    N1 = trapz(e,Φ1)
    
    if !isnothing(β)
        if β < β_c
            Ω1 = N1 * m_ncdm * GeV_to_g * Ω_dm * f_βm(β,m0) / m0
            if k > 0
                Ω1_approx = 2.85e-2 * m0 * m_ncdm*GeV_to_g / mp^2 * Ω_dm * f_βm(β,m0) * (1 - q^2)
            else
                Ω1_approx = 2.85e-2 * m0 * m_ncdm*GeV_to_g / mp^2 * Ω_dm * f_βm(β,m0)
            end
            #println("Deviation of Ω_ncdm ",round(100*abs(1-Ω1_approx/Ω1),sigdigits=4)," %")
        else
            Ω1 = Ω_dm * m_ncdm*1e3 / sqrt(m0/(1.1e7*mp)) * 2.85/3.2
            if k > 0
                Ω1 *= (1-q^2)
            end
        end
    elseif !isnothing(f)
        Ω1 = N1 * m_ncdm * GeV_to_g * Ω_dm * f / m0
        if k > 0
            Ω1_approx = 2.85e-2 * m0 * m_ncdm*GeV_to_g / mp^2 * Ω_dm * f * (1 - q^2)
        else
            Ω1_approx = 2.85e-2 * m0 * m_ncdm*GeV_to_g / mp^2 * Ω_dm * f
        end
        #println("Deviation of Ω_ncdm ",round(100*abs(1-Ω1_approx/Ω1),sigdigits=4)," %")
    elseif burst
        if k > 0
            Ω1 = Ω_ncdm * (1-q^2)
        else 
            Ω1 = Ω_ncdm
        end
        f_check = Ω1 * m0 / (N1 * m_ncdm * GeV_to_g * Ω_dm)
        if abs(1-f_check/f_comp) > 1e-3
            throw(ArgumentError("Something with the PBH abundance is not right..."))
        end
    else
        if k > 0
            Ω1 = Ω_ncdm * (1-q^2)/(1+q^2)
        else
            Ω1 = Ω_ncdm
        end
        f_check = Ω1 * m0 / (N1 * m_ncdm * GeV_to_g * Ω_dm)
        if abs(1-f_check/f_comp) > 1e-3
            throw(ArgumentError("Something with the PBH abundance is not right..."))
        end
    end

    if use_fit
        println("Fit cannot be used with multithreading at the moment...")
        fit = [[0,0,0,0,0]]
    else
        fit = [[0,0,0,0,0]]
    end
    
    writedlm(string("Input/psd_",i,"_sc.dat"),[e_class/T0 (T0^3/(mp*g_to_GeV)^2)*Φ1_class./ (e_class/T0).^2])
    
    if two_stages && k > 0
        
        T0 = T0/q

        t_end2 = bh_lifetime_estimate(m0,a_bh,k=k,burst=burst,q=q)

        a_end2 = t_to_a(t_end2)
        if md_comp
            a_end2 = a_end2*(9/16)^(1/4)
        end

        T_ncdm2 = T0*GeV_to_K/2.7255*a_end2

        Φ1, Φ2 = spectrum_redshift(e,m0,s,g,a_bh,
            t_start=t_start,t_end=t_end2,go=false,k=k,burst=burst,q=q,instant=false,n_steps=1000,split_flux=true,md=md_comp)
        Φ1_class, Φ2_class = spectrum_redshift(e_class,m0,s,g,a_bh,
            t_start=t_start,t_end=t_end2,go=false,k=k,burst=burst,q=q,instant=false,n_steps=1000,split_flux=true,md=md_comp)
        
        
        N2 = trapz(e,Φ2)

        if !isnothing(β)
            if β < β_c
                Ω2 = N2 * m_ncdm * GeV_to_g * 0.26 * f_βm(β,m0) / m0
                Ω2_approx = 2.85e-2 * m0 * m_ncdm*GeV_to_g / mp^2 * 0.26 * f_βm(β,m0) * q^2
                if !burst
                    Ω2_approx *= 2
                end
                #println("Deviation of Ω_ncdm ",round(100*abs(1-Ω1_approx/Ω1),sigdigits=4)," %")
            else
                Ω2 = Ω_dm * m_ncdm*1e3 / sqrt(m0/(1.1e7*mp)) * 2.85/3.2 * q^2
                if !burst
                    Ω2 *= 2
                end
            end
        elseif !isnothing(f)
            Ω2 = N2 * m_ncdm * GeV_to_g * 0.26 * f / m0
            Ω2_approx = 2.85e-2 * m0 * m_ncdm*GeV_to_g / mp^2 * 0.26 * f * q^2
            if !burst
                Ω2_approx *= 2
            end
            #println("Deviation of Ω_ncdm ",round(100*abs(1-Ω1_approx/Ω1),sigdigits=4)," %")
        elseif burst
            Ω2 = Ω_ncdm * q^2
            #Ω2 = Ω_ncdm
        else
            Ω2 = 2 * Ω_ncdm * q^2/(1+q^2)
            #Ω2 = Ω_ncdm
        end
        
        if use_fit
            println("Fit cannot be used with multithreading at the moment...")
            fit2 = [[0,0,0,0,0]]
        else
            fit2 = [[0,0,0,0,0]]
        end
        
        writedlm(string("Input/psd_",i,"_mb.dat"),[e_class/T0 (T0^3/(mp*g_to_GeV)^2)*Φ2_class ./(e_class/T0).^2])
    
        push!(parameters,"NNCDM"=>"2")
        push!(parameters,"DEGNCDM"=>"1,1")
        push!(parameters,"MNCDM"=>string(m_ncdm*1e9,",",m_ncdm*1e9))
        push!(parameters,"TNCDM"=>string(T_ncdm,",",T_ncdm2))
        push!(parameters,"OMEGANCDM"=>string(Ω1,",",Ω2))
        push!(parameters,"OMEGACDM"=>string(Ω_dm-Ω1-Ω2))
        push!(parameters,"FITPARAMS"=>string(fit[1][1],",",fit[1][2],",",fit[1][3],",",fit[1][4],",",fit[1][5],
    ",",fit2[1][1],",",fit2[1][2],",",fit2[1][3],",",fit2[1][4],",",fit2[1][5]))
        push!(parameters,"PSDFILE"=>string("../Input/psd_",i,"_sc.dat,../Input/psd_",i,"_mb.dat"))
        if use_fit
            push!(parameters,"USEFILE"=>"1,1")
        else
            push!(parameters,"USEFILE"=>"1,1")
        end
    else
        push!(parameters,"NNCDM"=>"1")
        push!(parameters,"DEGNCDM"=>"1")
        push!(parameters,"MNCDM"=>string(m_ncdm*1e9))
        push!(parameters,"TNCDM"=>string(T_ncdm))
        push!(parameters,"OMEGANCDM"=>string(Ω1))
        push!(parameters,"OMEGACDM"=>string(Ω_dm-Ω1))
        push!(parameters,"FITPARAMS"=>string(fit[1][1],",",fit[1][2],",",fit[1][3],",",fit[1][4],",",fit[1][5]))
        push!(parameters,"PSDFILE"=>string("../Input/psd_",i,"_sc.dat"))
        if use_fit
            push!(parameters,"USEFILE"=>"1")
        else
            push!(parameters,"USEFILE"=>"1")
        end
    end

    push!(parameters,"OUTPUTNAME"=>string(i))
    
    println("Running CLASS...")

    pbh_ini_string_mod = replace(pbh_ini_string,parameters...)
    open(string("Input/PBH_",i,".ini"), "w") do file
        write(file,pbh_ini_string_mod)
    end

    run(`./run.sh $(i)`)

    data_cdm = readdlm("Data/CDM00_pk.dat",skipstart=4)
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
end

writedlm("Output/area_crit.csv",[m0s ks qs bursts m_ncdms use_fits area])