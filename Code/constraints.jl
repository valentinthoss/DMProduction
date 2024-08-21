function β_fm(f,m;γ=1,g=106.75,h=0.67,Ω_dm=0.26)
    return 7.06e-18*γ^(-0.5)*(h/0.67)^2*(g/106.75)^0.25*f*(m/1e15)^0.5*Ω_dm #Equation 6 arXiv:2002.12778v2
end

function f_βm(β,m;γ=1,g=106.75,h=0.67,Ω_dm=0.26)
    return β/7.06e-18*γ^(0.5)*(h/0.67)^(-2)*(g/106.75)^(-0.25)*(m/1e15)^(-0.5)/Ω_dm #Equation 6 arXiv:2002.12778v2
end

function f_to_f0(f,m_init,m_today)
    if m_today > 0 && f > 0 && !isinf(f)
        return f*m_init/m_today 
    elseif isinf(f)
        return Inf
    elseif m_today == 0 && f == 0
        return Inf
    else
        print("Error f_to_f0")
    end
end

function f0_to_f(f,m_init,m_today)
    return f_to_f0(f,m_today,m_init)
end

function f0_to_f(f0,m,dndm0,dndm)
    return f0*trapz(m,dndm .* m)/trapz(m,dndm0 .* m)
end

function f_to_f0(f,m,dndm0,dndm)
    return f0_to_f(f,m,dndm,dndm0)
end


"""
    gc_from_los(r,l,b)

Computes the galactocentric distance from the line-of-sight distance and the galactic coordinates.

Input:

    r : line-of-sight distance in kpc
    l : Galactic longitude
    b : Galactic latitude
"""
function gc_from_los(r,l,b)
    R_sol = 8.2 #Hawthorn et. al. 2016
    return sqrt(r^2 + R_sol^2 - 2*r*R_sol*cos(l)*cos(b))
end

"""
    nfw(r,rs=28.2,rho0=0.0099/3.2)

Returns the dark matter density (in cgs units) following a NFW profile.

Input:

    r : distance from the galactic center in kpc
    rs : scale-length in kpc
    rho0 : density in cgs
"""
function nfw(r;rs=17,rho0=0.0125*Msol_to_g/(kpc_to_cm)^3*1e9)
    return rho0/(r/rs)/(1+r/rs)^2
end

"""
    density_los_integral(ρ,b_min,b_max,l_min,l_max,r_min_r_max)

Computes the line-of-sight integral over the given PBH density
ρ(r) in a region specified by galactic coordinates (l,b).

Units: 1/cm^2

Input:

    ρ(r) : Dark matter density profile in cgs units at the gc radius r in kpc
    b_min, b_max, l_min, l_max : Angular region in galactic coordinates
    r_min, r_max : Line-of-sight integration boundaries in kpc

"""
function density_los_integral(ρ,b_min,b_max,l_min,l_max,r_min,r_max)
    Ω = solid_angle(l_min,l_max,b_min,b_max)
    integrand(x,p) = ρ(gc_from_los(x[1],x[2],x[3]))*cos(x[3])/(cm_to_kpc)
    low = [Float64(r_min),Float64(l_min),Float64(b_min)]
    up = [Float64(r_max),Float64(l_max),Float64(b_max)]
    prob = IntegralProblem(integrand,low,up)
    sol = solve(prob,HCubatureJL())
    return sol.u/Ω
end
    

"""
    flux_pbh_gal(e,m,s,g,ρ,b_min,b_max,l_min,l_max,r_min_r_max,k=0,m_init=m)

Computes the flux from PBH in the galaxy in a region specified
by galactic coordinates with a specified density profile ρ(r).
Assumes a monochromatic mass spectrum.

Units: 1/(GeV*s*cm^2*sr)

Input:

    e : Energy in GeV
    m : Black hole mass in g
    s : Particle spin
    g : Particle degrees of freedom
    ρ(r) : Function returning the dark matter density of the galaxy at the galactocentric radius r in kpc
    b_min, b_max, l_min, l_max : Angular region in galactic coordinates
    r_min, r_max : Line-of-sight integration boundaries in kpc
    k : Power law index k (optional)
    m_init : Initial black hole mass in g (optional)

"""
function flux_pbh_gal(e,m,s,g,ρ,b_min,b_max,l_min,l_max,r_min,r_max,k=0,m_init=m)
    los_integral = density_los_integral(ρ,b_min,b_max,l_min,l_max,r_min,r_max)
    spectrum = spectrum_particle_e(e,m,s,g,k,m_init)
    return los_integral*spectrum/(4*pi*m)
end
    
"""
    flux_pbh_gal_integrate(e_min,e_max,m,s,g,ρ,b_min,b_max,l_min,l_max,r_min_r_max)

Computes the integrated flux  ∫ dE (dN/dEdtdΩdA) from PBH in
the galaxy in a region specified by galactic coordinates with
a specified density profile ρ(r).
Assumes a monochromatic mass spectrum.

Units: 1/(s*cm^2*sr)

Input:

    e_min : Lower integral bound in GeV
    e_max : Upper integral bound in GeV
    m : Black hole mass in g
    s : Particle spin
    g : Particle degrees of freedom
    ρ(r) : Function returning the dark matter density of the galaxy at the galactocentric radius r in kpc
    b_min, b_max, l_min, l_max : Angular region in galactic coordinates
    r_min, r_max : Line-of-sight integration boundaries in kpc
    k : Power law index k (optional)
    m_init : Initial black hole mass in g (optional)

"""
function flux_pbh_gal_integrate(e_min,e_max,m,s,g,ρ,b_min,b_max,l_min,l_max,r_min,r_max,k=0,m_init=m)
    los_integral = density_los_integral(ρ,b_min,b_max,l_min,l_max,r_min,r_max)
    spectrum = spectrum_integrate_e(e_min,e_max,m,s,g,k,m_init)
    return los_integral*spectrum/(4*pi*m)
end
    

function f_gal(flux_data,m_today,Ω;k=0,burst=true,q=0.5,m_init=nothing,sec=false,h="select",return_all=false)
    if m_today > 0
        if isnothing(m_init)
            m_init = get_m_init(m_today,t_0,k=k,q=q,burst=burst)
        end
        e_min = flux_data[:,1] .- flux_data[:,2]
        e_max = flux_data[:,1] .+ flux_data[:,3]
        Φ_pbh = [Ω*spectrum_integrate_fast(LinRange(e_min[l],e_max[l],100),m_today,1.0,2.0,
            k=k,burst=burst,q=q,m_init=m_init,sec=sec,h=h)/(4*pi*m_today) for l=1:length(e_min)]
        Φ_gc = flux_data[:,4].*(flux_data[:,3] .+ flux_data[:,2])
        
        if !return_all
            ind = argmin(Φ_gc ./ Φ_pbh)
            return (Φ_gc ./ Φ_pbh)[ind]
        else
            return Φ_gc ./ Φ_pbh
        end
    else
        return 0.0
    end
end

function constraint_f_gal_today(m_today;
            k=0,burst=true,q=0.5,m_init=nothing,sec=false,h="select")
    
    f_PBH = ones(length(flux_gc_data),length(m_today))
    
                    
    if isnothing(m_init)       
        m_init = zeros(length(m_today))
        Threads.@threads for i=1:length(m_init)
            m_init[i] = get_m_init(m_today[i],t_0,k=k,burst=burst,q=q)
        end
    end
                            
    Ω = [density_los_integral(nfw,flux_gc_b[i][1],flux_gc_b[i][2],flux_gc_l[i][1],flux_gc_l[i][2],0,20) for i=1:length(flux_gc_data)]
                            
    Threads.@threads for j=1:length(m_today)
        Threads.@threads for i=1:length(flux_gc_data)
            f_PBH[i,j] = f_gal(flux_gc_data[i],m_today[j],Ω[i],k=k,burst=burst,q=q,sec=sec,h=h,m_init=m_init[j])
        end
    end
                              
    return f_PBH
end

function constraint_f_gal_init(m_init;k=0,burst=true,q=0.5,m_today=nothing,sec=false,h="select")

    m_today = zeros(length(m_init))
    Threads.@threads for i=1:length(m_today)
        m_today[i] = bh_mass_t(m_init[i],t_0,k=k,burst=burst,q=q)
    end

    f_PBH = constraint_f_gal_today(m_today,k=k,burst=burst,q=q,m_init=m_init,sec=sec,h=h)

    return f_to_f0.(f_PBH,m_init',m_today')
end

function constraint_f_gal(m;time="f",args...)
    if time == "today" || time == "t"
        return constraint_f_gal_today(m;args...)
    elseif time == "formation" || time == "f"
        return constraint_f_gal_init(m;args...)
    end
end

function constraint_f_gal_ext(m;fast=true,args...)
    if fast
        return constraint_f_gal_ext_fast(m;args...)
    else
        return constraint_f_gal_ext_slow(m;args...)
    end
end

function constraint_f_gal_ext_fast(m;time="f",k=0,burst=true,q=0.5,sec=false,h="select",mass_function=lognormal(),m_range=(0,20,1000),rescale=true)

    if time == "f" || time == "formation"
        evolve = true
    else
        evolve = false
    end

    if q == 0.5
        m_k = 2.7102580628405025e14
    else
        ff(x) = bh_mass_t(x/q,t_0) - x
        m_k = find_zero(ff,(1e14,1e20))
    end

    if k == 0.0 || q == 1.0
        m_today = 10 .^ LinRange(m_range[1],m_range[2],m_range[3])
        i_t = nothing
    else
        n1 = Int(ceil((log10(m_k)-m_range[1])/(m_range[2]-m_range[1])*m_range[3]))
        n2 = Int(ceil((m_range[2]-log10(m_k))/(m_range[2]-m_range[1])*m_range[3]))
        m_today1 = 10 .^ LinRange(m_range[1],log10(m_k),n1)[1:end-1]
        if m_today1[end] >= m_k
            println("Error")
        end
        m_today2 = 10 .^ LinRange(log10(m_k),m_range[2],n2)
        m_today2[1] = m_k
        if m_today2[2] <= m_k
            println("Error2")
        end
        m_today = vcat(m_today1,m_today2)
        i_t = n1
    end

    f_PBH_mc = [ones(length(m_today),length(flux_gc_data[i][:,1])) for i=1:length(flux_gc_data)]
    f_PBH_tmp = [ones(length(m),length(flux_gc_data[i][:,1])) for i=1:length(flux_gc_data)]        


    m_init = zeros(length(m_today))
    Threads.@threads for i=1:length(m_init)
        m_init[i] = get_m_init(m_today[i],t_0,k=k,burst=burst,q=q)
    end

    if k > 0 && q < 1
        i_t2 = findfirst(x-> x>q,m_today ./ m_init)
    else
        i_t2 = nothing
    end

    if i_t !== i_t2
        println("Error3")
    end
                    
    Ω = [density_los_integral(nfw,flux_gc_b[i][1],flux_gc_b[i][2],flux_gc_l[i][1],flux_gc_l[i][2],0,20) for i=1:length(flux_gc_data)]
    
    Threads.@threads for j=1:length(m_today)
        Threads.@threads for i=1:length(flux_gc_data)
            f_PBH_mc[i][j,:] = f_gal(flux_gc_data[i],m_today[j],Ω[i],k=k,burst=burst,q=q,sec=sec,h=h,return_all=true,m_init=m_init[j])
        end 
    end

    if evolve
        dm = grad_m_m0(m_today,m_init,q,k,burst)
    end
    
    Threads.@threads for j=1:length(m)
        if !evolve
            dndm0 = pdf(mass_function,m_today,m[j])
            dndm = dndm0
            survival_frac = 1
            #m_avg = mean(mass_function,m[j])
        else
            #dndm0 = pdf(mass_function,m_today,m[j])
            dndm = pdf(mass_function,m_init,m[j]) .* dm .* m_today ./ m_init
            survival_frac = trapz(m_today,dndm)
            #survival_frac = f0_to_f(1,m_today,dndm0,dndm)
            #m_avg = mean(mass_function,m[j])
            #m_avg = trapz(m_today,m_today .* dndm)/trapz(m_today,dndm)
        end
        Threads.@threads for i=1:length(flux_gc_data)
            Threads.@threads for l=1:length(flux_gc_data[i][:,1])
                if survival_frac > 0.0
                    f_PBH_tmp[i][j,l] = ifelse(rescale,survival_frac,1)*convolve_constraints3(m_today,f_PBH_mc[i][:,l],dndm,i_t,m[j],k)
                else
                    f_PBH_tmp[i][j,l] = Inf
                end
            end
        end
    end

    f_PBH = [minimum(f_PBH_tmp[i][j,:]) for i=1:length(flux_gc_data), j=1:length(m)]
          
    return f_PBH
end

function constraint_f_gal_ext_slow(m;time="f",k=0,burst=true,q=0.5,sec=false,h="select",mass_function=lognormal(),m_range=(0,20,1000),rescale=true)


    if time == "f" || time == "formation"
        evolve = true
    else
        evolve = false
    end

    m_today = 10 .^ LinRange(m_range[1],m_range[2],m_range[3])
    m_init = get_m_init.(m_today,t_0,k=k,q=q,burst=burst)

    Ω = [density_los_integral(nfw,flux_gc_b[i][1],flux_gc_b[i][2],flux_gc_l[i][1],flux_gc_l[i][2],0,20) for i=1:length(flux_gc_data)]

    f_PBH = zeros(length(flux_gc_data),length(m))

    Threads.@threads for i=1:length(m)
        if !evolve
            dndm0 = pdf(mass_function,m_today,m[j])
            dndm = dndm0
            survival_frac = 1
            m_avg = mean(mass_function,m[j])
        else
            dndm0 = pdf(mass_function,m_today,m[j])
            dndm = pdf(mass_function,m_init,m[j]) .* dm
            survival_frac = f0_to_f(1,m_today,dndm0,dndm)
            m_avg = trapz(m_today,m_today .* dndm)/trapz(m_today,dndm)
        end
        Threads.@threads for j=1:length(flux_gc_data)
            flux_data = flux_gc_data[j]
            e_min = flux_data[:,1] .- flux_data[:,2]
            e_max = flux_data[:,1] .+ flux_data[:,3]
            Φ_pbh = zeros(length(e_min))
            for l=1:length(e_min)
                e = LinRange(e_min[l],e_max[l],100)
                spec = [spectrum_particle_e(e[i],m_today[j],1.0,2.0,k=k,burst=burst,q=q,m_init=m_init,sec=sec,h=h) for i=1:length(e), j=1:length(m_today)]
                spec_dndm = [trapz(m_today,spec[i,:] .* dndm) for i=1:length(e)]
                spec_dndm_e = trapz(e,spec_dndm)
                Φ_pbh[l] = Ω[j]*spec_dndm_e/(4*pi*m_avg)
            end
            
            Φ_gc = flux_data[:,4].*(flux_data[:,3] .+ flux_data[:,2])
            

            ind = argmin(Φ_gc ./ Φ_pbh)
            f_PBH[j,i] = (Φ_gc ./ Φ_pbh)[ind]
        end
    end

    return f_PBH
end

function constraint_k_gal_today(m_today,k_space;
        burst=true,q=0.5,m_init=nothing,Ω=nothing,return_f=true,return_k=true,maxiters=5,sec=false,h="select",f_target=1.0,f_PBH=nothing)
    if isnothing(m_init)
        m_init = zeros(length(m_today))
        Threads.@threads for i=1:length(m_init)
            m_init[i] = get_m_init(m_today[i],t_0)
        end
    end

    if isnothing(f_PBH)
        f_PBH = zeros(length(flux_gc_data),length(m_today),length(k_space))
        Threads.@threads for k=1:length(k_space)
            f_PBH[:,:,k] = constraint_f_gal(m_today,k=k_space[k],burst=burst,q=q,time="today",sec=sec,h=h)
        end
    end
    
    
    if isnothing(Ω)
        Ω = zeros(length(flux_gc_data))
        Threads.@threads for i=1:length(flux_gc_data)
            Ω[i] = density_los_integral(nfw,flux_gc_b[i][1],flux_gc_b[i][2],flux_gc_l[i][1],flux_gc_l[i][2],0,20)
        end
    end

    if return_k

        k_PBH = zeros(length(flux_gc_data),length(m_today))

        Threads.@threads for i=1:length(flux_gc_data)
            Threads.@threads for j=1:length(m_today)
                if m_init[j]*q < m_today[j] && q < 1.0
                    if f_PBH[i,j,1] < f_target
                        k_PBH[i,j] = Inf
                    else
                        k_PBH[i,j] = 0.0
                    end
                #elseif burst
                #    k_PBH[i,j] = maximum([0,-log(f_PBH[i,j,1])/log(bh_entropy(m_today[j]))])
                else
                    ind = findlast(x->x<f_target,f_PBH[i,j,:])
                    if !isnothing(ind)
                        if ind > 1 && ind < length(k_space)
                            k_min = k_space[ind]
                            k_max = k_space[ind+1]
                            f(x) = f_gal(flux_gc_data[i],m_today[j],Ω[i],k=x,burst=burst,q=q,m_init=get_m_init(m_today[j],t_0,k=x,burst=burst,q=q),sec=sec,h=h) - f_target
                            k_PBH[i,j] = find_zero(f,(k_min,k_max),maxiters=maxiters)    
                        end
                    end
                end
            end
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

function constraint_k_gal_init(m_init,k_space;
        burst=true,q=0.5,m_today=nothing,Ω=nothing,return_f=true,return_k=true,maxiters=5,sec=false,h="select",f_target=1,f_PBH=nothing)
    
    args = (q=q,burst=burst,sec=sec,h=h)
    
    if isnothing(m_today)
        m_today = zeros(length(m_init))
        Threads.@threads for i=1:length(m_today)
            m_today[i] = bh_mass_t(m_init[i],t_0)
        end
    end
    
    if isnothing(Ω)
        Ω = zeros(length(flux_gc_data))
        Threads.@threads for i=1:length(flux_gc_data)
            Ω[i] = density_los_integral(nfw,flux_gc_b[i][1],flux_gc_b[i][2],flux_gc_l[i][1],flux_gc_l[i][2],0,20)
        end
    end

    if isnothing(f_PBH)
        f_PBH = zeros(length(flux_gc_data),length(m_init),length(k_space))
        Threads.@threads for k=1:length(k_space)
            f_PBH[:,:,k] = constraint_f_gal(m_init;k=k_space[k],time="formation",args...)
        end
    end
    
    if return_k

        k_PBH = zeros(length(flux_gc_data),length(m_init))

        Threads.@threads for i=1:length(flux_gc_data)
            Threads.@threads for j=1:length(m_init)
                if m_init[j]*q < m_today[j] && q < 1.0
                    if f_PBH[i,j,1] < f_target
                        k_PBH[i,j] = Inf
                    else
                        k_PBH[i,j] = 0.0
                    end
                else
                    ind = findlast(x->x<f_target,f_PBH[i,j,:])
                    if !isnothing(ind)
                        if ind > 1 && ind < length(k_space)
                            k_min = k_space[ind]
                            k_max = k_space[ind+1]
                            m(x) = bh_mass_t(m_init[j],t_0,k=x,burst=burst,q=q)
                            f(x) = f_to_f0(f_gal(flux_gc_data[i],m(x),Ω[i];k=x,m_init=m_init[j],args...),m_init[j],m(x)) - f_target
                            k_PBH[i,j] = find_zero(f,(k_min,k_max),maxiters=maxiters)    
                        end
                    end
                end
            end
        end
    end
    
    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

function constraint_k_gal(m,k;time="today",args...)
    if time == "today" || time == "t"
        return constraint_k_gal_today(m,k;args...)
    elseif time == "formation" || time == "f"
        return constraint_k_gal_init(m,k;args...)
    end
end

function constraint_k_gal_ext(m,k_space;time="f",burst=true,q=0.5,m_today=nothing,Ω=nothing,return_f=true,return_k=true,maxiters=5,sec=false,h="select",f_target=1,mass_function=lognormal(),m_range=(0,20,1000),fast=true,rescale=true)

    args = (q=q,burst=burst,sec=sec,h=h,mass_function=mass_function,m_range=m_range,fast=fast,rescale=rescale)

    if time == "f" || time == "formation"
        evolve = true
    else
        evolve = false
    end

    f_PBH = zeros(length(flux_gc_data),length(m),length(k_space))
    Threads.@threads for k=1:length(k_space)
        f_PBH[:,:,k] = constraint_f_gal_ext(m;k=k_space[k],time=time,args...)
    end
    
    if return_k
        k_PBH = zeros(length(flux_gc_data),length(m))

        Threads.@threads for i=1:length(flux_gc_data)
            Threads.@threads for j=1:length(m)
                ind = findlast(x->x<f_target,f_PBH[i,j,:])
                if !isnothing(ind)
                    if ind > 1 && ind < length(k_space)
                        k_min = k_space[ind]
                        k_max = k_space[ind+1]
                        f(x) = constraint_f_gal_ext(m;time=time,k=x,args...)[i,j] - f_target
                        k_PBH[i,j] = find_zero(f,(k_min,k_max),maxiters=maxiters)    
                    end
                else
                    k_PBH[i,j] = Inf
                end
            end
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

function flux_exgb(e,m_init;k=0,burst=true,q=0.5,n_steps=1000,sec=false,h="select",a_rel=0.01)
    n = 1/7.99e-29*(m_init/2e33)^(-3/2)/(3.086e27)^3
    return n*c/(4*pi)*spectrum_redshift(e,m_init,1.0,2.0,0.0,
    t_start=t_rec,t_end=t_0,a_start=a_rec,a_end=1,k=k,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h,a_rel=a_rel)
end

function β_exgb(flux_data,m_init;k=0,burst=true,q=0.5,n_steps=1000,sec=false,h="select",return_all=false)
    e = 10 .^ LinRange(log10(minimum(flux_data[:,1]))-1,log10(maximum(flux_data[:,1]))+1,100)
    Φ = flux_exgb(e,m_init,k=k,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h)
    if sum(Φ) > 0.0
        interp = linear_interpolation(e,Φ)
        integrand(x,p) = interp(x)
        e_min = flux_data[:,1] .- flux_data[:,2]
        e_max = flux_data[:,1] .+ flux_data[:,3]
        prob = [IntegralProblem(integrand,e_min[k],e_max[k]) for k=1:length(e_min)]
        sol = [solve(prob[k],QuadGKJL()).u for k=1:length(e_min)]
        Φ_exgb = flux_data[:,4].*(flux_data[:,3] .+ flux_data[:,2])
        if !return_all
            return minimum([1,minimum(Φ_exgb ./ sol)])
        else
            return Φ_exgb ./ sol
        end
    else
        if !return_all
            return Inf
        else
            return Inf*ones(length(flux_data[:,1]))
        end
    end
end
            
function constraint_β_exgb_init(m_init;k=0,burst=true,q=0.5,n_steps=1000,sec=false,h="select")
    β_PBH = zeros(length(flux_exgb_data),length(m_init))

    Threads.@threads for j=1:length(m_init)
        for i=1:length(flux_exgb_data)
            β_PBH[i,j] = β_exgb(flux_exgb_data[i],m_init[j],k=k,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h)*m_init[j]
        end
    end
    return β_PBH
end
                    
function constraint_β_exgb_today(m_today;k=0,burst=true,q=0.5,n_steps=1000,m_init=nothing,sec=false,h="select")
    if isnothing(m_init)
        m_init = get_m_init.(m_today,t_0,k=k,burst=burst,q=q)
    end
    return constraint_β_exgb_init(m_init,k=k,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h)
end

function constraint_β_exgb_est(m_init;k=0,burst=true,q=0.5)
    args = (k=k,burst=burst,q=q)
    M_t = get_m_init(0,t_0;args...)
    M_r = get_m_init(0,t_rec;args...)
    if m_init < M_r
        return Inf
    elseif m_init < M_t
        m_rec = bh_mass_t(m_init,t_rec;args...)
        a = t_to_a_fast(bh_lifetime_estimate(m_init,k=k,burst=burst,q=q))
        return 1e-18/(a*m_rec^0.5)
    else
        dM = bh_mass_t(m_init,t_rec;args...) - bh_mass_t(m_init,t_0;args...)
        return 1e-33*m_init^(1.5)/dM
    end
end

function constraint_β_exgb(m;k=0,burst=true,q=0.5,n_steps=1000,m_init=nothing,time="today",sec=false,h="select")
    if time == "today" || time == "t"
        return constraint_β_exgb_today(m,k=k,burst=burst,q=q,n_steps=n_steps,m_init=m_init,sec=sec,h=h)
    elseif time == "formation" || time == "f"
        return constraint_β_exgb_init(m,k=k,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h)
    end
end

function constraint_f_exgb_init(m_init;args...)
    f_PBH = zeros(length(flux_exgb_data),length(m_init))

    for j=1:length(m_init)
        for i=1:length(flux_exgb_data)
            f_PBH[i,j] = f_βm(β_exgb(flux_exgb_data[i],m_init[j];args...),m_init[j])
        end
    end
    return f_PBH
end
                    
function constraint_f_exgb_today(m_today;m_init=nothing,args...)                                               
    if isnothing(m_init)
        m_init = get_m_init.(m_today,t_0,k=k,burst=burst,q=q)
    end
    
    f_PBH_0 = constraint_f_exgb_init(m_init;args...)
    return f0_to_f.(f_PBH_0,m_init',m_today')
end

function constraint_f_exgb(m;k=0,burst=true,q=0.5,n_steps=1000,m_init=nothing,time="today",sec=false,h="select")

    args = (k=k,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h)

    if time == "today" || time == "t"
        return constraint_f_exgb_today(m;m_init=m_init,args...)
    elseif time == "formation" || time == "f"
        return constraint_f_exgb_init(m;args...)
    end
end

function constraint_f_exgb_ext(m;k=0,burst=true,q=0.5,n_steps=1000,sec=false,h="select",mass_function=lognormal(),m_range=(0,20,1000),rescale=false)

    if q == 0.5
        m_k = 5.420516125681005e14
    else
        ff(x) = bh_mass_t(x,t_0) - x*q
        m_k = find_zero(ff,(1e14,1e20))
    end

    if k == 0.0 || q == 1.0
        m_init = 10 .^ LinRange(m_range[1],m_range[2],m_range[3])
        i_t = nothing
    else
        n1 = Int(ceil((log10(m_k)-m_range[1])/(m_range[2]-m_range[1])*m_range[3]))
        n2 = Int(ceil((m_range[2]-log10(m_k))/(m_range[2]-m_range[1])*m_range[3]))
        m_init1 = 10 .^ LinRange(m_range[1],log10(m_k),n1)[1:end-1]
        if m_init1[end] >= m_k
            println("Error")
        end
        m_init2 = 10 .^ LinRange(log10(m_k),m_range[2],n2)
        m_init2[1] = m_k
        if m_init2[2] <= m_k
            println("Error2")
        end
        m_init = vcat(m_init1,m_init2)
        i_t = n1
    end
    
    m_today = bh_mass_t.(m_init,t_0,k=k,q=q,burst=burst)
    if rescale
        m_init_init = get_m_init.(m_init,t_0,k=k,q=q,burst=burst)
        dm = grad_m_m0(m_init,m_init_init,q,k,burst)
    end

    f_PBH_mc = [ones(length(m_init),length(flux_exgb_data[i][:,1])) for i=1:length(flux_exgb_data)]
    f_PBH_tmp = [ones(length(m),length(flux_exgb_data[i][:,1])) for i=1:length(flux_exgb_data)]

    if k > 0 && q < 1
        i_t2 = findfirst(x-> x>q,m_today ./ m_init)
    else
        i_t2 = nothing
    end

    if i_t !== i_t2
        println("Error3")
    end
    
    Threads.@threads for j=1:length(m_init)
        Threads.@threads for i=1:length(flux_exgb_data)
            f_PBH_mc[i][j,:] = f_βm(β_exgb(flux_exgb_data[i],m_init[j],k=k,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h,return_all=true),m_init[j])
        end 
    end

    Threads.@threads for j=1:length(m)
        dndm0 = pdf(mass_function,m_init,m[j])
        #m_avg = mean(mass_function,m[j])
        if rescale
            dndm = pdf(mass_function,m_init_init,m[j]) .* dm .* m_init ./ m_init_init
            survival_frac = trapz(m_init,dndm)
        else
            survival_frac = 1
        end
        Threads.@threads for i=1:length(flux_exgb_data)
            Threads.@threads for l=1:length(flux_exgb_data[i][:,1])
                f_PBH_tmp[i][j,l] = survival_frac*convolve_constraints3(m_init,f_PBH_mc[i][:,l],dndm0,i_t,m[j],k)
            end
        end
    end

    f_PBH = [minimum(f_PBH_tmp[i][j,:]) for i=1:length(flux_exgb_data), j=1:length(m)]
          
    return f_PBH
end
                                                                            
function constraint_k_exgb_init(m_init,k_space;burst=true,q=0.5,return_f=true,return_k=true,maxiters=5,n_steps=100,sec=false,h="select",f_target=1,f_PBH=nothing)
    
    if isnothing(f_PBH)
        f_PBH = zeros(length(flux_exgb_data),length(m_init),length(k_space))

        Threads.@threads for k=1:length(k_space)
            f_PBH[:,:,k] = constraint_f_exgb(m_init,k=k_space[k],burst=burst,q=q,n_steps=n_steps,time="f",sec=sec,h=h)
        end
    end

    k_PBH = zeros(length(flux_exgb_data),length(m_init))
    Threads.@threads for j=1:length(m_init)
        m_today = bh_mass_t(m_init[j],t_0)
        Threads.@threads for i=1:length(flux_exgb_data)
            if m_init[j]*q < m_today
                if f_PBH[i,j,1] < f_target
                    k_PBH[i,j] = Inf
                else
                    k_PBH[i,j] = 0.0
                end
            else
                ind = findlast(x->x<f_target,f_PBH[i,j,:])
                if !isnothing(ind)
                    if ind < length(k_space) && ind > 1
                        k_min = k_space[ind]
                        k_max = k_space[ind+1]
                        f(x) = β_exgb(flux_exgb_data[i],m_init[j],k=x,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h)/β_fm(1,m_init[j]) - f_target
                        k_PBH[i,j] = find_zero(f,(k_min,k_max),maxiters=maxiters)    
                    elseif ind == length(k_space)
                        k_PBH[i,j] = Inf                                                                          
                    else
                        k_PBH[i,j] = 0.0
                    end
                end
            end
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

function constraint_k_exgb_today(m_today,k_space;
    burst=true,q=0.5,return_f=true,return_k=true,maxiters=5,n_steps=100,m_init=nothing,sec=false,h="select",f_target=1,f_PBH=nothing)

    if isnothing(m_init)
        m_init = zeros(length(m_today))
        Threads.@threads for i=1:length(m_init)
            m_init[i] = get_m_init(m_today[i],t_0)
        end
    end
    
    if isnothing(f_PBH)
        f_PBH = zeros(length(flux_exgb_data),length(m_today),length(k_space))

        Threads.@threads for k=1:length(k_space)
            f_PBH[:,:,k] = constraint_f_exgb(m_today,k=k_space[k],burst=burst,q=q,n_steps=n_steps,time="t",sec=sec,h=h)
        end
    end
    
    k_PBH = zeros(length(flux_exgb_data),length(m_today))
    Threads.@threads for j=1:length(m_today)
        Threads.@threads for i=1:length(flux_exgb_data)
            if m_init[j]*q < m_today[j]
                if f_PBH[i,j,1] < 1
                    k_PBH[i,j] = Inf
                else
                    k_PBH[i,j] = 0.0
                end
            else
                ind = findlast(x->x<f_target,f_PBH[i,j,:])
                if !isnothing(ind)
                    if ind < length(k_space) && ind > 1
                        k_min = k_space[ind]
                        k_max = k_space[ind+1]
                        m(x) = get_m_init(m_today[j],t_0,k=x,burst=burst,q=q)
                        f(x) = f0_to_f(β_exgb(flux_exgb_data[i],m(x),k=x,burst=burst,q=q,n_steps=n_steps,sec=sec,h=h)/β_fm(1,m(x)),m(x),m_today[j]) - f_target
                        k_PBH[i,j] = find_zero(f,(k_min,k_max),maxiters=maxiters)    
                    elseif ind == length(k_space)
                        k_PBH[i,j] = Inf                                                                          
                    else
                        k_PBH[i,j] = 0.0
                    end
                end
            end   
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end       
                                                                                                                         
function constraint_k_exgb(m,k;time="f",args...)
    if time == "today" || time == "t"
        return constraint_k_exgb_today(m,k;args...)
    elseif time == "formation" || time == "f"
        return constraint_k_exgb_init(m,k;args...)
    end
end


function constraint_k_exgb_ext(m,k_space;burst=true,q=0.5,return_f=true,return_k=true,maxiters=5,n_steps=100,sec=false,h="select",f_target=1,mass_function=lognormal(),m_range=(0,20,1000),rescale=false)

    args = (q=q,burst=burst,sec=sec,h=h,n_steps=n_steps,mass_function=mass_function,m_range=m_range,rescale=rescale)

    if time == "f" || time == "formation"
        evolve = true
    else
        evolve = false
    end

    f_PBH = zeros(length(flux_exgb_data),length(m),length(k_space))
    Threads.@threads for k=1:length(k_space)
        f_PBH[:,:,k] = constraint_f_exgb_ext(m;k=k_space[k],args...)
    end
    
    if return_k
        k_PBH = zeros(length(flux_exgb_data),length(m))

        Threads.@threads for i=1:length(flux_exgb_data)
            Threads.@threads for j=1:length(m)
                ind = findlast(x->x<f_target,f_PBH[i,j,:])
                if !isnothing(ind)
                    if ind > 1 && ind < length(k_space)
                        k_min = k_space[ind]
                        k_max = k_space[ind+1]
                        f(x) = constraint_f_exgb_ext(m;k=x,args...)[i,j] - f_target
                        k_PBH[i,j] = find_zero(f,(k_min,k_max),maxiters=maxiters)    
                    end
                else
                    k_PBH[i,j] = Inf
                end
            end
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

dir = @__DIR__
path_exoclass = string(dir,"/../Data/Constraint data/exoclass_data_001.csv")
exoclass_data = DelimitedFiles.readdlm(path_exoclass)

m_rec = [bh_mass_t(exoclass_data[i,1],t_rec;k=exoclass_data[i,2],q=1) for i=1:length(exoclass_data[:,1])]
m0 = 2.3e16
f0 = 5.1e-2
dmdt0 = mass_loss_rate(m0,q=1)
dmdt = [mass_loss_rate(m_rec[i];m_init=exoclass_data[i,1],k=exoclass_data[i,2],q=1) for i=1:length(m_rec)]
f_instant = f0*m_rec/m0 ./ (dmdt/dmdt0)
exoclass_bias_factor = exoclass_data[:,3] ./ f_instant
const exoclass_bias_interp = linear_interpolation(log.(exoclass_data[:,1]),log.(exoclass_bias_factor),extrapolation_bc=Flat())

path_stoecker = string(dir,"/../Data/Constraint data/cmb_exoclass.csv")
stoecker_data = DelimitedFiles.readdlm(path_stoecker,',')
const stoecker_interp = linear_interpolation(log.(stoecker_data[:,1]),log.(stoecker_data[:,2]),extrapolation_bc=Inf)
                                                                                                                            
function constraint_f_aniso_init(m_init;k=0,burst=true,q=0.5,method="instant",m0=2.3e16,f0=5.1e-2,z_ext=10)
    
    args = (k=k,burst=burst,q=q)
    f_PBH = zeros(length(m_init))
    M_r = get_m_init(0,t_rec;args...)
    
    dmdt0 = mass_loss_rate(m0)
    
    if method == "extended" || method == "rescaled"
        t_ext = z_to_t(z_ext)
    end
    
    Threads.@threads for i=1:length(m_init)
        if m_init[i] < M_r
            f_PBH[i] = Inf
        else
            m_rec = bh_mass_t(m_init[i],t_rec;args...)
            if method == "instant"
                #f_PBH[i] = 2e-25*m_rec / mass_loss_rate(m_rec;m_init=m_init[i],args...)
                dmdt = mass_loss_rate(m_rec;m_init=m_init[i],args...)
                f_PBH[i] = f0*m_rec/m0 / (dmdt/dmdt0)
            elseif method == "extended"
                dM = m_rec - bh_mass_t(m_init[i],t_ext;args...)
                if dM > 0
                    f_PBH[i] = f0*dmdt0*(t_ext-t_rec) / dM *m_rec/m0
                    #f_PBH[i] = 3e-9*m_rec / dM
                else
                    f_PBH[i] = Inf
                end
            elseif method == "instant_photon"
                #f_PBH[i] = 2e-25*m_rec / mass_loss_rate(m_rec;m_init=m_init[i],args...) * fM_int(m_rec)/fM_int(1e15)
                dmdt0 = mass_loss_rate(m0)
                dmdt = mass_loss_rate(m_rec;m_init=m_init[i],args...)
                f_PBH[i] = f0*m_rec/m0 / (dmdt/dmdt0) *fM_int(m_rec)/fM_int(m0)
            elseif method == "rescaled_instant"
                #dmdt0 = mass_loss_rate(m0)
                dmdt = mass_loss_rate(m_rec;m_init=m_init[i],args...)
                f_PBH[i] = f0*m_rec/m0 / (dmdt/dmdt0)*exp(exoclass_bias_interp(log(m_init[i])))
            elseif method == "rescaled"
                if m_rec <= q*m_init[i] && k > 0
                    dM = m_rec - bh_mass_t(m_init[i],t_ext;args...)
                    if dM > 0
                        f_PBH[i] = f0*dmdt0*(t_ext-t_rec) / dM *m_rec/m0*exp(exoclass_bias_interp(log(m_init[i])))
                    else
                        f_PBH[i] = Inf
                    end
                elseif m_init[i] >= stoecker_data[1,1]
                    f_PBH[i] = exp(stoecker_interp(log(m_init[i])))
                else
                    f_PBH[i] = stoecker_data[1,2]
                end
            else
                println("Method not understood")
                throw(error())
            end
        end
    end
    return f_PBH
end

function constraint_f_aniso_today(m_today;k=0,burst=true,q=0.5,method="instant",m_init=nothing,m0=2.3e16,f0=5.1e-2,z_ext=10)
    if isnothing(m_init)
        m_init = get_m_init.(m_today,t_0,k=k,burst=burst,q=q)
    end
    f_PBH_0 = constraint_f_aniso_init(m_init,k=k,burst=burst,q=q,method=method,m0=m0,f0=f0,z_ext=z_ext)
    return f0_to_f.(f_PBH_0,m_init,m_today)
end
    
function constraint_f_aniso(m;k=0,burst=true,q=0.5,method="instant",m_init=nothing,time="formation",m0=2.3e16,f0=5.1e-2,z_ext=10)
    if time == "today" || time == "t"
        return constraint_f_aniso_today(m,k=k,burst=burst,q=q,method=method,m_init=m_init,m0=m0,f0=f0,z_ext=z_ext)
    elseif time == "formation" || time == "f"
        return constraint_f_aniso_init(m,k=k,burst=burst,q=q,method=method,m0=m0,f0=f0,z_ext=z_ext)
    end
end

function constraint_f_aniso_ext(m;k=0,burst=true,q=0.5,method="instant",m0=2.3e16,f0=5.1e-2,z_ext=10,mass_function=lognormal(),m_range=(0,20,1000),rescale=false)

    if q == 0.5
        m_k = 5.420516125681005e14
    else
        ff(x) = bh_mass_t(x,t_0) - x*q
        m_k = find_zero(ff,(1e14,1e20))
    end

    if k == 0.0 || q == 1.0
        m_init = 10 .^ LinRange(m_range[1],m_range[2],m_range[3])
        i_t = nothing
    else
        n1 = Int(ceil((log10(m_k)-m_range[1])/(m_range[2]-m_range[1])*m_range[3]))
        n2 = Int(ceil((m_range[2]-log10(m_k))/(m_range[2]-m_range[1])*m_range[3]))
        m_init1 = 10 .^ LinRange(m_range[1],log10(m_k),n1)[1:end-1]
        if m_init1[end] >= m_k
            println("Error")
        end
        m_init2 = 10 .^ LinRange(log10(m_k),m_range[2],n2)
        m_init2[1] = m_k
        if m_init2[2] <= m_k
            println("Error2")
        end
        m_init = vcat(m_init1,m_init2)
        i_t = n1
    end

    m_today = bh_mass_t.(m_init,t_0,k=k,q=q,burst=burst)
    if rescale
        m_init_init = get_m_init.(m_init,t_0,k=k,q=q,burst=burst)
        dm = grad_m_m0(m_init,m_init_init,q,k,burst)
    end

    f_PBH_mc = constraint_f_aniso(m_init,k=k,burst=burst,q=q,method=method,m0=m0,f0=f0,z_ext=z_ext)
    f_PBH = zeros(length(m))

    if k > 0 && q < 1
        i_t2 = findfirst(x-> x>q,m_today ./ m_init)
    else
        i_t2 = nothing
    end

    if i_t !== i_t2
        println("Error3")
    end

    Threads.@threads for j=1:length(m)
        dndm0 = pdf(mass_function,m_init,m[j])
        #m_avg = mean(mass_function,m[j])
        if rescale
            dndm = pdf(mass_function,m_init_init,m[j]) .* dm .* m_init ./ m_init_init
            survival_frac = trapz(m_init,dndm)
        else
            survival_frac = 1
        end
        f_PBH[j] = survival_frac*convolve_constraints3(m_init,f_PBH_mc,dndm0,i_t,m[j],k)
    end
          
    return f_PBH
end

function constraint_β_aniso_init(m_init;k=0,burst=true,q=0.5,method="instant")
    return β_fm.(constraint_f_aniso_init(m_init,k=k,burst=burst,q=q,method=method),m_init)
end

function constraint_β_aniso_today(m_today;k=0,burst=true,q=0.5,method="instant",m_init=nothing)
    if isnothing(m_init)
        m_init = get_m_init.(m_today,t_0,k=k,burst=burst,q=q)
    end
    return constraint_β_aniso_init(m_init,k=k,burst=burst,q=q,method=method)
end
    
function constraint_β_aniso(m;k=0,burst=true,q=0.5,method="instant",m_init=nothing,time="formation")
    if time == "today" || time == "t"
        return constraint_β_aniso_today(m,k=k,burst=burst,q=q,method=method,m_init=m_init)
    elseif time == "formation" || time == "f"
        return constraint_β_aniso_init(m,k=k,burst=burst,q=q,method=method)
    end
end

function constraint_k_aniso_today(m_today,k_space;burst=true,q=0.5,method="instant",return_f=true,return_k=true,f_target=1)

    f_PBH = zeros(length(m_today),length(k_space))
    k_PBH = zeros(length(m_today))

    Threads.@threads for j=1:length(k_space)
        f_PBH[:,j] = constraint_f_aniso(m_today,k=k_space[j],burst=burst,q=q,time="t",method=method)
    end

    Threads.@threads for i=1:length(m_today)
        ind = findlast(x->x<f_target,f_PBH[i,:])
        if !isnothing(ind)
            if ind < length(k_space)
                k_min = k_space[ind]
                k_max = k_space[ind+1]
                f(x) = constraint_f_aniso(m_today[i],k=x,burst=burst,q=q,time="t",method=method)[1] - f_target
                k_PBH[i] = find_zero(f,(k_min,k_max))    
            else
                k_PBH[i] = Inf
            end
        else
            k_PBH[i] = 0.0
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

function constraint_k_aniso_init(m_init,k_space;burst=true,q=0.5,method="instant",return_f=true,return_k=true,f_target=1)
    
    f_PBH = zeros(length(m_init),length(k_space))
    k_PBH = zeros(length(m_init))

    Threads.@threads for j=1:length(k_space)
        f_PBH[:,j] = constraint_f_aniso(m_init,k=k_space[j],burst=burst,q=q,time="f",method=method)
    end

    Threads.@threads for i=1:length(m_init)
        ind = findlast(x->x<f_target,f_PBH[i,:])
        if !isnothing(ind)
            if ind < length(k_space)
                k_min = k_space[ind]
                k_max = k_space[ind+1]
                f(x) = constraint_f_aniso(m_init[i],k=x,burst=burst,q=q,time="f",method=method)[1] - f_target
                k_PBH[i] = find_zero(f,(k_min,k_max))    
            else
                k_PBH[i] = Inf
            end
        else
            k_PBH[i] = 0.0
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

function constraint_k_aniso(m,k;time="f",args...)
    if time == "today" || time == "t"
        return constraint_k_aniso_today(m,k;args...)
    elseif time == "formation" || time == "f"
        return constraint_k_aniso_init(m,k;args...)
    end
end

function constraint_k_aniso_ext(m,k_space;burst=true,q=0.5,method="instant",return_f=true,return_k=true,f_target=1,mass_function=lognormal(),m_range=(0,20,1000),rescale=false,maxiters=5)

    args = (q=q,burst=burst,method=method,mass_function=mass_function,m_range=m_range,rescale=rescale)

    f_PBH = zeros(length(m),length(k_space))
    Threads.@threads for k=1:length(k_space)
        f_PBH[:,k] = constraint_f_aniso_ext(m;k=k_space[k],args...)
    end
    
    if return_k
        k_PBH = zeros(length(m))
        Threads.@threads for j=1:length(m)
            ind = findlast(x->x<f_target,f_PBH[j,:])
            if !isnothing(ind)
                if ind > 1 && ind < length(k_space)
                    k_min = k_space[ind]
                    k_max = k_space[ind+1]
                    f(x) = constraint_f_aniso_ext(m;k=x,args...)[j] - f_target
                    k_PBH[j] = find_zero(f,(k_min,k_max),maxiters=maxiters)    
                end
            else
                k_PBH[j] = Inf
            end
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

path_bbn_mb = string(dir,"/../Data/Constraint data/bbn_mb.csv")
bbn_mb_data = DelimitedFiles.readdlm(path_bbn_mb,',')
const bbn_mb_interp = linear_interpolation(log.(bbn_mb_data[:,1]),log.(bbn_mb_data[:,2]),extrapolation_bc=Inf)
 
path_bbn_sc = string(dir,"/../Data/Constraint data/bbn_sc.csv")
bbn_sc_data = DelimitedFiles.readdlm(path_bbn_sc,',')
const bbn_sc_interp = linear_interpolation(log.(bbn_sc_data[:,1]),log.(bbn_sc_data[:,2]),extrapolation_bc=Inf)

function constraint_β_bbn(m_init;burden=false)
    if burden
        return exp.(bbn_mb_interp(log.(m_init)))
    else
        return exp.(bbn_sc_interp(log.(m_init)))
    end
end
    
function constraint_f_bbn(m_init;burden=false)
    return f_βm.(constraint_β_bbn(m_init,burden=burden),m_init)
end

function constraint_f_bbn_ext(m_init;burden=false,mass_function=lognormal(),m_int=(-2,2,100),fast=true)
    return f_βm.(constraint_β_bbn(m_init,burden=burden),m_init)
end

function constraint_f_bbn_ext(m;burden=false,mass_function=lognormal(),m_range=(0,20,1000),rescale=false)

    m_init = 10 .^ LinRange(m_range[1],m_range[2],m_range[3])
    i_t = nothing
    
    if rescale
        m_init_init = get_m_init.(m_init,t_0)
        dm = grad_m_m0(m_init,m_init_init,0.5,0,true)
    end

    f_PBH_mc = constraint_f_bbn(m_init,burden=burden)
    f_PBH = zeros(length(m))

    Threads.@threads for j=1:length(m)
        dndm0 = pdf(mass_function,m_init,m[j])
        #m_avg = mean(mass_function,m[j])
        if rescale
            dndm = pdf(mass_function,m_init_init,m[j]) .* dm .* m_init ./ m_init_init
            survival_frac = trapz(m_init,dndm)
        else
            survival_frac = 1
        end
        f_PBH[j] = survival_frac*convolve_constraints3(m_init,f_PBH_mc,dndm0,i_t,m[j],0.0)
    end
          
    return f_PBH
end

function constraint_k_bbn(m_init,k_space;f_target=1,return_f=true,return_k=true,burden=false)

    f_PBH = constraint_f_bbn(m_init,burden=burden) * ones(length(k_space))'
    k_PBH = zeros(length(m_init))

    Threads.@threads for i=1:length(m_init)
        if f_PBH[i,1] < f_target
            k_PBH[i] = Inf
        else
            k_PBH[i] = 0.0
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end

function constraint_k_bbn_ext(m,k_space;burden=false,return_f=true,return_k=true,f_target=1,mass_function=lognormal(),m_range=(0,20,1000),rescale=false)

    f_PBH = constraint_f_bbn_ext(m,burden=burden,mass_function=mass_function,m_range=m_range,rescale=rescale) * ones(length(k_space))'
    k_PBH = zeros(length(m))
    
    Threads.@threads for i=1:length(m)
        if f_PBH[i,1] < f_target
            k_PBH[i] = Inf
        else
            k_PBH[i] = 0.0
        end
    end

    if return_f && return_k
        return k_PBH, f_PBH
    elseif return_f
        return f_PBH
    elseif return_k
        return k_PBH
    end
end