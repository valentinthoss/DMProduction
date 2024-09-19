function flux_exgb2(e,m_init;k=0,burst=true,q=0.5,n_steps=1000,sec=false,h="select",e_t=nothing)
    
    args = (m_init=m_init,k=k,burst=burst,q=q,sec=sec,h=h)
    
    n = 1/7.99e-29*(m_init/2e33)^(-3/2)/(3.086e27)^3
    m_rec = bh_mass_t(m_init,t_rec,k=k,burst=burst,q=q)
    m_fin = bh_mass_t(m_init,t_0,k=k,burst=burst,q=q)
    
    if m_rec == 0.0
        if isa(e,Array)
            return zeros(length(e))
        else
            return 0.0
        end
    end
    
    if q < 1
        a_half = minimum([1,t_to_a_fast(minimum([bh_lifetime_estimate(m_init,m_fin=m_init*q),2*t_0]))])
    else
        a_half = 0.0
    end
    a_fin = minimum([1,t_to_a_fast(minimum([bh_lifetime_estimate(m_init,m_fin=1,k=k,burst=burst,q=q),2*t_0]))])
    
    Φ = zeros(length(e))
    
    if k==0.0
        a_half = a_fin
        q = 0.0
    end

    if a_half > a_rec

        a_space = 10 .^ LinRange(log10(a_rec),log10(a_half),n_steps)
        t_space = a_to_t_fast.(a_space)
        M_space = bh_mass_t_vec(m_init,t_space)
        a_i = [searchsortedfirst(a_space./(7.527e-15*e[i]*M_space),10) for i=1:length(e)]
        e_i = minimum([searchsortedfirst(a_i,length(a_space)),length(e)])
        
        if !isnothing(e_t)
            e_i = minimum([searchsortedfirst(e,e_t),length(e)])
        end

        Threads.@threads for i=1:e_i
            spec = spectrum_particle_e(e[i]./a_space,M_space,1.0,2.0,sec=sec,h=h)
            H = H_a.(a_space)
            y = n*spec./H*c/(4*pi)./a_space.^2
            Φ[i] += abs(integrate(a_space,y,Trapezoidal()))
        end

        M_space = 10 .^ LinRange(log10(Float64(maximum([m_init*q,m_fin+1]))),log10(Float64(m_rec)),n_steps)
        M_space[1] = maximum([m_init*q,m_fin+1])
        M_space[end] = m_rec
        if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
            println("Error")
        end

        t_space = broadcast(x->bh_lifetime_estimate(m_init,m_fin=x),M_space)
        a_space = t_to_a_fast.(t_space)

        Threads.@threads for i=e_i+1:length(e)
            spec = spectrum_particle_e(e[i]./a_space,M_space,1.0,2.0,sec=sec,h=h)./mass_loss_rate.(M_space)
            y = n*c/(4*pi)*spec./a_space
            Φ[i] += abs(integrate(M_space,y,Trapezoidal()))
        end
        
        m_half = m_init*q

    else
        a_half = a_rec
        m_half = m_rec
    end

    if a_half < 1.0 && k > 0.0
        
        if abs(m_half-m_fin)/(m_half) < 1e-3
        
            a_space = 10 .^ LinRange(log10(a_half),log10(a_fin),n_steps)
            for i=1:length(e)
                spec = spectrum_particle_e(e[i]./a_space,m_half,1.0,2.0;args...)
                H = H_a.(a_space)
                y = n*spec./H*c/(4*pi)./a_space.^2
                Φ[i] += abs(integrate(a_space,y,Trapezoidal()))
            end
        else

            a_space = 10 .^ LinRange(log10(a_half),log10(a_fin),n_steps)
            t_space = a_to_t_fast.(a_space)
            M_space = bh_mass_t_vec(m_init,t_space,k=k,burst=burst,q=q)
            if burst
                a_i = [searchsortedfirst(a_space./x_em(e[i],M_space),10) for i=1:length(e)]
                e_i = minimum([searchsortedfirst(a_i,length(a_space)),length(e)])
            else
                a_i = [searchsortedfirst(a_space/x_em(e[i],m_init*q),10) for i=1:length(e)]
                e_i = minimum([searchsortedfirst(a_i,length(a_space)),length(e)])
            end

            if M_space[1] > m_half
                i = findfirst(x->x<=m_half,M_space)
                M_space = M_space[i:end]
                a_space = a_space[i:end]
            end
            
            if M_space[end] == 0.0
                i = findfirst(x->x==0.0,M_space)
                M_space = M_space[1:i-1]
                a_space = a_space[1:i-1]
            end
            
            Threads.@threads for i=1:e_i
                spec = spectrum_particle_e(e[i]./a_space,M_space,1.0,2.0;args...)
                H = H_a.(a_space)
                y = n*spec./H*c/(4*pi)./a_space.^2
                Φ[i] += abs(integrate(a_space,y,Trapezoidal()))
            end

            M_space = 10 .^ LinRange(log10(Float64(m_fin+1)),log10(Float64(m_half)),n_steps)
            M_space[1] = m_fin+1
            M_space[end] = m_half
            if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
                println("Error")
            end
            
            t_space = broadcast(x->bh_lifetime_estimate(m_init,m_fin=x,k=k,burst=burst,q=q),M_space)
            a_space = t_to_a_fast.(t_space)
            
            Threads.@threads for i=e_i+1:length(e)
                spec = spectrum_particle_e(e[i]./a_space,M_space,1.0,2.0;args...)./mass_loss_rate.(M_space,k=k,burst=burst,q=q,m_init=m_init)
                y = n*c/(4*pi)*spec./a_space
                Φ[i] += abs(integrate(M_space,y,Trapezoidal()))
            end
        end
    end


    if isa(e,Array)
        return Φ
    else
        return Φ[1]
    end
end

function constraint_f_gal_ext_old(m;time="f",k=0,burst=true,q=0.5,sec=false,h="select",mass_function=lognormal,m_range=(0,20,1000),fast=true,rescale=true)

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
    if fast
        f_PBH_mc = [ones(length(m_today),length(flux_gc_data[i][:,1])) for i=1:length(flux_gc_data)]
        f_PBH_tmp = [Inf*ones(length(m),length(flux_gc_data[i][:,1])) for i=1:length(flux_gc_data)]        
    else
        f_PBH = ones(length(flux_gc_data),length(m))
    end

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
            if fast
                f_PBH_mc[i][j,:] = f_gal(flux_gc_data[i],m_today[j],Ω[i],k=k,burst=burst,q=q,sec=sec,h=h,return_all=true,m_init=m_init[j])
            else
                f_PBH[i,j] = f_gal_ext(flux_gc_data[i],m[j],m_today,Ω[i],k=k,burst=burst,q=q,sec=sec,h=h,mass_function=mass_function,evolve=evolve)
            end
        end 
    end

    if fast
        dm = grad_m_m0(m_today,m_init,q,k,burst)
        Threads.@threads for j=1:length(m)
            if !evolve
                dndm0 = mass_function.(m_today,m[j])
                dndm = dndm0
                survival_frac = 1
            else
                dndm0 = mass_function.(m_today,m[j])
                dndm = mass_function.(m_init,m[j]) .* dm
                survival_frac = f0_to_f(1,m_today,dndm0,dndm)
            end
            Threads.@threads for i=1:length(flux_gc_data)
                Threads.@threads for l=1:length(flux_gc_data[i][:,1])
                    if survival_frac > 0.0
                        f_PBH_tmp[i][j,l] = ifelse(rescale,survival_frac,1)*convolve_constraints3(m_today,f_PBH_mc[i][:,l],dndm,i_t,m[j],k)
                    else
                        f_PBH_tmp[i][j,l] = 0.0
                    end
                end
            end
        end

        f_PBH = [minimum(f_PBH_tmp[i][j,:]) for i=1:length(flux_gc_data), j=1:length(m)]

    end
                    
    return f_PBH
end


const cosmo_x = log.(10 .^ LinRange(-20,20,100000))
const cosmo_y = log.(a_to_t.(exp.(cosmo_x)))
a_to_t_interp = linear_interpolation(exp.(cosmo_x),exp.(cosmo_y),extrapolation_bc=Line())
t_to_a_interp = linear_interpolation(exp.(cosmo_y),exp.(cosmo_x),extrapolation_bc=Line())

function a_to_t_fast(a)
    s = log(a)
    if s > cosmo_x[1] && s < cosmo_x[end]
        i = searchsortedfirst(cosmo_x,s)
        return exp(cosmo_y[i-1]+(s-cosmo_x[i-1])/(cosmo_x[i]-cosmo_x[i-1])*(cosmo_y[i]-cosmo_y[i-1]))
    elseif s <= cosmo_x[1]
        return exp(cosmo_y[1])*(a/exp(cosmo_x[1]))^2
    else
        return (s - cosmo_x[end])/H0_fid/sqrt(Ω_Λ_fid) + exp(cosmo_y[end])
    end
end

function t_to_a_fast(t)
    s = log(t)
    if s > cosmo_y[1] && s < cosmo_y[end]
        i = searchsortedfirst(cosmo_y,s)
        return exp(cosmo_x[i-1]+(s-cosmo_y[i-1])/(cosmo_y[i]-cosmo_y[i-1])*(cosmo_x[i]-cosmo_x[i-1]))
    elseif s <= cosmo_y[1]
        return exp(cosmo_x[1])*(t/exp(cosmo_y[1]))^(1/2)
    elseif s >= cosmo_y[end]
        return exp.((t-exp(cosmo_y[end]))*H0_fid*sqrt(Ω_Λ_fid)+cosmo_x[end])
    end
end

function flux_exgb2(e,m_init;k=0,burst=true,q=0.5,n_steps=1000,sec=false,h="select",a_rel=0.01)
    
    args = (m_init=m_init,k=k,burst=burst,q=q,sec=sec,h=h)
    
    n = 1/7.99e-29*(m_init/2e33)^(-3/2)/(3.086e27)^3
    m_rec = bh_mass_t(m_init,t_rec,k=k,burst=burst,q=q)
    m_fin = bh_mass_t(m_init,t_0,k=k,burst=burst,q=q)
    
    if m_rec == 0.0
        if isa(e,Array)
            return zeros(length(e))
        else
            return 0.0
        end
    end
    
    if q < 1
        t_half = bh_lifetime_estimate(m_init,m_fin=m_init*q)
        if t_half > t_rec/2
            a_half = minimum([1,t_to_a_fast(minimum([t_half,2*t_0]))])
        else
            a_half = 0.0
        end
    else
        a_half = 0.0
    end
    a_fin = minimum([1,t_to_a_fast(minimum([bh_lifetime_estimate(m_init,m_fin=1,k=k,burst=burst,q=q),2*t_0]))])
    
    Φ = zeros(length(e))
    
    if k==0.0
        a_half = a_fin
        q = 0.0
    end

    if a_half > a_rec
        
        if abs(m_rec-maximum([m_init*q,m_fin]))/(m_rec) < 1e-1

            a_space = 10 .^ LinRange(log10(a_rec),log10(a_half),n_steps)
            t_space = a_to_t_fast.(a_space)
            M_space = bh_mass_t_vec(m_init,t_space)

            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i]./a_space,M_space,1.0,2.0,sec=sec,h=h)
                H = H_a.(a_space)
                y = n*spec./H*c/(4*pi)./a_space.^2
                Φ[i] += abs(trapz(a_space,y))
            end

        else    
        
            M_space = 10 .^ LinRange(log10(Float64(maximum([m_init*q,m_fin+1]))),log10(Float64(m_rec)),n_steps)
            M_space[1] = maximum([m_init*q,m_fin+1])
            M_space[end] = m_rec
            if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
                println("Error")
            end

            t_space = broadcast(x->bh_lifetime_estimate(m_init,m_fin=x),M_space)
            a_space = t_to_a_fast.(t_space)
            ind = searchsortedfirst(a_space[1:end-1]./a_space[2:end] .- 1 , a_rel)
            
            if ind <= n_steps
                n_ext = Int(ceil(log(a_space[ind]/a_space[end])/log(1 + a_rel)+1))
                a_ext = 10 .^ LinRange(log10(a_space[ind]),log10(a_space[end]),n_ext)
                a_space = vcat(a_space[1:ind-1],a_ext)
                M_ext = bh_mass_t_vec(m_init,a_to_t_fast.(a_ext),k=k,burst=burst,q=q)
                M_space = vcat(M_space[1:ind-1],M_ext)
            end


            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i]./a_space,M_space,1.0,2.0,sec=sec,h=h)./mass_loss_rate.(M_space)
                y = n*c/(4*pi)*spec./a_space
                Φ[i] += abs(trapz(M_space,y))
            end
        end
        m_half = m_init*q

    else
        a_half = a_rec
        m_half = m_rec
    end

    if a_half < 1.0 && k > 0.0
        
        if abs(m_half-m_fin)/(m_half) < 1e-3
        
            a_space = 10 .^ LinRange(log10(a_half),log10(a_fin),n_steps)
            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i]./a_space,m_half,1.0,2.0;args...)
                H = H_a.(a_space)
                y = n*spec./H*c/(4*pi)./a_space.^2
                Φ[i] += abs(trapz(a_space,y))
            end
            
        elseif abs(m_half-m_fin)/(m_half) < 1e-1

            a_space = 10 .^ LinRange(log10(a_half),log10(a_fin),n_steps)
            t_space = a_to_t_fast.(a_space)
            M_space = bh_mass_t_vec(m_init,t_space,k=k,q=q,burst=burst)

            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i]./a_space,M_space,1.0,2.0;args...)
                H = H_a.(a_space)
                y = n*spec./H*c/(4*pi)./a_space.^2
                Φ[i] += abs(trapz(a_space,y))
            end

        else   

            M_space = 10 .^ LinRange(log10(Float64(m_fin+1)),log10(Float64(m_half)),n_steps)
            M_space[1] = m_fin+1
            M_space[end] = m_half
            if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
                println("Error")
            end
            
            t_space = broadcast(x->bh_lifetime_estimate(m_init,m_fin=x,k=k,burst=burst,q=q),M_space)
            a_space = t_to_a_fast.(t_space)
            ind = searchsortedfirst(a_space[1:end-1]./a_space[2:end] .- 1 , a_rel)
            
            if ind <= n_steps
                n_ext = Int(ceil(log(a_space[ind]/a_space[end])/log(1 + q)+1))
                a_ext = 10 .^ LinRange(log10(a_space[ind]),log10(a_space[end]),n_ext)
                a_space = vcat(a_space[1:ind-1],a_ext)
                M_ext = bh_mass_t_vec(m_init,a_to_t_fast.(a_ext),k=k,burst=burst,q=q)
                M_space = vcat(M_space[1:ind-1],M_ext)                
            end
            
            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i]./a_space,M_space,1.0,2.0;args...)./mass_loss_rate.(M_space,k=k,burst=burst,q=q,m_init=m_init)
                y = n*c/(4*pi)*spec./a_space
                Φ[i] += abs(trapz(M_space,y))
            end
        end
    end


    if isa(e,Array)
        return Φ
    else
        return Φ[1]
    end
end

function spectrum_redshift_ext(e,m,s=1.0,g=2.0,a_bh=0;t_start=t_rec,t_end=t_0,
    a_start=nothing,a_end=nothing,k=0,burst=true,q=0.5,n_steps=1000,sec=false,
    h="select",a_rel=0.1,go=false,instant=false,split_flux=false,rd=false,
    mass_function=lognormal(),m_range=(0,20,1000))

#args = (m_init=m_init,k=k,burst=burst,q=q,sec=sec,h=h,go=go)
                        
#tail = false

#if instant
#    return spectrum_redshift_instant(e,m_init,s,g,a_bh;t_start=t_start,t_end=t_end,a_start=a_start,a_end=a_end,n_steps=n_steps,a_rel=a_rel,k=k,burst=burst,q=q,sec=sec,h=h,go=go,split_flux=split_flux,rd=rd)
#end

#m_start, a_bh_start = bh_mass_t(m_init,t_start,a_bh,k=k,burst=burst,q=q,return_spin=true)
#m_fin, a_bh_fin = bh_mass_t(m_init,t_end,a_bh,k=k,burst=burst,q=q,return_spin=true)

#if m_fin < m_init*1e-5
#    t_fin = bh_lifetime_estimate(m_init,a_bh,k=k,burst=burst,q=q)
#    a_fin = t_to_a_fast(t_fin)
#    if !sec && burst
#        tail = true
#        e_tail = 1e17/(q*m_init)/(t_to_a_fast(t_end)/a_fin)
#        m_fin = 1e13/(e_tail*t_to_a_fast(t_end)/a_fin)
#    else
#        m_fin = minimum([1e10/(maximum(e)*t_to_a_fast(t_end)/a_fin),m_init*1e-5])
#    end
#else
t_fin = t_end
#end

if isnothing(a_start)
    a_start = t_to_a_fast(t_start)
end

if isnothing(a_end)
    a_end = t_to_a_fast(t_end)
    a_fin = t_to_a_fast(t_fin)
    
end

#if m_start == 0.0
#    if isa(e,Array)
#        return zeros(length(e))
#    else
#        return 0.0
#    end
#end

#if k > 0 && m_fin <= m_init*q
#    if m_start > m_init*q
#        sc_phase = true
#        burden_phase = true
#        m_half = m_init*q
#        t_half = bh_lifetime_estimate(m_init,a_bh,m_fin=m_init*q)
#        if a_bh > 0.0
#            a_bh_half = bh_mass_t(m_init,t_half,a_bh,return_spin=true)[2]
#        else
#            a_bh_half = 0.0
#        end
#        a_half = t_to_a_fast(t_half)
#    else
#        m_half = m_start
#        a_bh_half = a_bh_start
#        a_half = a_start
#        t_half = t_start
#        sc_phase = false
#        burden_phase = true
#    end
#else
#    m_half = m_fin
#    a_bh_half = a_bh_fin
#    t_half = t_fin
#    sc_phase = true
#    burden_phase = false
#end

sc_phase = true
burden_phase = false

#m_today = 10 .^ LinRange(m_range[1],m_range[2],m_range[3])

Φ1_temp = zeros(length(e),n_steps)
Φ1 = zeros(length(e))
#Φ2 = zeros(length(e))
        
#dndm0 = pdf(mass_function,m_init,m)

if sc_phase

    if true
        t_space = 10 .^ LinRange(log10(t_start),log10(t_end),n_steps)
        a_space = t_to_a_fast.(t_space)
        Threads.@threads for j=1:length(t_space)
            ff(x) = bh_mass_t(x/q,t_space[j]) - x
            m_k = find_zero(ff,(1e10,1e20))
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
            m_init = get_m_init.(m_today,t_space[j],k=k,q=q,burst=burst)
            dm = grad_m_m0(m_today,m_init,q,k,burst)
            dndm = pdf(mass_function,m_init,m) .* dm
            #println(trapz(m_today,dndm))
            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i] * a_end / a_space[j],m_today,s,g,m_init=m_init,sec=sec,h=h,go=go,k=k,q=q,burst=burst)
                if isnothing(i_t)
                    spec_int = trapz(m_today,dndm.*spec)
                else
                    spec_int = trapz(m_today[1:i_t-1],dndm[1:i_t-1].*spec[1:i_t-1]) + trapz(m_today[i_t:end],dndm[i_t:end].*spec[i_t:end])
                end
                Φ1_temp[i,j] += abs(spec_int*a_end/a_space[j])
            end
        end
        
        Threads.@threads for i=1:length(e)
            Φ1[i] = trapz(t_space,Φ1_temp[i,:])
        end
        
        return Φ1
    else

        if tail
            i_tail = minimum([maximum([2,searchsortedfirst(e,e_tail)]),length(e)])
        else
            i_tail = length(e)
        end

        t_crit = bh_lifetime_estimate(m_init,a_bh,m_fin=m_start-0.1*(m_start-m_half))
        t_space = 10 .^ LinRange(log10(t_start),log10(t_crit),n_steps)
        a_space = t_to_a_fast.(t_space)

        if a_bh > 0.0
            M_space, a_bh_space = bh_mass_t_vec(m_init,t_space,a_bh,return_spin=true)
        else
            M_space = bh_mass_t_vec(m_init,t_space)
            a_bh_space = zeros(n_steps)
        end

        Threads.@threads for i=1:i_tail
            spec = spectrum_particle_e(e[i] * a_end ./ a_space,M_space,s,g,a_bh_space,sec=sec,h=h,go=go)
            y = spec*a_end./a_space 
            #Φ1[i] += abs(trapz(t_space,y))
        end
        
        M_space = 10 .^ LinRange(log10(m_half),log10(M_space[end]),n_steps)
        t_space = broadcast(x->bh_lifetime_estimate(m_init,a_bh,m_fin=x),M_space)
        a_space = t_to_a_fast.(t_space)

        Threads.@threads for i=1:i_tail
            spec = spectrum_particle_e(e[i] * a_end ./ a_space,M_space,s,g,a_bh_space,sec=sec,h=h,go=go)./mass_loss_rate.(M_space,a_bh_space)
            y = spec*a_end./a_space 
            #Φ1[i] += abs(trapz(M_space,y))
        end

        Threads.@threads for i=i_tail+1:length(e)
            #Φ1[i] = Φ1[i_tail]*(e[i]/e[i_tail])^(-3)
        end
            
    end
        
end

if burden_phase
    if abs(m_half-m_fin)/(m_half) < 1e-3

        a_space = 10 .^ LinRange(log10(a_half),log10(a_end),n_steps)
        Threads.@threads for i=1:length(e)
            spec = spectrum_particle_e(e[i] * a_end ./ a_space,m_half,s,g,a_bh_half;args...)
            H = H_a.(a_space)
            y = spec*a_end./H./a_space.^2
            Φ2[i] += abs(trapz(a_space,y))
        end
        
    elseif abs(m_half-m_fin)/(m_half) < 1e-1

        a_space = 10 .^ LinRange(log10(a_half),log10(a_end),n_steps)
        t_space = a_to_t_fast.(a_space)
        M_space, a_bh_space = bh_mass_t_vec(m_init,t_space,a_bh,m_init=m_init,k=k,q=q,burst=burst,return_spin=true)

        Threads.@threads for i=1:length(e)
            spec = spectrum_particle_e(e[i] * a_end ./ a_space,M_space,s,g,a_bh_space;args...)
            H = H_a.(a_space)
            y = spec*a_end./H./a_space.^2
            Φ2[i] += abs(trapz(a_space,y))
        end

    else                        
                    
        M_space = 10 .^ LinRange(log10(m_fin),log10(m_half),n_steps)
        M_space[1] = m_fin
        M_space[end] = m_half
        if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
            println("Error2")
        end
        
        t_space = broadcast(x->bh_lifetime_estimate(m_init,a_bh,m_fin=x,k=k,burst=burst,q=q),M_space)
        if a_bh_half > 0.0
            M_space, a_bh_space = bh_mass_t_vec(m_init,reverse(t_space),a_bh,k=k,burst=burst,q=q,return_spin=true)
            M_space = reverse(M_space)
            a_bh_space = reverse(a_bh_space)
            M_space[1] = m_fin
            M_space[end] = m_half
            if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
                println("Error3")
            end
        else
            a_bh_space = zeros(n_steps)
        end
                            
        a_space = t_to_a_fast.(t_space)
        ind = searchsortedfirst(a_space[1:end-1]./a_space[2:end] .- 1 , a_rel)
        
        if ind <= n_steps
            n_ext = Int(ceil(log(a_space[ind]/a_space[end])/log(1 + a_rel)+1))
            a_ext = 10 .^ LinRange(log10(a_space[ind]),log10(a_space[end]),n_ext)
            a_space = vcat(a_space[1:ind-1],a_ext)
            t_ext = reverse(a_to_t_fast.(a_ext))
            t_ext[1] = t_half
            M_ext, a_ext = bh_mass_t_vec(m_init,t_ext,a_bh,k=k,burst=burst,q=q,return_spin=true)
            M_space = vcat(M_space[1:ind-1],reverse(M_ext))
            a_bh_space = vcat(a_bh_space[1:ind-1],reverse(a_ext))
        end
                                                    
        if tail
            i_tail = minimum([maximum([2,searchsortedfirst(e,e_tail)]),length(e)])
        else
            i_tail = length(e)
        end
        
        Threads.@threads for i=1:i_tail
            spec = spectrum_particle_e(e[i] * a_end ./ a_space,M_space,s,g,a_bh_space;args...)./mass_loss_rate.(M_space,a_bh_space,k=k,burst=burst,q=q,m_init=m_init)
            y = spec*a_end./a_space
            Φ2[i] += abs(trapz(M_space,y))
        end
                                                        
        Threads.@threads for i=i_tail+1:length(e)
            Φ2[i] = Φ2[i_tail]*(e[i]/e[i_tail])^(-3)
        end
                                                        
    end
end

if isa(e,Array)
    if !split_flux
        return Φ1 .+ Φ2
    else
        return Φ1, Φ2
    end
else
    if !split_flux
        return Φ1[1] + Φ2[1]
    else
        return Φ1[1], Φ2[1]
    end
end
end

function f_gal_ext(flux_data,m_today,m_spectrum,Ω;k=0,burst=true,q=0.5,sec=false,h="select",mass_function=rect_function,evolve=false)
    e_min = flux_data[:,1] .- flux_data[:,2]
    e_max = flux_data[:,1] .+ flux_data[:,3]
    Φ_gc = flux_data[:,4].*(flux_data[:,3] .+ flux_data[:,2])
    Φ_pbh = zeros(length(m_spectrum),length(e_min))

    Threads.@threads for i=1:length(m_spectrum)
        m_init = get_m_init(m_spectrum[i],t_0,k=k,burst=burst,q=q)
        Threads.@threads for l=1:length(e_min)
            Φ_pbh[i,l] = Ω*spectrum_integrate_fast(LinRange(e_min[l],e_max[l],100),m_spectrum[i],1.0,2.0,
            k=k,burst=burst,q=q,m_init=m_init,sec=sec,h=h)/(4*pi*m_spectrum[i])
        end
    end
    
    if !evolve
        mass_func = mass_function.(m_spectrum,m_today)
    else
        mass_func = mass_function_evolution(m_spectrum,t_0,mass_function,m_today,k=k,burst=burst,q=q)
    end
    Φ_pbh_ext = [trapz(m_spectrum,Φ_pbh[:,i].*mass_func) for i=1:length(e_min)]

    return minimum(Φ_gc ./ Φ_pbh_ext)
end