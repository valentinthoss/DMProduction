"""
    bh_entropy(m)

Returns the dimensionless black hole entropy S

Input:

    m : Black hole mass in g
"""
function bh_entropy(m)
    return 2.653e10*m^2
end

"""
    bh_temp(m)

Returns the temperature of black hole in GeV

Input:
    m : Black hole mass in g
"""
function bh_temp(m)
    return 1.0579e13/m
end

"""
    bh_formation(m)

Returns the formation time of a primordial black hole in s.

Input:

    m : Black hole mass in g
"""
function bh_formation(m,γ=1)
    return m*G/γ/c^3 #Carr et. al. 2021, Eq. 4
end


"""
    spectrum_hawking_x(x,s,g)

Returns the interpolated value of dN/dEdt = g/(2*pi*hbar)*Γ/(exp(x)-(-1)^s) in units of 1/(GeV*s)

Input:

    x : E/(4*pi*k_B*T) = 2*G*M*E/(hbar*c^3)
    s : spin
    g : degrees of freedom
"""
function spectrum_hawking_x(x,s,g,a=0;go=false)
    if !go
        return greybody_int(x,s,a)*g/(2*pi*hbar)
    else
        return g/(2*pi*hbar)*27/4*x^2 / (exp(4*pi*x) - (-1)^(2*s))
    end
end

function x_em(e,m)
    return 7.527e-15*m*e
end

"""
    spectrum_hawking_e(e,m,s,g)

Returns the interpolated value of dN/dEdt = g/(2*pi*hbar)*Γ/(exp(x)-(-1)^s) in units of 1/(GeV*s)

Input:

    e : Energy in GeV
    m : Black hole mass in g
    s : spin
    g : degrees of freedom

"""
function spectrum_hawking_e(e,m,s,g,a=0)
    x = x_em(e,m)
    return spectrum_hawking_x(x,s,g,a)
end

"""
    spectrum_particle_e(e,m,s,g,a=0;m_init=0,k=0,burst=true,q=0.5,sec=false,h="select",part_fin="Photon",go=false)

Returns the interpolated value of dN/dEdt = g/(2*pi*hbar)*Γ/(exp(x)-(-1)^s) in units of 1/(GeV*s)
Assumes instantenous memory burden effect.

Input:

    e : Energy in GeV
    m : Black hole mass in g
    s : Particle spin
    g : Particle degrees of freedom
    a : Dimensionless black hole spin parameter (default: 0)
    m_init : Initial black hole mass in g (optional)
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
    sec : If true, secondary emission is included (optional)
    h : Choice of hadronization scheme (optional)
    part_fin : particle species (optional)
    go : If true, then the geometric optics approximation is used (optional)
"""
function spectrum_particle_e(e,m,s,g,a=0;m_init=0,k=0,burst=true,q=0.5,sec=false,h="select",part_fin="Photon",go=false)
    if !burst && k > 0.0
        m = ifelse.(m .> m_init*q,m,m_init*q)
    end
    if !sec
        x = x_em.(e,m)
        Φ = spectrum_hawking_x.(x,s,g,a,go=go)
    else
        x = x_em.(e,m)
        Φ0 = spectrum_hawking_x.(x,s,g,a)
        if isa(m,Array) || isa(m,LinRange)
            Φ1 = spectrum_particle_e_sec.(e,m,h=h,part_fin=part_fin)*g
        else
            Φ1 = spectrum_particle_e_sec(e,m,h=h,part_fin=part_fin)*g
        end
        
        if length(x) > 1
            i = searchsortedfirst(x,1e-6)
        elseif x > 1e-6
            i = 1
        else
            i = 2
        end
        
        if i==1
            Φ = Φ1
        elseif i>= length(x)
            Φ = Φ0
        else
            Φ = vcat(Φ0[1:i-1],Φ1[i:end])
        end
    end

    if k == 0
        return Φ
    end

    if isa(m,Array) || isa(m,LinRange)
        mask = m .<= m_init*q
        if burst
            Φ[mask] .*= 1 ./ bh_entropy.(m[mask]).^k
        elseif length(m_init) > 1
            Φ[mask] .*= 1 ./ bh_entropy.(m_init[mask]*q).^k
        else
            Φ[mask] *= 1 / bh_entropy(m_init*q)^k
        end
        return Φ
    else
        if m > m_init*q
            return Φ
        elseif burst
            return Φ/bh_entropy(m)^k
        else
            return Φ/bh_entropy(m_init*q)^k
        end
    end
end

"""
    spectrum_integrate_e(e_min,e_max,m,s,g,m_init=m,k=0,burst=true,q=0.5)

Returns the integrated spectrum ∫ dE (dN/dEdt) in units of 1/s.
Assumes instantenous memory burden effect.

Input:

    e_min : Lower integral bound in GeV
    e_max : Upper integral bound in GeV
    m : Black hole mass in g
    s : spin
    g : degrees of freedom
    m_init : Initial black hole mass in g (optional)
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function spectrum_integrate_e(e_min,e_max,m,s,g,a=0;m_init=m,k=0,burst=true,q=0.5,sec=false,h="select")
    integrand(e,p) = spectrum_particle_e(e,m,s,g,a,m_init=m,k=0,burst=true,q=0.5,sec=false,h="select")
    prob = IntegralProblem(integrand, BigFloat(e_min), BigFloat(e_max))
    sol = solve(prob, QuadGKJL())
    return sol.u
end

"""
    spectrum_integrate_fast(e_min,e_max,m,s,g,m_init=m,k=0,burst=true,q=0.5)

Returns the integrated spectrum ∫ dE (dN/dEdt) in units of 1/s.
Assumes instantenous memory burden effect.

Input:

    e_space : Energie samples for integration in GeV
    m : Black hole mass in g
    s : spin
    g : degrees of freedom
    m_init : Initial black hole mass in g (optional)
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function spectrum_integrate_fast(e_space,m,s,g,a=0;m_init=m,k=0,burst=true,q=0.5,sec=false,h="select")
    Φ = spectrum_particle_e(e_space,m,s,g,a,m_init=m_init,k=k,burst=burst,q=q,sec=sec,h=h)
    return trapz(e_space,Φ)
end

"""
    generate_spectrum_interpolation(data)

Returns an array which can be used to interpolate spectra with the function spectrum_interpolate(int)
Input array with dimensions (# of masses, # of particles, # of energies).
"""
function generate_spectrum_interpolation(data,m_space;col=2)
    dM = m_space[2]/m_space[1]
    int = []
    e_space = data[1,1,:]
    e1 = searchsortedfirst(e_space,e_space[1]*dM)
    e_space_red = e_space[e1:end]
    for i=1:length(m_space)-1
        int1 = linear_interpolation(e_space,loginf.(data[i,col,:]),extrapolation_bc=NaN)
        int2 = linear_interpolation(e_space,loginf.(data[i+1,col,:]),extrapolation_bc=NaN)
        push!(int,linear_interpolation((log.([m_space[i],m_space[i+1]]),e_space_red),[int1(e_space_red) int2(e_space_red/dM)]',extrapolation_bc=NaN))
    end
    return int
end

function spectrum_interpolate(int,m,e;s=1.0,g=2.0,a=0.0)
    if m > m_space[1] && m < m_space[end]
        ind = searchsortedlast(m_space,m)
    elseif m==m_space[1]
        ind = 1
    elseif m==m_space[end]
        ind = length(int)
    elseif m < m_space[1] || m > m_space[end]
        return spectrum_particle_e(e,m,s,g,a)
    end
    
    if isa(e,Array) || isa(e,LinRange)
        Φ = exp.(int[ind]([log(m)],e*m/m_space[ind]))[:]
        for i=1:length(Φ)
            if Φ[i]==0.0 || isnan(Φ[i])
                Φ[i] = spectrum_particle_e(e[i],m,s,g,a)
            end
        end
        return Φ
    else
        Φ = exp(int[ind](log(m),e*m/m_space[ind]))
        if Φ == 0.0 || isnan(Φ)
            Φ = spectrum_particle_e(e,m,s,g,a)
        end
        return Φ
    end
end

function hdm_secondary(nij,e,x,e_sec,m,part_fin="Photon",a=0)
    sec = zeros(length(e_sec))
    for i=1:size(nij)[3]
        nij_int = linear_interpolation((e,x),nij[:,:,i],extrapolation_bc=0)
        for j=1:length(sec)
            e_prim = 2 * e_sec[j] ./ x
            prim = spectrum_particle_e(e_prim/2,m,s_part[i],g_part[i],a)
            sec[j] += abs(trapz(log.(e_prim),prim.*nij_int.(e_prim,x)))
        end
    end
    ind = findfirst(x->x=="Photon",part)
    sec .+= spectrum_particle_e(e_sec,m,s_part[ind],g_part[ind],a)
    
    return sec
end

function spectrum_redshift(e,m_init,s=1.0,g=2.0,a_bh=0;t_start=t_rec,t_end=t_0,a_start=nothing,a_end=nothing,k=0,burst=true,q=0.5,n_steps=1000,sec=false,h="select",a_rel=0.1,go=false,instant=false,split_flux=false,md=false)
    
    args = (m_init=m_init,k=k,burst=burst,q=q,sec=sec,h=h,go=go)
                            
    tail = false
    
    if instant
        return spectrum_redshift_instant(e,m_init,s,g,a_bh;t_start=t_start,t_end=t_end,a_start=a_start,a_end=a_end,n_steps=n_steps,a_rel=a_rel,k=k,burst=burst,q=q,sec=sec,h=h,go=go,split_flux=split_flux,md=md)
    end
    
    m_start, a_bh_start = bh_mass_t(m_init,t_start,a_bh,k=k,burst=burst,q=q,return_spin=true)
    m_fin, a_bh_fin = bh_mass_t(m_init,t_end,a_bh,k=k,burst=burst,q=q,return_spin=true)

    if m_fin < m_init*1e-5
        t_fin = bh_lifetime_estimate(m_init,a_bh,k=k,burst=burst,q=q)
        a_fin = t_to_a_fast(t_fin)
        if (!sec && burst) || (!sec && k == 0)
            tail = true
            e_tail = 1e17/(m_init)/(t_to_a_fast(t_end)/a_fin)
            m_fin = 1e13/(e_tail*t_to_a_fast(t_end)/a_fin)
        else
            m_fin = minimum([1e10/(maximum(e)*t_to_a_fast(t_end)/a_fin),m_init*1e-5])
        end
    else
        t_fin = t_end
    end
    
    if isnothing(a_start)
        a_start = t_to_a_fast(t_start)
    end
    
    if isnothing(a_end)
        a_end = t_to_a_fast(t_end)
        a_fin = t_to_a_fast(t_fin)
        
    end
    
    if m_start == 0.0
        if isa(e,Array)
            return zeros(length(e))
        else
            return 0.0
        end
    end

    if k > 0 && m_fin <= m_init*q
        if m_start > m_init*q
            sc_phase = true
            burden_phase = true
            m_half = m_init*q
            t_half = bh_lifetime_estimate(m_init,a_bh,m_fin=m_init*q)
            if a_bh > 0.0
                a_bh_half = bh_mass_t(m_init,t_half,a_bh,return_spin=true)[2]
            else
                a_bh_half = 0.0
            end
            a_half = t_to_a_fast(t_half)
        else
            m_half = m_start
            a_bh_half = a_bh_start
            a_half = a_start
            t_half = t_start
            sc_phase = false
            burden_phase = true
        end
    else
        m_half = m_fin
        a_bh_half = a_bh_fin
        t_half = t_fin
        sc_phase = true
        burden_phase = false
    end
    
    Φ1 = zeros(length(e))
    Φ2 = zeros(length(e))

    if sc_phase

        if abs(m_start-m_half)/(m_start) < 0.1
            t_space = 10 .^ LinRange(log10(t_start),log10(t_half),n_steps)
            if !md
                a_space = t_to_a_fast.(t_space)
            else
                a_space = (t_space/t_end).^(2/3)*a_fin
            end
            M_space, a_bh_space = bh_mass_t_vec(m_init,t_space,a_bh,return_spin=true)

            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i] * a_end ./ a_space,M_space,s,g,a_bh_space,sec=sec,h=h,go=go)
                y = spec*a_end./a_space
                Φ1[i] += abs(trapz(t_space,y))
            end

        else

            if tail
                i_tail = minimum([maximum([2,searchsortedfirst(e,e_tail)]),length(e)])
            else
                i_tail = length(e)
            end

            t_crit = bh_lifetime_estimate(m_init,a_bh,m_fin=m_start-0.1*(m_start-m_half))
            t_space = 10 .^ LinRange(log10(t_start),log10(t_crit),n_steps)
            if !md
                a_space = t_to_a_fast.(t_space)
            else
                a_space = (t_space/t_end).^(2/3)*a_fin
            end

            if a_bh > 0.0
                M_space, a_bh_space = bh_mass_t_vec(m_init,t_space,a_bh,return_spin=true)
            else
                M_space = bh_mass_t_vec(m_init,t_space)
                a_bh_space = zeros(n_steps)
            end

            Threads.@threads for i=1:i_tail
                spec = spectrum_particle_e(e[i] * a_end ./ a_space,M_space,s,g,a_bh_space,sec=sec,h=h,go=go)
                y = spec*a_end./a_space 
                Φ1[i] += abs(trapz(t_space,y))
            end
            
            M_space = 10 .^ LinRange(log10(m_half),log10(M_space[end]),n_steps)
            t_space = broadcast(x->bh_lifetime_estimate(m_init,a_bh,m_fin=x),M_space)
            if !md
                a_space = t_to_a_fast.(t_space)
            else
                a_space = (t_space/t_end).^(2/3)*a_fin
            end

            Threads.@threads for i=1:i_tail
                spec = spectrum_particle_e(e[i] * a_end ./ a_space,M_space,s,g,a_bh_space,sec=sec,h=h,go=go)./mass_loss_rate.(M_space,a_bh_space)
                y = spec*a_end./a_space 
                Φ1[i] += abs(trapz(M_space,y))
            end

            Threads.@threads for i=i_tail+1:length(e)
                Φ1[i] = Φ1[i_tail]*(e[i]/e[i_tail])^(-3)
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
                        
function spectrum_redshift_instant(e,m_init,s=1.0,g=2.0,a_bh=0;t_start=t_rec,t_end=t_0,a_start=nothing,a_end=nothing,k=0,burst=true,q=0.5,n_steps=1000,sec=false,h="select",a_rel=0.01,go=false,split_flux=false,md=false)
    
    args = (m_init=m_init,k=k,burst=burst,q=q,sec=sec,h=h,go=go)
    
    m_start, a_bh_start = bh_mass_t(m_init,t_start,a_bh,k=k,burst=burst,q=q,return_spin=true)
    m_fin, a_bh_fin = bh_mass_t(m_init,t_end,a_bh,k=k,burst=burst,q=q,return_spin=true)
    
    if m_fin < 1.0
        m_fin = 1.0
        t_fin = bh_lifetime_estimate(m_init,a_bh,k=k,burst=burst,q=q)
        a_fin = t_to_a_fast(t_fin)
    else
        t_fin = t_end
    end
    
    if isnothing(a_start)
        a_start = t_to_a_fast(t_start)
    end
    
    if isnothing(a_end)
        a_end = t_to_a_fast(t_end)
        a_fin = t_to_a_fast(t_fin)
        
    end
    
    if m_start == 0.0
        if isa(e,Array)
            return zeros(length(e))
        else
            return 0.0
        end
    end
    
    if k > 0 && m_fin <= m_init*q
        if m_start > m_init*q
            sc_phase = true
            burden_phase = true
            m_half = m_init*q
            t_half = bh_lifetime_estimate(m_init,a_bh,m_fin=m_init*q)
            if a_bh > 0.0
                a_bh_half = bh_mass_t(m_init,t_half,a_bh,return_spin=true)[2]
            else
                a_bh_half = 0.0
            end
            a_half = t_to_a_fast(t_half)
        else
            m_half = m_start
            a_bh_half = a_bh_start
            a_half = a_start
            t_half = t_start
            sc_phase = false
            burden_phase = true
        end
    else
        m_half = m_fin
        a_bh_half = a_bh_fin
        t_half = t_fin
        sc_phase = true
        burden_phase = false
    end
    
    Φ1 = zeros(length(e))
    Φ2 = zeros(length(e))

    if sc_phase

        if abs(m_start-m_half)/(m_start) < 1e-1
            t_space = 10 .^ LinRange(log10(t_start),log10(t_half),n_steps)
            M_space, a_bh_space = bh_mass_t_vec(m_init,t_space,a_bh,return_spin=true)
        
            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i],M_space,s,g,a_bh_space,sec=sec,h=h,go=go)
                Φ1[i] += abs(trapz(t_space,spec))
            end

        else
            M_space = 10 .^ LinRange(log10(m_half),log10(Float64(m_start)),n_steps)
            M_space[1] = m_half
            M_space[end] = m_start
            if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
                println("Error")
            end
            
            t_space = broadcast(x->bh_lifetime_estimate(m_init,a_bh,m_fin=x),M_space)
            if a_bh > 0.0
                M_space, a_bh_space = bh_mass_t_vec(m_init,reverse(t_space),a_bh,return_spin=true)
                M_space = reverse(M_space)
                a_bh_space = reverse(a_bh_space)
            else
                a_bh_space = zeros(n_steps)
            end

            a_space = t_to_a_fast.(t_space)
            ind = searchsortedfirst(a_space[1:end-1]./a_space[2:end] .- 1 , a_rel)
            
            if ind <= n_steps
                n_ext = Int(ceil(log(a_space[ind]/a_space[end])/log(1 + a_rel)+1))
                a_ext = 10 .^ LinRange(log10(a_space[ind]),log10(a_space[end]),n_ext)
                a_space = vcat(a_space[1:ind-1],a_ext)
                M_ext, a_bh_ext = bh_mass_t_vec(m_init,reverse(a_to_t_fast.(a_ext)),a_bh,return_spin=true)
                M_space = vcat(M_space[1:ind-1],reverse(M_ext))
                a_bh_space = vcat(a_bh_space[1:ind-1],reverse(a_bh_ext))
            end

            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i],M_space,s,g,a_bh_space,sec=sec,h=h,go=go)./mass_loss_rate.(M_space,a_bh_space)
                Φ1[i] += abs(trapz(M_space,spec))
            end
                
        end
            
    end

    if burden_phase
        if abs(m_half-m_fin)/(m_half) < 1e-3

            a_space = 10 .^ LinRange(log10(a_half),log10(a_end),n_steps)
            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i],m_half,s,g,a_bh_half;args...)
                H = H_a.(a_space)
                y = spec./H./a_space
                Φ2[i] += abs(trapz(a_space,y))
            end
            
        elseif abs(m_half-m_fin)/(m_half) < 1e-1

            a_space = 10 .^ LinRange(log10(a_half),log10(a_end),n_steps)
            t_space = a_to_t_fast.(a_space)
            M_space, a_bh_space = bh_mass_t_vec(m_init,t_space,a_bh,m_init=m_init,k=k,q=q,burst=burst,return_spin=true)

            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i],M_space,s,g,a_bh_space;args...)
                H = H_a.(a_space)
                y = spec./H./a_space
                Φ2[i] += abs(trapz(a_space,y))
            end

        else
                        
            M_space = 10 .^ LinRange(log10(m_fin),log10(m_half),n_steps)
            M_space[1] = m_fin
            M_space[end] = m_half
            if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
                println("Error")
            end
            
            t_space = broadcast(x->bh_lifetime_estimate(m_init,a_bh,m_fin=x,k=k,burst=burst,q=q),M_space)
            if a_bh_half > 0.0
                M_space, a_bh_space = bh_mass_t_vec(m_init,reverse(t_space),a_bh,k=k,burst=burst,q=q,return_spin=true)
                M_space = reverse(M_space)
                a_bh_space = reverse(a_bh_space)
                M_space[1] = m_fin
                M_space[end] = m_half
                if M_space[2] < M_space[1] || M_space[end-1] > M_space[end]
                    println("Error")
                end
            else
                a_bh_space = zeros(n_steps)
            end
            
                            
                                
            a_space = t_to_a_fast.(t_space)
            ind = searchsortedfirst(a_space[1:end-1]./a_space[2:end] .- 1 , a_rel)
            
            if ind <= n_steps
                n_ext = Int(ceil(log(a_space[ind]/a_space[end])/log(1 + q)+1))
                a_ext = 10 .^ LinRange(log10(a_space[ind]),log10(a_space[end]),n_ext)
                a_space = vcat(a_space[1:ind-1],a_ext)
                t_ext = reverse(a_to_t_fast.(a_ext))
                t_ext[1] = t_half
                M_ext, a_ext = bh_mass_t_vec(m_init,t_ext,a_bh,k=k,burst=burst,q=q,return_spin=true)
                M_space = vcat(M_space[1:ind-1],reverse(M_ext))
                a_bh_space = vcat(a_bh_space[1:ind-1],reverse(a_ext))
            end
            
            Threads.@threads for i=1:length(e)
                spec = spectrum_particle_e(e[i],M_space,s,g,a_bh_space;args...)./mass_loss_rate.(M_space,a_bh_space,k=k,burst=burst,q=q,m_init=m_init)
                Φ2[i] += abs(trapz(M_space,spec))
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

#function spectrum_redshift_approximation(x,t_start,t_end,a=1.39077e-5,b=-0.0873,c=0.05660,d=0.71195)
#    x0 = 3.7*sqrt(t_start/t_end)
#    if x > x0
#        return  a*(x + b*x^2 + c*x^3) / (exp(d*x) + 1)
#    else
#        return (x/x0)^1.93 * 8e-6 * x0
#    end
#end

function spectrum_redshift_approximation(x,t_start,t_end,a,b,c,d,e)
    x0 = 0*sqrt(t_start/t_end)
    if x > x0
        return a*(x + b*x^2 + c*x^3) / (exp(d*x) + e)
    else
        return (x/x0)^1.93 * 1.808e-6 * x0
    end
end

function spectrum_redshift_approximation_burst(x,t_start,t_end,a,b,c,d,e)
    x0 = 0*sqrt(t_start/t_end)
    if x > x0
        return a*(x + b*x^2 + c*x^3) / (d + e*x^6)
    else
        return (x/x0)^1.93 * 1.808e-6 * x0
    end
end


"""
    burden_factor(k)

Returns the dimensionless pre-factor of the memory burden effect

Input:

    k : Power law index
"""
function burden_factor(k)
    if k == 0
        return 1
    else
        return 2.653e10^k
    end
end

"""
    mass_loss_rate(m,m_init=m,k=0,burst=true,q=0.5)

Returns the mass-loss rate of a black hole with instantaneous memory burden effect

Input:

    m : Black hole mass in g
    a : Dimensionless black hole spin parameter
    m_init : Initial black hole mass in g (optional)
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function mass_loss_rate(m,a=0;m_init=m,k=0,burst=true,q=0.5)
    if m>m_init*q || k == 0.0
        return fM_int(m,a)/m^2
    elseif burst
        return fM_int(m,a)/m^(2+2*k)/burden_factor(k)
    else
        return fM_int(m_init*q,a)/(m_init*q)^(2+2*k)/burden_factor(k)
    end
end

"""
    bh_lifetime_integrate(m_init,m_fin=0,k=0,burst=true,q=0.5)

Returns the exact timescale for the black hole mass to decrease from m_init to m_fin.
Assumes instantenous memory burden effect. Only valid for non-spinning black holes.

Input:

    m_init : Initial black hole mass in g
    m_fin : Final black hole mass in g (optional)
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function bh_lifetime_integrate(m_init;m_fin=0,k=0,burst=true,q=0.5)
    integrand(x,p) = 1/mass_loss_rate(x,m_init=m_init,k=k,burst=burst,q=q)
    prob = IntegralProblem(integrand, BigFloat(m_fin), BigFloat(m_init))
    sol = solve(prob, QuadGKJL())
    return sol.u
end


"""
    bh_lifetime_fast(m_init,m_fin=0,k=0,burst=true,q=0.5)

Quick estimate of the timescale for the black hole mass to decrease from m_init to m_fin.
Assumes f(M) is constant over the evolution. Assumes instantenous memory burden effect.
Only valid for non-spinning black holes.

Input:

    m_init : Initial black hole mass in g
    m_fin : Final black hole mass in g (optional)
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function bh_lifetime_fast(m_init;m_fin=0,k=0,burst=true,q=0.5)
    if k==0 || m_fin>m_init*q
        return (m_init^3-m_fin^3)/3/fM_int(m_init)
    elseif burst
        return (m_init^3-(m_init*q)^3)/3/fM_int(m_init)+((m_init*q)^(3+2*k)-m_fin^(3+2*k))/(3+2*k)/fM_int(m_init*q)*burden_factor(k) 
    else
        return (m_init^3-(m_init*q)^3)/3/fM_int(m_init)+((m_init*q)-m_fin)*(m_init*q)^(2+2*k)*burden_factor(k)/fM_int(m_init*q)
    end
end

"""
    bh_lifetime_estimate(m_init;m_fin=0,k=0,burst=true,q=0.5)

Estimates the timescale for the black hole mass to decrease from m_init to m_fin.
Assumes instantenous memory burden effect. Only valid for non-spinning black holes.

Input:

    m_init : Initial black hole mass in g
    m_fin : Final black hole mass in g (optional)
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function bh_lifetime_estimate(m;m_init=m,m_fin=0,k=0,burst=true,q=0.5)
    args = (m_init=m_init,k=k,burst=burst,q=q)
    τ = 0
    if m > m_tab[1]
        i1 = searchsortedlast(m_tab,m)
        i2 = searchsortedfirst(m_tab,m_fin)
        if i1>i2
            τ += t_step_m(m,m_tab[i1],fM_int(m);args...)
            for j=i1:-1:i2+1
                τ += t_step_m(m_tab[j],m_tab[j-1],fM_tab[j];args...)
            end
            if i2>1
                τ += t_step_m(m_tab[i2],m_fin,fM_tab[i2];args...)
                return τ
            end
        else
            return t_step_m(m,m_fin,fM_int(m);args...)
        end
        m = m_tab[1]
    end
    τ += t_step_m(m,m_fin,fM_int(m);args...)
    return τ
end
                
"""
    bh_lifetime_estimate(m_init;m_fin=0,k=0,burst=true,q=0.5)

Estimates the timescale for the black hole mass to decrease from m_init to m_fin.
Assumes instantenous memory burden effect.

Input:

    m_init : Initial black hole mass in g
    a : Dimensionless black hole spin parameter
    m_fin : Final black hole mass in g (optional)
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function bh_lifetime_estimate(m,a;m_init=m,m_fin=0,k=0,burst=true,q=0.5)
    if a <= a_tab[2]
        return bh_lifetime_estimate(m,m_init=m_init,m_fin=m_fin,k=k,burst=burst,q=q)
    end
    args = (m_init=m_init,k=k,q=q,burst=burst)
    τ = 0
    i = searchsortedlast(a_tab,a)
    dt = t_step_a(a,a_tab[i],m,fM_int(m,a),gM_int(m,a);args...)
    m_temp = M_step(m,dt,fM_int(m,a);args...)
    if m_temp > m_fin
        m = m_temp
        τ += dt
        while i>2
            dt = t_step_a(a_tab[i],a_tab[i-1],m,fM_int(m,a_tab[i]),gM_int(m,a_tab[i]);args...)
            m_temp = M_step(m,dt,fM_int(m,a_tab[i]);args...)
            if m_temp > m_fin
                m = m_temp
                τ += dt
            else
                return τ + t_step_m(m,m_fin,fM_int(m,a_tab[i]);args...)
            end
            i -= 1
        end
        if i == 2
            return τ + bh_lifetime_estimate(m;m_fin=m_fin,args...)
        end
    else
        return t_step_m(m,m_fin,fM_int(m,a);args...)
    end
end
    

"""
    bh_mass_t(m,t;k=0,burst=true,q=0.5,dtype=Float64)

Computes the mass of a non-spinning black hole with initial mass M after a time t.

Input:

    m : Initial black hole mass in g
    t : Elapsed time in seconds
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function bh_mass_t(m,t;m_init=m,k=0,burst=true,q=0.5,dtype=Float64)
    m_init = dtype(m_init)
    m = dtype(m)
    t = dtype(t)
    if m > m_tab[1]
        i = searchsortedlast(m_tab,m)
        dt0 = t_step_m(m,m_tab[i],fM_int(m),m_init=m_init,k=k,burst=burst,q=q)
        if t > dt0
            m = dtype(m_tab[i])
            t_remain = t - dt0
            if i==1
                return M_step(m,t_remain,fM_tab[1],m_init=m_init,k=k,burst=burst,q=q)
            end
            dt = t_step_m(m,m_tab[i-1],fM_tab[i],m_init=m_init,k=k,burst=burst,q=q)
            while t_remain > dt
                t_remain = t_remain - dt
                i -= 1
                if i>1
                    m = dtype(m_tab[i])
                    dt = t_step_m(m,m_tab[i-1],fM_tab[i],m_init=m_init,k=k,burst=burst,q=q)
                else
                    return M_step(m_tab[1],t_remain,fM_tab[1],m_init=m_init,k=k,burst=burst,q=q)
                end
            end
            return M_step(m,t_remain,fM_tab[i],m_init=m_init,k=k,burst=burst,q=q)
        else
            return M_step(m,t,fM_int(m),m_init=m_init,k=k,burst=burst,q=q)
        end    
    else
        return M_step(m,t,fM_tab[1],m_init=m_init,k=k,burst=burst,q=q)
    end
end

"""
    bh_mass_t(m,t,a;k=0,burst=true,q=0.5,dtype=Float64)

Computes the mass of a black hole with initial mass M and spin a after a time t.

Input:

    m : Initial black hole mass in g
    t : Elapsed time in seconds
    a : Dimensionless black hole spin parameter
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
    spin : If true, the value of the spin is returned as well (optional)
"""
function bh_mass_t(m,t,a;m_init=m,k=0,burst=true,q=0.5,dtype=Float64,return_spin=false)
    if a <= a_tab[2]
        if !return_spin
            return bh_mass_t(m,t,k=k,burst=burst,q=q,dtype=dtype)
        else
            return bh_mass_t(m,t,k=k,burst=burst,q=q,dtype=dtype), 0.0
        end
    end
    m_init = dtype(m_init)
    m = dtype(m)
    t = dtype(t)
    a = dtype(a)
    t_elapsed = dtype(0)
    i = searchsortedlast(a_tab,a)
    dt0 = t_step_a(a,a_tab[i],m,fM_int(m,a),gM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
    if t > dt0
        m = M_step(m,dt0,fM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
        a = dtype(a_tab[i])
        t_remain = t - dt0
        if i==2
            if !return_spin
                return bh_mass_t(m,t_remain,m_init=m_init,k=k,burst=burst,q=q,dtype=dtype)
            else
                return bh_mass_t(m,t_remain,m_init=m_init,k=k,burst=burst,q=q,dtype=dtype), 0.0
            end
        end
        dt = t_step_a(a,a_tab[i-1],m,fM_int(m,a),gM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
        while t_remain > dt
            m = M_step(m,dt,fM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
            t_remain = t_remain - dt
            i -= 1
            if i>2
                a = dtype(a_tab[i])
                dt = t_step_a(a,a_tab[i-1],m,fM_int(m,a),gM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
            else
                if !return_spin
                    return bh_mass_t(m,t_remain,m_init=m_init,k=k,burst=burst,q=q,dtype=dtype)
                else
                    return bh_mass_t(m,t_remain,m_init=m_init,k=k,burst=burst,q=q,dtype=dtype), 0.0
                end
            end
        end
        if !return_spin
            return M_step(m,t_remain,fM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
        else
            return M_a_step(m,a,t_remain,fM_int(m,a),gM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
        end
    else
        if !return_spin
            return M_step(m,t,fM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
        else
            return M_a_step(m,a,t,fM_int(m,a),gM_int(m,a),m_init=m_init,k=k,burst=burst,q=q)
        end
    end    
end

"""
    bh_mass_t_vec(m,t,k=0,burst=true,q=0.5)

Computes the mass of a non-spinning black hole with initial mass m at different points in time t.

Input:

    m : Initial black hole mass in g
    t : Sorted array of times in seconds
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function bh_mass_t_vec(m,t;m_init=m,k=0,burst=true,q=0.5,dtype=Float64)
    args = (m_init=m_init,k=k,burst=burst,q=q)
    mt = zeros(dtype,length(t))
    tt = 0.0
    if m > m_tab[1]
        i = searchsortedlast(m_tab,m)
        tt += t_step_m(m,m_tab[i],fM_int(m);args...)
        i1 = searchsortedlast(t,tt)
        for j = 1:i1
            mt[j] = M_step(m,t[j],fM_int(m);args...)
        end
        while i>1 && mt[end] == 0.0
            dt = t_step_m(m_tab[i],m_tab[i-1],fM_tab[i];args...)
            i_temp = searchsortedlast(t,tt+dt)
            if i_temp > i1
                for j=i1+1:i_temp
                    mt[j] = M_step(m_tab[i],t[j]-tt,fM_tab[i];args...)
                end
                i1 = i_temp
            end
            tt += dt
            i -= 1
        end
        if mt[end] == 0.0 && i==1
            for j=i1+1:length(mt)
                mt[j] = M_step(m_tab[1],t[j]-tt,fM_tab[1];args...)
            end
        end    
    else
        for j=1:length(mt)
            mt[j] = M_step(m,t[j],fM_tab[1];args...)
        end
    end
    return mt
end

"""
    bh_mass_t_vec(m,t,a,k=0,burst=true,q=0.5)

Computes the mass of a black hole with initial mass m at different points in time t.

Input:

    m : Initial black hole mass in g
    t : Sorted array of times in seconds
    k : Power law index k (optional)
    burst : If true then T~1/m during memory burden (optional)
    q : Fraction of mass loss after which memory burden starts (optional)
"""
function bh_mass_t_vec(m,t,a;m_init=m,k=0,burst=true,q=0.5,dtype=Float64,return_spin=false)
    args = (m_init=m_init,k=k,burst=burst,q=q)
    if a <= a_tab[2]
        if !return_spin
            return bh_mass_t_vec(m,t,m_init=m_init,k=k,burst=burst,q=q,dtype=Float64)
        else
            return bh_mass_t_vec(m,t,m_init=m_init,k=k,burst=burst,q=q,dtype=Float64), zeros(length(t))
        end
    end
    mt = zeros(dtype,length(t))
    if return_spin
        at = zeros(dtype,length(t))
    end
    tt = 0.0
    i = searchsortedlast(a_tab,a)
    tt += t_step_a(a,a_tab[i],m,fM_int(m,a),gM_int(m,a);args...)
    i1 = searchsortedlast(t,tt)
    Threads.@threads for j = 1:i1
        if !return_spin
            mt[j] = M_step(m,t[j],fM_int(m,a);args...)
        else
            mt[j], at[j] = M_a_step(m,a,t[j],fM_int(m,a),gM_int(m,a);args...)
        end
    end
    m = M_step(m,tt,fM_int(m,a);args...)
    while i>2 && i1 < length(t)
        dt = t_step_a(a_tab[i],a_tab[i-1],m,fM_int(m,a_tab[i]),gM_int(m,a_tab[i]);args...)
        i_temp = searchsortedlast(t,tt+dt)
        if i_temp > i1
            for j=i1+1:i_temp
                if !return_spin
                    mt[j] = M_step(m,t[j]-tt,fM_int(m,a_tab[i]);args...)
                else
                    mt[j], at[j] = M_a_step(m,a_tab[i],t[j]-tt,fM_int(m,a_tab[i]),gM_int(m,a_tab[i]);args...)
                end
            end
            i1 = i_temp
        end
        m = M_step(m,dt,fM_int(m,a_tab[i]);args...)
        tt += dt
        i -= 1
    end
    if i1 < length(t)
        mt[i1+1:end] = bh_mass_t_vec(m,t[i1+1:end] .- tt,m_init=m_init,k=k,burst=burst,q=q,dtype=dtype)
        if return_spin
            at[i1+1:end] .= 0.0
        end
    end
    if !return_spin
        return mt
    else
        return mt, at
    end
end


function M_step(m,dt,f;m_init=m,k=0,burst=true,q=0.5)
    if k==0.0
        return relu(m^3-3*f*dt)^(1/3)
    end
    if m <= m_init*q
        if burst
            λ = 3 + 2*k
            a = m
            b = λ*f*dt/burden_factor(k)
            if a^λ < 1e10*b
                return relu(a^λ-b)^(1/λ)
            else
                return a-b/λ/a^(λ-1)
            end
        else
            return relu(m-dt*mass_loss_rate(m,m_init=m_init,k=k,burst=burst,q=q))
        end
    else
        m_temp = relu(m^3-3*f*dt)^(1/3)
        if m_temp > m_init*q
            return m_temp
        else
            if burst
                λ = 3 + 2*k
                dt1 = (m^3-(m_init*q)^3)/(3*f)
                dt2 = dt-dt1
                a = m_init*q
                b = λ*f*dt2/burden_factor(k)
                if a^λ < 1e10*b
                    return relu(a^λ-b)^(1/λ)
                else
                    return a-b/λ/a^(λ-1) #Taylor approximation to avoid catastrophic cancellation
                end
            else
                dt1 = (m^3-(m_init*q)^3)/(3*f)
                dt2 = dt-dt1
                return relu(m_init*q-dt2*mass_loss_rate(m_init*q,m_init=m_init,k=k,burst=burst,q=q))
            end
        end
    end
end
        
function M_a_step(m,a,dt,f,g;m_init=m,k=0,burst=true,q=0.5)
    if k==0.0
        return relu(m^3-3*f*dt)^(1/3), a*exp((2*f-g)/m^3*dt)
    end
    if m <= m_init*q
        if burst
            λ = 3 + 2*k
            α = m
            β = λ*f*dt/burden_factor(k)
            if α^λ < 1e10*β
                return relu(α^λ-β)^(1/λ), a*exp((2*f-g)/burden_factor(k)/m^λ*dt)
            else
                return α-β/λ/α^(λ-1), a*exp((2*f-g)/burden_factor(k)/m^λ*dt)
            end
        else
            return relu(m-dt*mass_loss_rate(m,m_init=m_init,k=k,burst=burst,q=q)), a*exp(dt/m/burden_factor(k)*(2*fM_int(m_init*q,a)/(q*m_init)^(2+2*k)-gM_int(m_init*q,a)/m^(2+2*k)))
        end
    else
        m_temp = relu(m^3-3*f*dt)^(1/3)
        if m_temp > m_init*q
            return m_temp, a*exp((2*f-g)/m^3*dt)
        else
            if burst
                λ = 3 + 2*k
                dt1 = (m^3-(m_init*q)^3)/(3*f)
                a_mid = a*exp((2*f-g)/m^3*dt1)
                dt2 = dt-dt1
                α = m_init*q
                β = λ*f*dt2/burden_factor(k)
                if α^λ < 1e10*β
                    return relu(α^λ-β)^(1/λ), a_mid*exp((2*f-g)/burden_factor(k)/m^λ*dt2)
                else
                    return α-β/λ/α^(λ-1), a_mid*exp((2*f-g)/burden_factor(k)/m^λ*dt2) #Taylor approximation to avoid catastrophic cancellation
                end
            else
                dt1 = (m^3-(m_init*q)^3)/(3*f)
                a_mid = a*exp((2*f-g)/m^3*dt1)
                dt2 = dt-dt1
                return relu(m_init*q-dt2*mass_loss_rate(m_init*q,m_init=m_init,k=k,burst=burst,q=q)), a_mid*exp(dt2/m/burden_factor(k)*(2*fM_int(m_init*q,a)/(q*m_init)^(2+2*k)-gM_int(m_init*q,a)/m^(2+2*k)))
            end
        end
    end
end

function t_step_m(m1,m2,f;m_init=m1,k=0,burst=true,q=0.5)
    if k==0.0
        return (m1^3-m2^3)/(3*f)
    end
    if m1<=m_init*q
        if burst
            λ = 3 + 2*k
            return (m1^λ-m2^λ)/(f*λ)*burden_factor(k)
        else
            return (m1-m2)/mass_loss_rate(m1,m_init=m_init,k=k,burst=burst,q=q)
        end
    elseif m2>m_init*q
        return (m1^3-m2^3)/(3*f)
    else
        if burst
            λ = 3 + 2*k
            a = (m_init*q)
            b = m2
            if abs(a/b-1) > 1e-5
                return (m1^3-(m_init*q)^3)/(3*f)+(a^λ-b^λ)/(f*λ)*burden_factor(k)
            else
                return (m1^3-(m_init*q)^3)/(3*f)+λ*a^(λ-1)*(a-b)/(f*λ)*burden_factor(k) #Taylor approximation to avoid catastrophic cancellation
            end
        else
            return (m1^3-(m_init*q)^3)/(3*f)+(m_init*q-m2)/mass_loss_rate(m2,m_init=m_init,k=k,burst=burst,q=q)
        end
    end
end
        
function t_step_a(a1,a2,m,f,g;m_init=m,k=0,burst=true,q=0.5)
    if k==0.0
        return log(a2/a1)*m^3/(2*f-g)
    end
    if m<=m_init*q
        if burst
            λ = 3 + 2*k
            return log(a2/a1)*m^λ/(2*f-g)*burden_factor(k)
        else
            return log(a2/a1)*m*burden_factor(k)/(2*fM_int(m_init*q,a1)/(q*m_init)^(2+2*k)-gM_int(m_init*q,a1)/m^(2+2*k))
        end
    else
        dt = log(a2/a1)*m^3/(2*f-g)
        m_temp = relu(m^3-3*f*dt)^(1/3)
        if m_temp > m_init*q
            return log(a2/a1)*m^3/(2*f-g)
        else
            if burst
                λ = 3 + 2*k
                dt = (m^3-(m_init*q)^3)/(3*f)
                a_mid = a1*exp((2*f-g)/m^3*dt)
                return dt + log(a2/a_mid)*(m_init*q)^λ/(2*f-g)*burden_factor(k)
            else
                dt = (m^3-(m_init*q)^3)/(3*f)
                a_mid = a1*exp((2*f-g)/m^3*dt)
                return dt + log(a2/a_mid)*m*burden_factor(k)/(2*fM_int(m_init*q,a_mid)/(q*m_init)^(2+2*k)-gM_int(m_init*q,a_mid)/m^(2+2*k))
            end
        end
    end         
end

function m_init_estimate(m_fin,t;k=0,burst=true,q=0.5)
    m = (m_fin^3+3*fM_int(Inf)*t)^(1/3)
    if k==0
        return m
    elseif m*q < m_fin
        return m
    elseif burst
        λ = 3 + 2*k
        m = (m_fin^λ+λ*fM_int(0)*t/burden_factor(k))^(1/λ)
        m_space = LinRange(m,m/q,200)
        res = abs.(t .- bh_lifetime_fast.(m_space,m_fin=m_fin,k=k,burst=burst,q=q))
        return m_space[argmin(res)]
    else
        m_space = LinRange(m_fin,m,200)
        res = abs.(t .- bh_lifetime_fast.(m_space,m_fin=m_fin,k=k,burst=burst,q=q))
        return m_space[argmin(res)]
    end 
end
            
function get_m_init(m_fin,t;k=0,burst=true,q=0.5,dtype=Float64,maxiters=100)
            
    if m_fin == 0
        if k == 0
            m_min = (3*fM_int(Inf)*t)^(1/3)
            m_max = (3*fM_int(0)*t)^(1/3)
        else
            m_min = (fM_int(Inf)*t/burden_factor(k))^(1/(3+2*k))
            m_max = (3*fM_int(0)*t)^(1/3)
        end
        #println(m_min," ",m_max," ",k," ",t)
        g(x) = t - bh_lifetime_estimate(x,k=k,burst=burst,q=q)
        return find_zero(g,(m_min*0.99,m_max*1.01),maxiters=maxiters)
    elseif k == 0
        dm = 3*fM_int(Inf)*t
        if dm/m_fin^3 > 1e-14
            m_min = (m_fin^3 + 3*fM_int(Inf)*t)^(1/3)
            m_max = (m_fin^3 + 3*fM_int(0)*t)^(1/3)
            if m_max < 1e9
                return m_max
            end
        else
            return m_fin
        end            
    else
        m = bh_mass_t(m_fin/q,t,k=k,burst=burst,q=q)
        if m > m_fin
            m_min = m_fin
            m_max = m_fin/q
        else
            m_min = m_fin/q
            m_max = (m_fin^3 + 3*fM_int(0)*t)^(1/3)
        end
    end
    f(x) = dtype(m_fin) - bh_mass_t(x,t,k=k,burst=burst,q=q,dtype=dtype)
    return find_zero(f,(m_min*0.99,m_max*1.01),maxiters=maxiters)
end        


"""
    bh_mass_t_smooth(m,t,k,p)

Computes the mass of a black hole with initial mass m after a time t.
Takes a smooth function k(m,m_init,p...).

Input:

    m : Initial black hole mass in g
    t : Elapsed time in seconds
    k : Power law function k(m,m_init,p...)
    p : Parameters passed to k(m,m_init,p...)
"""
function bh_mass_t_smooth(m,t,k,p)
    m_init = BigFloat(m)
    m = BigFloat(m)
    t = BigFloat(t)
    if m > m_tab[2]
        i = searchsortedlast(m_tab,m)
        dt0 = t_step_smooth(m,m_tab[i],fM_int(m),k,p,m_init)
        if t > dt0
            m = BigFloat(m_tab[i])
            t_remain = t - dt0
            if i==1
                return M_step_smooth(m,t,fM_tab[1],k,p,m_init)
            end
            dt = t_step_smooth(m_tab[i],m_tab[i-1],fM_tab[i],k,p,m_init)
            while t_remain > dt
                t_remain = t_remain - dt
                i -= 1
                if i>1
                    m = BigFloat(m_tab[i])
                    dt = t_step_smooth(m_tab[i],m_tab[i-1],fM_tab[i],k,p,m_init)
                else
                    return M_step_smooth(m,t_remain,fM_tab[1],k,p,m_init)
                end
            end
            return M_step_smooth(m,t_remain,fM_tab[i],k,p,m_init)
        else
            return M_step_smooth(m,t,fM_int(m_init),k,p,m_init)
        end    
    else
        return M_step_smooth(m,t,fM_tab[1],k,p,m_init)
    end
end

function k_poly(m,m_init,k,n)
    if m > m_init/2
        if n == Inf
            return 0
        else
            return (2-2*m/m_init)^n*k
        end
    else
        return k
    end
end


# NOT CORRECT
function M_step_smooth(m,dt,f,k,p,m_init)
    k_m = k(m,m_init,p...)
    λ = 3 + 2*k_m
    a = m
    b = λ*f*dt*burden_factor(k_m)
    if a^λ < 1e10*b
        return relu(a^λ-b)^(1/λ)
    else
        return a-b/λ/a^(λ-1) #Taylor approximation to avoid catastrophic cancellation
    end
end

# NOT CORRECT
function t_step_smooth(m1,m2,f,k,p,m_init)
    k_m = k(m1,m_init,p...)
    λ = 3 + 2*k_m
    return (m1^λ-m2^λ)/(f*λ*burden_factor(k_m))
end
        
function m_evap_k(k;burst=true,q=0.5,f=fM_tab[1])
    if burst
        return ((3+2*k)*t_0*f/burden_factor(k))^(1/(3+2*k))/q
    else
        return (t_0*f/burden_factor(k))^(1/(3+2*k))/q
    end
end

function k_evap_m(m_init;burst=true,q=0.5)
    g(x) = m_evap_k(x,burst=burst,q=q,f=fM_int(m_init)) - m_init
    if g(0) > 0
        return find_zero(g,(0,5))
    else
        return 0.0
    end
end
