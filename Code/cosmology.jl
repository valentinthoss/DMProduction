dir = @__DIR__

#Cosmological parameters
const H0_fid = 2.1927e-18 #Hubble constant in 1/s
const t_H_fid = 1/H0_fid #Hubble time in s
const Ω_m_fid = 0.3111
const Ω_Λ_fid = 0.6889
const Ω_r_fid = 9.4e-5

function H_z(z;Ω_m = Ω_m_fid,Ω_Λ = Ω_Λ_fid,Ω_r = Ω_r_fid,H0 = H0_fid)
    return H0*sqrt(Ω_m*(1+z)^3 + Ω_Λ + rad_correct(1/(1+z))*Ω_r*(1+z)^4)
end

function H_a(a;Ω_m = Ω_m_fid,Ω_Λ = Ω_Λ_fid,Ω_r = Ω_r_fid,H0 = H0_fid)
    return H0*sqrt(Ω_m*a^(-3) + Ω_Λ + rad_correct(a)*Ω_r*a^(-4))
end
    
function z_to_t(z;Ω_m = Ω_m_fid,Ω_Λ = Ω_Λ_fid,Ω_r = Ω_r_fid,H0 = H0_fid)
    return a_to_t(1/(1+z),Ω_m=Ω_m,Ω_Λ=Ω_Λ,Ω_r=Ω_r,H0=H0)
end

function a_to_t(a;Ω_m = Ω_m_fid,Ω_Λ = Ω_Λ_fid,Ω_r = Ω_r_fid,H0 = H0_fid)
    integrand(x,p) = 1/x/H_a(x;Ω_m=Ω_m,Ω_Λ=Ω_Λ,Ω_r=Ω_r,H0=H0)
    prob = IntegralProblem(integrand,Float64(0),Float64(a))
    sol = solve(prob, QuadGKJL())
    return sol.u
end

function t_to_a_estimate(t;Ω_m = Ω_m_fid,Ω_Λ = Ω_Λ_fid,Ω_r = Ω_r_fid,H0 = H0_fid)
    args = (Ω_m=Ω_m,Ω_Λ=Ω_Λ,Ω_r=Ω_r,H0=H0)
    if t<z_to_t(Ω_m/Ω_r-1;args...)
        a = sqrt(2*H0*sqrt(Ω_r)*t)
    elseif t<z_to_t((Ω_Λ/Ω_m)^(1/3)-1;args...)
        a = (3*H0*sqrt(Ω_m)*(t-z_to_t(Ω_m/Ω_r-1;args...))/2+(Ω_r/Ω_m)^(3/2))^(2/3)
    else
        a = (Ω_Λ/Ω_m)^(-1/3)*exp(H0*sqrt(Ω_Λ)*(t-z_to_t((Ω_Λ/Ω_m)^(1/3)-1;args...)))
    end
    return a
end

function t_to_a(t;Ω_m = Ω_m_fid,Ω_Λ = Ω_Λ_fid,Ω_r = Ω_r_fid,H0 = H0_fid,n_steps=5)
    args = (Ω_m=Ω_m,Ω_Λ=Ω_Λ,Ω_r=Ω_r,H0=H0)
    a0 = t_to_a_estimate(t;args...)
    a = a0
    for i=1:n_steps
        a += -a*H_a(a;args...)*(a_to_t(a;args...)-t)
    end
    return a
end

function t_to_z(t;Ω_m = Ω_m_fid,Ω_Λ = Ω_Λ_fid,Ω_r = Ω_r_fid,H0 = H0_fid,n_steps = 5)
    a =  t_to_a(t,Ω_m=Ω_m,Ω_Λ=Ω_Λ,Ω_r=Ω_r,H0=H0,n_steps=n_steps)
    if a > 1-1e-8
        return (z_to_t(0,Ω_m=Ω_m,Ω_Λ=Ω_Λ,Ω_r=Ω_r,H0=H0)-t)*H_z(0,Ω_m=Ω_m,Ω_Λ=Ω_Λ,Ω_r=Ω_r,H0=H0)
    else
        return 1/a - 1
    end
end

"""
    bh_mass_z(m,z,k=0)

Computes the mass of a black hole with initial mass m at redshift z.
Assumes the black hole is formed at infinite redshift.
Uses cosmological parameters from Planck 2018.

Input:

    m : Initial black hole mass in g
    t : Elapsed time in seconds
    k : Power law index k (optional)
"""
function bh_mass_z(m,z;k=0,burst=burst,q=q,Ω_m = Ω_m_fid,Ω_Λ = Ω_Λ_fid,Ω_r = Ω_r_fid,H0 = H0_fid)
    return bh_mass_t(m,z_to_t(z,Ω_m=Ω_m,Ω_Λ=Ω_Λ,Ω_r=Ω_r,H0=H0),k=k,burst=burst,q=q)
end

const z_rec = 1089
const a_rec = 1/(1+z_rec)
const t_rec = a_to_t(a_rec)
const t_0 = a_to_t(1.0)

const cosmo_x = 10 .^ LinRange(-20,20,10000)
const cosmo_y = a_to_t.(cosmo_x)
a_to_t_interp = linear_interpolation(cosmo_x,cosmo_y,extrapolation_bc=Line())
t_to_a_interp = linear_interpolation(cosmo_y,cosmo_x,extrapolation_bc=Line())

function a_to_t_fast(a)
    if a > cosmo_x[1] && a < cosmo_x[end]
        i = searchsortedfirst(cosmo_x,a)
        return cosmo_y[i-1]+(a-cosmo_x[i-1])/(cosmo_x[i]-cosmo_x[i-1])*(cosmo_y[i]-cosmo_y[i-1])
    elseif a <= cosmo_x[1]
        return cosmo_y[1]*(a/cosmo_x[1])^2
    else
        return log(a/cosmo_x[end])/H0_fid/sqrt(Ω_Λ_fid) + cosmo_y[end]
    end
end

function t_to_a_fast(t)
    if t > cosmo_y[1] && t < cosmo_y[end]
        i = searchsortedfirst(cosmo_y,t)
        return cosmo_x[i-1]+(t-cosmo_y[i-1])/(cosmo_y[i]-cosmo_y[i-1])*(cosmo_x[i]-cosmo_x[i-1])
    elseif t <= cosmo_y[1]
        return cosmo_x[1]*(t/cosmo_y[1])^(1/2)
    else
        return exp.((t-cosmo_y[end])*H0_fid*sqrt(Ω_Λ_fid))*cosmo_x[end]
    end
end
    
function z_to_t_fast(z)
    return a_to_t_fast.(1 ./ (1 .+ z))
end

function t_to_z_fast(t)
    return 1 ./t_to_a_fast.(t) .- 1
end

function g_e_integral(u,z,sign)
    return u.^2 .* sqrt.(u.^2 .- z.^2) ./ (exp.(u) .+ sign)
end


function g_e_fermion(x)
    integrand(u,p) = g_e_integral(u,1/x,1)
    prob = IntegralProblem(integrand,1/x,Inf)
    sol = solve(prob,QuadGKJL())
    return 15/pi^4*sol.u
end

function g_e_boson(x)
    integrand(u,p) = g_e_integral(u,1/x,-1)
    prob = IntegralProblem(integrand,1/x,Inf)
    sol = solve(prob,QuadGKJL())
    return 15/pi^4*sol.u
end

function g_p_integral(u,z,sign)
    return (u.^2 .- z.^2).^(3/2) ./ (exp.(u) .+ sign)
end


function g_p_fermion(x)
    integrand(u,p) = g_p_integral(u,1/x,1)
    prob = IntegralProblem(integrand,1/x,Inf)
    sol = solve(prob,QuadGKJL())
    return 15/pi^4*sol.u
end

function g_p_boson(x)
    integrand(u,p) = g_p_integral(u,1/x,-1)
    prob = IntegralProblem(integrand,1/x,Inf)
    sol = solve(prob,QuadGKJL())
    return 15/pi^4*sol.u
end

function g_s_fermion(x)
    return (3*g_e_fermion(x) + g_p_fermion(x))/4
end

function g_s_boson(x)
    return (3*g_e_boson(x) + g_p_boson(x))/4
end

