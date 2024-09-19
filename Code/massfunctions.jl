function confinement_get_λΛ3(f,M)
    return 3*f/(32*pi)*(M/meq)^(1/2)
end

function confinement_get_f(λΛ3,M)
    return 32*pi/3*λΛ3*(M/meq)^(-1/2)
end

function confinement_get_spin(M,Λ)
    M = M/mp
    return 1 / (M * Λ) * log(M/Λ)^(1/2)
end

function confinement_get_dndm(M,M0)
    if M > M0
        return 3/2*M0^(3/2) / M^(5/2)
    else
        return 0.0
    end
end

struct lognormal
    σ
end

lognormal() = lognormal(1.0)
mean(d::lognormal,x0) = x0*exp(d.σ^2/2)
pdf(d::lognormal,x,x0) = 1 ./ (sqrt(2*pi)*d.σ*x) .* exp.(-(log.(x) .- log(x0)).^2/(2*d.σ^2))

struct skew_lognormal
    σ
    α
end

skew_lognormal() = skew_lognormal(1.0,0.0)
pdf(d::skew_lognormal,x,x0) = 1 ./ (sqrt(2*pi)*d.σ*x) .* exp.(-(log.(x) .- log(x0)).^2/(2*d.σ^2)) .* (1 .+ erf.(d.α/sqrt(2)/d.σ*log.(x/x0)))

struct critical_collapse
    n
end

critical_collapse() = critical_collapse(1/0.36)
mean(d::critical_collapse,x0) = x0*gamma(d.n+1)/gamma(d.n)
pdf(d::critical_collapse,x,x0) = d.n ./ (x0*gamma(1 + 1/d.n)) .* (x/x0).^(d.n) .* exp.(-(x/x0).^(d.n))

struct critical_collapse_n
    n
end

critical_collapse_n() = critical_collapse_n(1/0.36)
mean(d::critical_collapse_n,x0) = x0*gamma(1 + 1/d.n)
pdf(d::critical_collapse_n,x,x0) = d.n ./ x0 .* (x/x0).^(d.n-1) .* exp.(-(x/x0).^(d.n))

struct generalized_critical_collapse
    α
    β
end

generalized_critical_collapse() = generalized_critical_collapse(1/0.36,1/0.36)
function pdf(d::generalized_critical_collapse,x,x0)
    if (d.α+1)/d.β < 10
        return d.β/x0/gamma((d.α+1)/d.β)*(d.α/d.β)^((d.α+1)/d.β)*(x/x0).^(d.α) .* exp.(-d.α/d.β*(x/x0).^d.β)
    else
        return 1/x0*(d.β*(d.α+1-d.β)/(2*pi))^(1/2)*exp((d.α+1)/d.β - 1)*(d.α/(d.α+1-d.β))^((d.α+1)/d.β)*(x/x0).^(d.α) .* exp.(-d.α/d.β*(x/x0).^d.β)
    end
end

struct generalized_critical_collapse_n
    α
    β
end

generalized_critical_collapse_n() = generalized_critical_collapse_n(1/0.36,1/0.36)
mean(d::generalized_critical_collapse_n,x0) = x0*gamma((1+d.α)/d.β)/gamma(d.α/d.β)*(d.β/d.α)^(1/d.β)
pdf(d::generalized_critical_collapse_n,x,x0) = d.α/x0/gamma(d.α/d.β)*(d.α/d.β)^(d.α/d.β-1)*(x/x0).^(d.α-1) .* exp.(-d.α/d.β*(x/x0).^d.β)

struct powerlaw
    n
    xmax
end

powerlaw(n) = powerlaw(n,Inf)

function pdf(d::powerlaw,x,x0)
    if d.n > -1.0
        return ifelse.(x .< x0,(d.n+1.0)/float(x0)^(d.n+1.0) * x.^float(d.n),0.0)
    elseif d.n < -1.0
        return ifelse.(d.xmax .> x .> x0,-(d.n+1.0)/float(x0)^(d.n+1.0) * x.^float(d.n),0.0)
    elseif d.n == -1.0
        if d.xmax != Inf 
            return ifelse.(d.xmax .> x .> x0,1/log(d.xmax/x0) * x.^float(d.n),0.0)
        else
            throw(DomainError(d.n, "Slope of -1 chosen but no xmax defined. Distribution cannot be normalized."))
        end
    end
end
    
function mean(d::powerlaw,x0)
    if d.n > -1.0
        return x0 * (d.n+1)/(d.n+2)
    elseif d.n !== -2.0 && d.n !== -1.0
        if d.xmax != Inf
            return x0 * (d.n+1)/(d.n+2) - float(d.xmax)^(d.n+2.0)/float(x0)^(d.n+1.0) * (d.n+1)/(d.n+2)
        elseif d.n < -2.0
            return x0 * (d.n+1)/(d.n+2)
        else
            throw(DomainError(d.n, "Slope between -2 and -1 chosen but no xmax defined. Mean is not finite."))
        end
    elseif d.n == -2.0
        if d.xmax != Inf
            return -(d.n+1)/float(x0)^(d.n+1.0) * log(d.xmax/x0)
        else
            throw(DomainError(d.n, "Slope of -2 chosen but no xmax defined. Mean is not finite."))
        end
    elseif d.n == -1.0
        if d.xmax !== Inf
            return 1/log(d.xmax/x0) * (d.xmax - x0)
        else
            throw(DomainError(d.n, "Slope of -1 chosen but no xmax defined. Mean is not finite."))
        end
    end
end

function mass_function_evolution(m,t,mass_function,m_c;k=0,burst=true,q=0.5,fargs...)
    m0 = get_m_init.(m,t,k=k,burst=burst,q=q)
    if k == 0.0
        return mass_function.(m0,m_c;fargs...) .* mass_loss_rate.(m0) ./ mass_loss_rate.(m)
    else
        return mass_function.(m0,m_c,fargs...) .* grad_m_m0(m,m0,q,k,burst)
    end
end

function convolve_constraints(m_mc,f_mc,m_spectrum,mass_func)
    lf = log.(1 ./ f_mc)
    f_int = linear_interpolation(log.(m_mc),lf,extrapolation_bc=-Inf)
    integrand = mass_func.* exp.(f_int(log.(m_spectrum)))
    if argmax(integrand) == 1 || argmax(integrand) == length(integrand)
        #printlognormal("The mass function could not be properly integrated. powerlawease extend the integration range!")
    end
    mask = .!isnan.(integrand)
    return 1 / trapz(m_spectrum[mask],integrand[mask])
end

function convolve_constraints3(m_mc,f_mc,mass_func,i_t,m,k)
    lf = log.(mass_func ./ f_mc)
    if argmax(lf) == 1 || argmax(lf) == length(lf) || lf[1] > maximum(lf)*0.01
        #printlognormal("The mass function could not be properly integrated for m=$(m) and k=$(k). powerlawease extend the integration range!")
    end
    if isnothing(i_t)             
        integrand = mass_func ./ f_mc                  
        mask = .!isnan.(integrand)
        return 1 / trapz(m_mc[mask],integrand[mask])
    else
        integrand1 = mass_func[1:i_t-1] ./ f_mc[1:i_t-1]              
        mask = .!isnan.(integrand1)    
        sol1 = trapz(m_mc[1:i_t-1][mask],integrand1[mask])
        integrand2 = mass_func[i_t:end] ./ f_mc[i_t:end]              
        mask = .!isnan.(integrand2)    
        sol2 = trapz(m_mc[i_t:end][mask],integrand2[mask])
        return 1 / (sol1 + sol2)
    end
end

function convolve_constraints_integrand(x,interp)
    y = exp(interp(x)) * exp(x)
    if !isnan(y)
        return y
    else
        return 0.0
    end
end

function convolve_constraints2(m_mc,f_mc,mass_func,i_t,m,k)
    lf = log.(mass_func ./ f_mc)
    lf[isnan.(lf)] .= -Inf
    if argmax(lf) == 1 || argmax(lf) == length(lf) || lf[1] > maximum(lf)*0.01
        printlognormal("The mass function could not be properly integrated for m=$(m) and k=$(k). powerlawease extend the integration range!")
    end
    interp = linear_interpolation(log.(m_mc),lf,extrapolation_bc=-Inf)
    integrand(x,p) = convolve_constraints_integrand(x,interp)
    if isnothing(i_t)                                   
        prob = IntegralProblem(integrand,log(m_mc[1]),log(m_mc[end]))
        sol = solve(prob,QuadGKJL())
        return 1 / sol.u
    else
        prob1 = IntegralProblem(integrand,log(m_mc[i_t]),log(m_mc[end]))
        sol1 = solve(prob1,QuadGKJL())
        prob2 = IntegralProblem(integrand,log(m_mc[1]),log(m_mc[i_t-1]))
        sol2 = solve(prob2,QuadGKJL())
        return 1 / (sol1.u + sol2.u)
    end
end

function grad_m_m0(m,m0,q,k,burst=true)
    if k == 0
        return mass_loss_rate.(m0) ./ mass_loss_rate.(m)
    elseif burst
        return ifelse.(m .<= q*m0,1/q*fM_int.(q*m0)./fM_int.(m).*(m./(q*m0)).^(2 + 2*k),fM_int.(m0)./fM_int.(m).*(m./(m0)).^2)
    else
        return grad(m0,m)
    end
end

function grad(y,x)
    g = (y[3:end] .- y[1:end-2]) ./ (x[3:end] .- x[1:end-2])
    g_start = (y[2] - y[1])/(x[2]-x[1])
    g_end = (y[end] - y[end-1])/(x[end]-x[end-1])
    return vcat(g_start,g,g_end)
end