#Interpolation routine for rad_correct(a)

function rad_correct(a)
    if a > rad_correct_x[1] && a < rad_correct_x[end]
        i = searchsortedfirst(rad_correct_x,a)
        return rad_correct_y[i-1]+(a-rad_correct_x[i-1])/(rad_correct_x[i]-rad_correct_x[i-1])*(rad_correct_y[i]-rad_correct_y[i-1])
    elseif a <= rad_correct_x[1]
        return rad_correct_y[1]
    else
        return rad_correct_y[end]
    end
end


#Interpolation routines for f(M,a) and g(M,a)

"""
    fM_int(m)

Returns the interpolated value of f(m,a=0)

Input:

    m : Black hole mass in g
"""
function fM_int(m)
    if m > m_tab[1] && m < m_tab[end]
        i = searchsortedfirst(m_tab,m)
        return fM_tab[i-1,1]+(m-m_tab[i-1])/(m_tab[i]-m_tab[i-1])*(fM_tab[i,1]-fM_tab[i-1,1])
    elseif m <= m_tab[1]
        return fM_tab[1,1]
    else
        return fM_tab[end,1]
    end
end

"""
    fM_int(m,a)

Returns the interpolated value of f(m,a)

Input:

    m : Black hole mass in g
    a : Dimensionless black hole spin parameter
"""
function fM_int(m,a)
    
    if a == 0
        return fM_int(m)
    end
    
    if m > m_tab[1] && m < m_tab[end]
        i_m = searchsortedfirst(m_tab,m)
    elseif m <= m_tab[1]
        i_m = 2
        m = m_tab[1]
    else
        i_m = length(m_tab)
        m = m_tab[end]
    end
    
    if a > a_tab[1] && a < a_tab[end]
        i_a = searchsortedfirst(a_tab,a)
    elseif a <= a_tab[1]
        i_a = 2
        a = a_tab[1]
    else
        i_a = length(a_tab)
        a = a_tab[end]
    end

    norm = 1/(a_tab[i_a]-a_tab[i_a-1])/(m_tab[i_m]-m_tab[i_m-1])
    f1 = fM_tab[i_m-1,i_a-1]*(m_tab[i_m]-m)*(a_tab[i_a]-a)
    f2 = fM_tab[i_m,i_a-1]*(m-m_tab[i_m-1])*(a_tab[i_a]-a)
    f3 = fM_tab[i_m-1,i_a]*(m_tab[i_m]-m)*(a-a_tab[i_a-1])
    f4 = fM_tab[i_m,i_a]*(m-m_tab[i_m-1])*(a-a_tab[i_a-1])
    return norm*(f1+f2+f3+f4)
end

"""
    gM_int(m)

Returns the interpolated value of g(m,a=0)

Input:

    m : Black hole mass in g
"""
function gM_int(m)
    if m > m_tab[1] && m < m_tab[end]
        i = searchsortedfirst(m_tab,m)
        return gM_tab[i-1,1]+(m-m_tab[i-1])/(m_tab[i]-m_tab[i-1])*(gM_tab[i,1]-gM_tab[i-1,1])
    elseif m <= m_tab[1]
        return gM_tab[1,1]
    else
        return gM_tab[end,1]
    end
end

"""
    gM_int(m,a)

Returns the interpolated value of g(m,a)

Input:

    m : Black hole mass in g
    a : Dimensionless black hole spin parameter
"""
function gM_int(m,a)
    
    if a == 0
        return fM_int(m)
    end
    
    if m > m_tab[1] && m < m_tab[end]
        i_m = searchsortedfirst(m_tab,m)
    elseif m <= m_tab[1]
        i_m = 2
        m = m_tab[1]
    else
        i_m = length(m_tab)
        m = m_tab[end]
    end
    
    if a > a_tab[1] && a < a_tab[end]
        i_a = searchsortedfirst(a_tab,a)
    elseif a <= a_tab[1]
        i_a = 2
        a = a_tab[1]
    else
        i_a = length(a_tab)
        a = a_tab[end]
    end

    norm = 1/(a_tab[i_a]-a_tab[i_a-1])/(m_tab[i_m]-m_tab[i_m-1])
    f1 = gM_tab[i_m-1,i_a-1]*(m_tab[i_m]-m)*(a_tab[i_a]-a)
    f2 = gM_tab[i_m,i_a-1]*(m-m_tab[i_m-1])*(a_tab[i_a]-a)
    f3 = gM_tab[i_m-1,i_a]*(m_tab[i_m]-m)*(a-a_tab[i_a-1])
    f4 = gM_tab[i_m,i_a]*(m-m_tab[i_m-1])*(a-a_tab[i_a-1])
    return norm*(f1+f2+f3+f4)
end


#Interpolation routines for the greybody factors γ(x,s,a)

"""
    greybody_int(x)

Returns the interpolated value of Γ/(exp(x)-(-1)^s)

Input:

    x : e/(4*pi*k_B*T) = 2*G*m*e/(hbar*c^3)
    s : spin
"""
function greybody_int(x,s)
    if s == 0
        γ = γ_0_tab
        γ_fit = γ_0_fit
    elseif s == 0.5
        γ = γ_05_tab
        γ_fit = γ_05_fit
    elseif s == 1
        γ = γ_1_tab
        γ_fit = γ_1_fit
    elseif s == 1.5
        γ = γ_15_tab
        γ_fit = γ_15_fit
    elseif s == 2
        γ = γ_2_tab
        γ_fit = γ_2_fit
    end
    
    if x > x_tab[1] && x < x_tab[end]
        i_x = searchsortedfirst(x_tab,x)
        return γ[1,i_x-1]+(x-x_tab[i_x-1])/(x_tab[i_x]-x_tab[i_x-1])*(γ[1,i_x]-γ[1,i_x-1])
    elseif x <= x_tab[1]
        return 10.0 ^ (γ_fit[1,1]*log10(x) + γ_fit[1,2]) 
    else
        return 10.0 ^ (x*γ_fit[1,3] + γ_fit[1,4] + γ_fit[1,5]*cos(x*γ_fit[1,7]) + γ_fit[1,6]*sin(x*γ_fit[1,7]))
    end
end

"""
    greybody_int(x)

Returns the interpolated value of Γ/(exp(x)-(-1)^s)

Input:

    x : e/(4*pi*k_B*T) = 2*G*m*e/(hbar*c^3)
    s : spin
    a : Dimensionless black hole spin parameter
"""
function greybody_int(x,s,a)
    
    if a == 0
        return greybody_int(x,s)
    end
    
    if s == 0
        γ = γ_0_tab
        γ_fit = γ_0_fit
    elseif s == 0.5
        γ = γ_05_tab
        γ_fit = γ_05_fit
    elseif s == 1
        γ = γ_1_tab
        γ_fit = γ_1_fit
    elseif s == 1.5
        γ = γ_15_tab
        γ_fit = γ_15_fit
    elseif s == 2
        γ = γ_2_tab
        γ_fit = γ_2_fit
    end
    
    if a < a_t[end]
        i_a = searchsortedfirst(a_t,a)
    else
        i_a = length(a_t)
    end
    
    
    if x > x_tab[1] && x < x_tab[end]
        i_x = searchsortedfirst(x_tab,x)
        if a > a_t[1] && a < a_t[end]
            i_a = searchsortedfirst(a_t,a)
        elseif a <= a_t[1]
            i_a = 2
            a = a_t[1]
        else
            i_a = length(a_t)
            a = a_t[end]
        end
    
    
        norm = 1/(x_tab[i_x]-x_tab[i_x-1])/(a_t[i_a]-a_t[i_a-1])
        f1 = γ[i_a-1,i_x-1]*(a_t[i_a]-a)*(x_tab[i_x]-x)
        f2 = γ[i_a,i_x-1]*(a-a_t[i_a-1])*(x_tab[i_x]-x)
        f3 = γ[i_a-1,i_x]*(a_t[i_a]-a)*(x-x_tab[i_x-1])
        f4 = γ[i_a,i_x]*(a-a_t[i_a-1])*(x-x_tab[i_x-1])
        return norm*(f1+f2+f3+f4)
        
    elseif x <= x_tab[1]
        fit = γ_fit[i_a-1,1:2] .+ (a-a_t[i_a-1])/(a_t[i_a]-a_t[i_a-1])*(γ_fit[i_a,1:2] .- γ_fit[i_a-1,1:2])
        return 10.0 ^ (fit[1]*log10(x) + fit[2]) 
    else
        fit = γ_fit[i_a-1,:] .+ (a-a_t[i_a-1])/(a_t[i_a]-a_t[i_a-1])*(γ_fit[i_a,:] .- γ_fit[i_a-1,:])
        return 10.0 ^ (x*fit[3] + fit[4] + fit[5]*cos(x*fit[7]) + fit[6]*sin(x*fit[7]))
    end
end