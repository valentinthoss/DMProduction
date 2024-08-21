function loginf(x)
    if x>0
        return log(x)
    else
        return -Inf
    end
end

function relu(x)
    return ifelse(x>0,x,0)
end

"""
    solid_angle(l_min,l_max,b_min,b_max)

Returns the solid angle (in radians) covered by galactic coordinates (l,b).
"""
function solid_angle(l_min,l_max,b_min,b_max)
    integrand(angles,p) = sin(angles[2])
    low = [Float64(l_min),Float64(b_min)+pi/2]
    up = [Float64(l_max),Float64(b_max)+pi/2]
    prob = IntegralProblem(integrand,low,up)
    sol = solve(prob, HCubatureJL())
    return sol.u
end