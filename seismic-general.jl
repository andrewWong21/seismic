function c(r)
    # c = r^n, where n > 0:
    # violates Herglotz conditions (r/c is always increasing)
    # r/c = r/r^n = 1/r^(n-1), always decreasing

    # c = 1:
    # r/c = r, always increasing

    return r^-1
end

function f(p, r)
    valC = c(r)
    return p / (r * sqrt( (r/valC)^2 - p^2 ) )
end

# sanity check
# original equation:
# p / ( r * sqrt((r/c)^2 - p^2) )
# c = 1/r
# then r/c = r/(1/r) = r^2
# then (r/c)^2 = (r^2)^2 = r^4
# then (r/c)^2 - p^2 = r^4 - p^2
# r^4 - p^2 = 0, r^4 = p^2
# (r^2)^2 = p^2, r^2 = p
# r_s = sqrt(p)


# general case
# (r/c)^2 - p^2 = 0
# r^2/c^2 = p^2
# r/c = p, solve for r
# r = pc
# find zeroes for r/c - p = 0 to get turning radius
# p is given, c=c(r), r=?

# da/dp measures curvature of wavefront
# alpha in degrees, p is angular momentum

function getRadius(p)
    func(r) = r/(c(r)) - p
    print("turning radius: ")
    println(find_zero(func, 0.5))
    return find_zero(func, 0.5)
end

# Herglotz condition: r/c is always increasing
# max value 1/c, c evaluated at 1

function checkBoundary(p, R)
    # determines if path is turning or reflecting off inner boundary
    # if p > R/c(R), path is turning smoothly
    # if p < R/c(R), path is reflecting off inner boundary

    # general case
    # grazingP = R / c(R)

    # no inner boundary
    if (R == 0)
        return p
    end

    # threshold/grazing p-value
    grazingP = R / c(R)

    print("p: ")
    println(p)
    print("grazingP: ")
    println(grazingP)

    if (p < grazingP)
        println("p < grazingP, reflect\n")
        return grazingP
    elseif (p > grazingP)
        println("p > grazingP, smooth turn\n")
    else
        println("p = grazingP\n")
    end
    return p
end

function getAlpha(p)
    minP = checkBoundary(p, R)
    a = getRadius(minP) # lower bound of integral may change
    ans = midpoint(a, 1, p)
    return ans
end

function midpoint(a, b, p)
    #print("integral goes from ")
    #print(a)
    #print(" to ")
    #println(b)
    sum = 0
    h = (b - a)/30 # step size
    for num in 0:29
        #println("Step: " * string(num))
        sum += f(p, a + num * h + h/2)
    end
    return sum * h
end

function fillRadiusValues(p, minRadius)
    h = (1 - minRadius)/30
    radius_values = collect(minRadius:h:1)
    radius_path = vcat(reverse(radius_values), radius_values[2:31])
    return radius_path
end

function fillThetaValues(p, alpha, radius_values, minRadius)
    turningRadius = getRadius(p)
    theta_values = zeros(Float64, 31)
    for num in 1:30
        setindex!(theta_values, midpoint(minRadius, getindex(radius_values, num), p) + (pi/2 - alpha), num)
    end
    setindex!(theta_values, (pi/2 - alpha), 31) # pi/2
    
    theta_path = vcat((theta_values), zeros(Float64, 30))
    for num in 2:31
        setindex!(theta_path, -1 * getindex(reverse(theta_values), num) + 2 * (pi/2 - alpha), num + 30)
    end
    return theta_path
end

function fillXValues(radius_path, theta_path)
    x_values = zeros(Float64, 61)
    for num in 1:61
        setindex!(x_values, getindex(radius_path, num) * cos(getindex(theta_path, num)), num)
    end
    return x_values
end

function fillYValues(radius_path, theta_path)
    y_values = zeros(Float64, 61)
    for num in 1:61
        setindex!(y_values, getindex(radius_path, num) * sin(getindex(theta_path, num)), num)
    end
    return y_values
end

function reflect(radius_path, theta_path, alpha, r)
    if (r > 0)
        xshifted = zeros(Float64, length(radius_path))
        yshifted = zeros(Float64, length(radius_path))
        for num in 1:r
            for ind in 1:length(radius_path)
                setindex!(xshifted, getindex(radius_path, ind) * cos(getindex(theta_path, ind) - (2 * num * alpha)), ind)
                setindex!(yshifted, getindex(radius_path, ind) * sin(getindex(theta_path, ind) - (2 * num * alpha)), ind)
            end
            plot!(xshifted, yshifted, color=:blue, label="Reflection")
        end
    end
end

using Roots
R = 0.5 # inner boundary radius

findGrazing = false
compute = true
showGraph = true

if (findGrazing)
    windingNum = 2   # n
    numReflects = 35 # m

    # a(p) <= a(grazingP)
    # m * a(p) = pi * n
    # m = (pi * n) / (a(p)) >= (pi * n) / (a(grazingP))
    # (pi * n) / (a(grazingP)) is minimum value for m (# of reflections)

    # numReflects and windingNum should be relatively prime

    closing(p) = numReflects/pi * getAlpha(p, 1) - windingNum

    println("Enter initial guess: ")
    initialGuess = readline()
    initialGuess = parse(Float64, initialGuess)
    println(find_zero(closing, initialGuess))
end

if (compute)
    println("Enter number of p values to plot (n >= 1): ")
    numvalues = readline()
    numvalues = parse(Int64, numvalues)

    p_values = zeros(Float64, numvalues)
    for num in 1:length(p_values)
        println("Enter p (0 < p < 1): ")
        p = readline()
        p = parse(Float64, p)
        setindex!(p_values, p, num)
    end

    alpha_values = Array{Float64}(undef, numvalues)
    radius_values_array = Array{Array{Float64}}(undef, numvalues)
    theta_values_array = Array{Array{Float64}}(undef, numvalues)
    x_values_array = Array{Array{Float64}}(undef, numvalues)
    y_values_array = Array{Array{Float64}}(undef, numvalues)
    
    for num in 1:length(radius_values_array)
        p_value = getindex(p_values, num)
        minRadius = getRadius(checkBoundary(p_value, R))

        setindex!(alpha_values, getAlpha(p_value), num)
        #println("Alpha OK (changed)")
        setindex!(radius_values_array, fillRadiusValues(p_value, minRadius), num)
        #println("Radii OK (changed)")
        setindex!(theta_values_array, fillThetaValues(p_value, getindex(alpha_values, num), getindex(radius_values_array, num), minRadius), num)
        #println("Theta OK (changed)")
        setindex!(x_values_array, fillXValues(getindex(radius_values_array, num), getindex(theta_values_array, num)), num)
        setindex!(y_values_array, fillYValues(getindex(radius_values_array, num), getindex(theta_values_array, num)), num)
    end

    for num in 1:length(p_values)
        print("p value: ")
        println(getindex(p_values, num))
        println(getindex(alpha_values, num))
        println()
    end

end

if (showGraph)
    println("Enter number of additional reflections for first path (r >= 0): ")
    r = readline()
    r = parse(Int64, r)
    #r = 

    using Plots
    plotly()
    plot(xlims=(-1, 1), ylims=(-1, 1), aspect_ratio=:equal, size=(1600, 900))
    # plot unit circle
    m(t) = cos(t)
    n(t) = sin(t)
    plot!(m, n, 0, 2pi, label="Surface", color=:black)
    if (R != 0)
         mini_m(t) = R * sin(t)
         mini_n(t) = R * cos(t)
        plot!(mini_m, mini_n, 0, 2pi, label="Inner Boundary", color=:red)
    end

    for num in 1:length(radius_values_array)
        plot!(getindex(x_values_array, num), getindex(y_values_array, num), label="p = "*string(getindex(p_values, num)))
    end

    if (r > 0)
        reflect(getindex(radius_values_array, 1), getindex(theta_values_array, 1), getindex(alpha_values, 1), r)
    end

    gui()
end