function cOuter(r)
    return 0.5 * r^-2
end

function cInner(r)
    return 0.5 * r^-1
end
#c1 and c2 functions depending on inside/outside boundary

function f(p, r, isOuterLayer)  # add layer input to determine which c-function to use
    valC = 0;
    if (isOuterLayer)
        valC = cOuter(r)
    else
        valC = cInner(r)
    end
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

function getRadius(p, isOuterLayer)
    # minimum turning radius is given by the grazing p value
    if (isOuterLayer)
        graze1(r) = r/(cOuter(r)) - p
        turningRadius = find_zero(graze1, 0.5)
    else
        graze2(r) = r/(cInner(r)) - p
        turningRadius = find_zero(graze2, 0.5)
    end
    return turningRadius
end

# Herglotz condition: r/c is always increasing
# max value 1/c, c evaluated at 1

function checkBoundary(p, R, isOuterLayer)
    # determines if path is turning or reflecting off inner boundary
    # grazingP = R / c(R)
    # if p > grazingP, path is turning smoothly
    # if p < grazingP, path is reflecting off inner boundary

    # two layers, two grazing P values outerGrazingP and innerGrazingP, innerGrazingP < outerGrazingP
    # if p >= outerGrazingP, turning in outer layer (assuming it starts at outer layer)
    # if innerGrazingP < p < outerGrazingP, total internal reflection (does not enter inner layer)
    # if p < innerGrazingP, transmission from outer layer to inner layer
    # if p = innerGrazingP, headwave: critical angle, glides along interface

    # no inner boundary, indicates inner layer
    if (R == 0)
        println("inner boundary is 0 (inner layer)\n")
        return p
    end
    # threshold/grazing p-value
    if (isOuterLayer)
        println("using outer grazingP")
        grazingP = R / cOuter(R)
    #else
        #println("using inner grazingP")
        #grazingP = R / cInner(R)
    end
    
    print("p: ")
    println(p)
    print("grazingP: ")
    println(grazingP)
    
    # for reflections
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

function getAlpha(p, R, outerRadius, isOuterLayer)
    println("in getAlpha")
    minP = checkBoundary(p, R, isOuterLayer) # 
    innerRadius = getRadius(minP, isOuterLayer) # lower bound of integral
    print("innerRadius: ")
    println(innerRadius)
    println()
    return midpoint(innerRadius, outerRadius, p, isOuterLayer)
end

function midpoint(a, b, p, isOuterLayer) # add layer input to determine which c-function to use
    sum = 0
    h = (b - a)/30 # step size
    for num in 0:29
        sum += f(p, a + num * h + h/2, isOuterLayer)
    end
    return sum * h
end

function fillRadiusValues(p, minRadius, maxRadius)
    h = (maxRadius - minRadius)/30
    radius_values = collect(minRadius:h:maxRadius)
    radius_path = vcat(reverse(radius_values), radius_values[2:31])
    return radius_path
end

function fillThetaValues(p, alpha, radius_values, minRadius, isOuterLayer, offset)
    #turningRadius = getRadius(p, true)
    theta_values = zeros(Float64, 31)
    for num in 1:30
        setindex!(theta_values, midpoint(minRadius, getindex(radius_values, num), p, isOuterLayer) + (pi/2 - alpha), num)
    end
    setindex!(theta_values, (pi/2 - alpha), 31) # pi/2
    
    theta_path = vcat((theta_values), zeros(Float64, 30))
    for num in 2:31
        setindex!(theta_path, -1 * getindex(reverse(theta_values), num) + 2 * (pi/2 - alpha), num + 30)
    end

    for num in 1:61
        setindex!(theta_path, getindex(theta_path, num) - offset, num)
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

findGrazing = false
compute = false
showGraph = false

if (findGrazing)
    #windingNum = 2   # n
    #numReflects = 35 # m

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
        setindex!(radius_values_array, fillRadiusValues(p_value, minRadius), num)
        setindex!(theta_values_array, fillThetaValues(p_value, getindex(alpha_values, num), getindex(radius_values_array, num), minRadius), num)
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

# -------------------------------------------------------------------------------------------------

function seismicPath(p, drawLayers, reflectionNumber) # 
    
    # innerGrazingP = 0.25, outerGrazingP = 0.5
    # if p >= outerGrazingP, turning in outer layer (assuming it starts at outer layer)            0.75  ooooo
    # if innerGrazingP < p < outerGrazingP, total internal reflection (does not enter inner layer) 0.3   ooooo
    # if p = innerGrazingP, headwave: critical angle, glides along interface                       0.25  
    # if p < innerGrazingP, transmission from outer layer to inner layer                           0.2   oiiio

    # if path reflects in outer layer, 
    # iterate through drawLayers string
    # update offset value using either 1 alpha (half path) or 2 alpha (full path)

    # use parity of O's - count occurrences
    # draw first half if odd, second half if even (assuming path is reflecting)

    # recursive option
    # seismicPath returns an offset value
    # draw first half of reflecting path for outer layer
    # calls itself for the inner layer, uses offset value

    interface = R
    surface = 1
    fullPathRadiusVals = Float64[]
    fullPathThetaVals = Float64[]
    innerOffset = 0
    if (cmp(drawLayers[1], 'o') == 0) # first path in outer layer
        outerOffset = 0
        outerParity = 0
    else # first path in inner layer, measures for fixing exiting half of outer layer path
        outerOffset = -1 * getAlpha(p, interface, surface, true)
        outerParity = 1
    end

    for layer in drawLayers
        alpha = 0
        radius_values = Array{Float64}(undef, 61)
        theta_values = Array{Float64}(undef, 61)
        x_values = Array{Float64}(undef, 61)
        y_values = Array{Float64}(undef, 61)

        isOuterLayer = (layer == 'o')
        if (isOuterLayer)
            outerParity += 1
            #println("drawing in outer layer")
            
            #println("getting minRadius")
            pValue = checkBoundary(p, interface, isOuterLayer)
            minRadius = getRadius(pValue, isOuterLayer)
            maxRadius = surface

            #println("getting alpha")
            alpha = getAlpha(p, interface, maxRadius, isOuterLayer)
            radius_values = fillRadiusValues(p, minRadius, maxRadius)
            theta_values = fillThetaValues(p, alpha, radius_values, minRadius, isOuterLayer, outerOffset)
        else # inner layer
            #println("drawing in inner layer")

            #println("getting minRadius")
            minRadius = getRadius(p, isOuterLayer)
            maxRadius = interface

            #println("getting alpha")
            alpha = getAlpha(p, 0, maxRadius, isOuterLayer)
            radius_values = fillRadiusValues(p, minRadius, maxRadius)
            theta_values = fillThetaValues(p, alpha, radius_values, minRadius, isOuterLayer, innerOffset)
        end

        x_values = fillXValues(radius_values, theta_values)
        y_values = fillYValues(radius_values, theta_values)

        if (!isOuterLayer) # drawing in inner layer
            println("drawing inner")
            fullPathRadiusVals = [fullPathRadiusVals; radius_values]
            fullPathThetaVals = [fullPathThetaVals; theta_values]
            plot!(x_values, y_values, label="inner", color=:darkorange)
            innerOffset += 2 * alpha
            outerOffset += 2 * alpha

        elseif (p < outerGrazingP) #drawing in outer layer, reflecting
            println("outer reflecting")
            if (outerParity % 2 == 0) # only drawing second half of path
                println("")
                fullPathRadiusVals = [fullPathRadiusVals; radius_values[31:61]]
                fullPathThetaVals = [fullPathThetaVals; theta_values[31:61]]
                plot!(x_values[31:61], y_values[31:61], label="outer, second half", color=:red)
                innerOffset += alpha
                outerOffset += 2 * alpha

            else # only drawing first half of path
                println("outer first half")
                fullPathRadiusVals = [fullPathRadiusVals; radius_values[1:31]]
                fullPathThetaVals = [fullPathThetaVals; theta_values[1:31]]
                plot!(x_values[1:31], y_values[1:31], label="outer, first half", color=:green)
                innerOffset += alpha
            end
        else # drawing in outer layer, turning
            println("outer turning")
            fullPathRadiusVals = [fullPathRadiusVals; radius_values]
            fullPathThetaVals = [fullPathThetaVals; theta_values]
            pathLabel = "p = " * string(p)
            plot!(x_values, y_values, label=pathLabel, color=:gray)
            outerOffset += 2 * alpha
        end # path plotting
    end # for loop

    reflectTotalPath(fullPathRadiusVals, fullPathThetaVals, getEpicentralDistance(p, drawLayers), reflectionNumber)
end

function reflectTotalPath(radius_path, theta_path, alpha, r)
    if (r > 0)
        xshifted = zeros(Float64, length(radius_path))
        yshifted = zeros(Float64, length(radius_path))
        for num in 1:r
            for ind in 1:length(radius_path)
                setindex!(xshifted, getindex(radius_path, ind) * cos(getindex(theta_path, ind) - (num * alpha)), ind)
                setindex!(yshifted, getindex(radius_path, ind) * sin(getindex(theta_path, ind) - (num * alpha)), ind)
            end
            plot!(xshifted, yshifted, color=:blue, label="Reflection")
        end
    end
end

    # determine total offset for one full path (epicentral distance)
    # m * d_p1 = 2pi * n
    # d_p1 = 2*a_o + 2*a_i*r
        # r: number of inner paths
    # m : number of full paths
    # n : winding number
    # use with findZeroes to get p value that results in a closing path

    function getEpicentralDistance(p, drawLayers)
        println("in getEpicentralDistance")
        print("current p-value: ")
        println(p)
        numInnerPaths = 0
        numOuterPaths = 0
        for layer in drawLayers
            if (layer == 'o')
                numOuterPaths += 1
            elseif (layer == 'i')
                numInnerPaths += 1
            end
        end

        innerAlpha = 0
        outerAlpha = 0
        interface = R
        surface = 1
    
        #= possibilities
        outerGrazingP = 0.5
        p > outerGrazingP (p = 0.6)
            o
        p < outerGrazingP
            if innerGrazingP < outerGrazingP (innerGrazingP = 0.2)
                p > innerGrazingP (0.2 < p < 0.5, p = 0.4)
                    o
                p < innerGrazingP (p < 0.2 < 0.5, p = 0.1)
                    o
                    i
                    oio

            if innerGrazingP > outerGrazingP (innerGrazingP = 0.7)
                p < innerGrazingP (p < 0.5, p < 0.7 in all cases)
                    o
                    i
                    oio
        =#

        if (p > outerGrazingP && numInnerPaths == 0) # smooth turning
            println("smooth turning outer paths")
            outerAlpha = getAlpha(p, interface, surface, true)
            return 2 * numOuterPaths * outerAlpha
        else # p < outerGrazingP
            if (p > innerGrazingP)
                println("total internal reflection")
                outerAlpha = getAlpha(p, interface, surface, true)
                return 2 * numOuterPaths * outerAlpha
            end

            if (numInnerPaths > 0) # skip if no inner paths will be drawn
                println("inner paths")
                innerAlpha = getAlpha(p, 0, interface, false)
                #print("innerAlpha: ")
                #println(innerAlpha)
            end
            if (numOuterPaths > 0)
                println("reflecting outer paths")
                outerAlpha = getAlpha(p, interface, surface, true)
                #print("outerAlpha: ")
                #println(outerAlpha)
            end
            return 2 * numInnerPaths * innerAlpha + numOuterPaths * outerAlpha
        end

        #=
        if (numOuterPaths > 0 && p >= PG1) # smooth turning
            println("passing true to isOuterLayer")
            println("in smooth turning")
            #=
            print("p: ")
            println(p)
            print("PG1: ")
            println(PG1)
            =#
            outerAlpha = getAlpha(p, interface, surface, true)
            #print("outerAlpha: ")
            #println(outerAlpha)
            return 2 * numOuterPaths * outerAlpha
        else
            if (numInnerPaths > 0) # skip if no inner paths will be drawn
                println("inner")
                println("passing false to isOuterLayer")
                innerAlpha = getAlpha(p, 0, interface, false)
                print("innerAlpha: ")
                println(innerAlpha)
            end
            println("checking outers")
            outerAlpha = getAlpha(p, interface, surface, true)
            #print("outerAlpha: ")
            #println(outerAlpha)
            return 2 * numInnerPaths * innerAlpha + numOuterPaths * outerAlpha
        end
    =#
    end

# ----- START OF ENCAPSULATED PROGRAM -----

using Roots
R = 0.5
outerGrazingP = 0.5 / cOuter(0.5) # 0.25
innerGrazingP = 0.5 / cInner(0.5) # 0.5
print("outerGrazingP: ")
println(outerGrazingP)
print("innerGrazingP: ")
println(innerGrazingP)

findClosing = true
drawingPaths = true

if (findClosing)
    windingNum = 1    # n
    numPatterns = 15  # m, total number of times to draw pattern
    # windingNum should be coprime with numPatterns * drawPath.length()

    # windingNum 2, 3
    # all o's turning
    # all o's reflecting
    # all i's

    # m * d_p = 2pi * n
    # d_p = 2*a_o + 2*a_i*r
        # r: number of inner paths
    # m : number of full paths
    # n : winding number

    # turning
    # m * d_p = 2pi * n
    # d_p = 2 * outerAlpha * numOuterPaths
    # m * 2 * outerAlpha * numOuterPaths = 2pi * n
    # m * outerAlpha * numOuterPaths = pi * n

    # reflecting
    # m * d_p = 2pi * n
    # d_p = 2 * outerAlpha * numOuterPaths + innerAlpha * numInnerPaths
    # m * (2 * outerAlpha * numOuterPaths + innerAlpha * numInnerPaths) = 2pi * n

    # n = 1, m = 8,  p ~= 0.165, initialGuess 0.2 "oio"
    # n = 1, m = 6,  p ~= 0.210, initialGuess 0.2 "oiio"
    # n = 1, m = 4,  p ~= 0.195, initialGuess 0.2 "oiiio"

    # innerGrazingP = 0.25

    # n = 1, m = 16, p ~= 0.420, initialGuess 0.4 "o"
    # n = 1, m = 8,  p ~= 0.420, initialGuess 0.4 "oo"
    # n = 1, m = 4,  p ~= 0.420, initialGuess 0.4 "oooo"
    
    # outerGrazingP = 0.5

    # n = 2, m = 6,  p ~= 0.858, initialGuess 0.7 "o"
    # n = 2, m = 11, p ~= 0.652, initialGuess 0.7 "o"
    #println("Enter layers path goes through, in order:")
    #drawPath = lowercase(readline())
    drawPath = "i"
    closingPath(p) = numPatterns/(2 * pi) * getEpicentralDistance(p, drawPath) - windingNum

    println("Enter initial guess: ")
    initialGuess = readline()
    initialGuess = parse(Float64, initialGuess)
    # println("calculated p-value: ")
    println(find_zero(closingPath, initialGuess))
end

if (drawingPaths)

    # wave may start in outer or inner layer
    # limit to 4 paths I/O

    # maximum p-value is surface / cOuter(surface)
    println("Enter a value for p:")
    pValue = readline()
    pValue = parse(Float64, pValue)
    #=
    println("Enter layers to draw path in (outer \"O\" or inner \"I\"):")
    println("I's must be separated by an even number of O's.")
    drawLayers = lowercase(readline())
=#
    using Plots
    plotly()
    graphTitle = "p = " * string(pValue)
    plot(title=graphTitle, xlims=(-1, 1), ylims=(-1, 1), aspect_ratio=:equal, size=(1600, 900))
    m1(t) = cos(t)
    n1(t) = sin(t)
    m2(t) = R * cos(t)
    n2(t) = R * sin(t)
    plot!(m1, n1, 0, 2pi, label="Surface", color=:black)
    plot!(m2, n2, 0, 2pi, label="Interface", color=:brown)
    seismicPath(pValue, "i", numPatterns - 1)

    gui()
end