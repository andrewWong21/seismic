function cOuter(r)
    return r^-2
end

function cInner(r) 
    return 0.5 * r^-1
end

function cCore(r)
    return 0.5 * r^-1
end

function cOuterSecondary(r)
    return 0.5 * cOuter(r)
end

function cInnerSecondary(r)
    return 0.5 * cInner(r)
end

function cCoreSecondary(r)
    return 0.5 * cCore(r)
end

function f(p, r, wave, layer)
    valC = 0
    if (wave == 'p' && layer == 'o')
        valC = cOuter(r)
    elseif (wave == 'p' && layer == 'i')
        valC = cInner(r)
    elseif (wave == 'p' && layer == 'k')
        valC = cCore(r)
    elseif (wave == 's' && layer == 'o')
        valC = cOuterSecondary(r)
    elseif (wave == 's' && layer == 'i')
        valC = cInnerSecondary(r)
    elseif (wave == 's' && layer == 'k')
        valC = cCoreSecondary(r)
    else
        println("ERROR: Invalid inputs.")
        print("Current Wave Type: ")
        println(wave)
        print("Current Layer: ")
        println(layer)
    end
    return p / (r * sqrt((r / valC)^2 - p^2) )
end

function midpoint(a, b, p, wave, layer)
    sum = 0
    h = (b - a) / 30 # step size
    for num in 0:29
        sum += f(p, a + num * h + h / 2, wave, layer)
    end
    return sum * h
end

function checkBoundary(p, innerBoundary, wave, layer)
    # determines if path is turning or reflecting off inner boundary
    # grazingP = R / c(R)
    # if p > grazingP, path is turning smoothly
    # if p < grazingP, path is reflecting off inner boundary

    # two layers, two grazing P values outerGrazingP and innerGrazingP, innerGrazingP < outerGrazingP
    # if p >= outerGrazingP, turning in outer layer (assuming it starts at outer layer)
    # if innerGrazingP < p < outerGrazingP, total internal reflection (does not enter inner layer)
    # if p < innerGrazingP, transmission from outer layer to inner layer
    # if p = innerGrazingP, headwave: critical angle, glides along interface

    # threshold/grazing p-value
    if (wave == 'p' && layer == 'o')
        grazingP = innerBoundary / cOuter(innerBoundary)
    elseif (wave == 'p' && layer == 'i')
        grazingP = innerBoundary / cInner(innerBoundary)
    elseif (wave == 's' && layer == 'o')
        grazingP = innerBoundary / cOuterSecondary(innerBoundary)
    elseif (wave == 's' && layer == 'i')
        grazingP = innerBoundary / cInnerSecondary(innerBoundary)
    elseif (layer == 'k')
        return p
    else
        println("ERROR: Invalid inputs.")
        print("Current Wave Type: ")
        println(wave)
        print("Current Layer: ")
        println(layer)
    end
    #= 
    print("p: ")
    println(p)
    print("grazingP: ")
    println(grazingP) =#
    # for reflections
    if (p < grazingP)
        # println("p < grazingP, reflect")
        return grazingP
    elseif (p > grazingP)
        # println("p > grazingP, smooth turn")
        return p
    else
        # println("p = grazingP")
        return p
    end
end

function getRadius(p, wave, layer)
    # minimum turning radius is given by the grazing p value
    if (wave == 'p' && layer == 'o')
        grazeOuter(r) = r / (cOuter(r)) - p
        turningRadius = find_zero(grazeOuter, 0.5)
    elseif (wave == 'p' && layer == 'i')
        grazeInner(r) = r / (cInner(r)) - p
        turningRadius = find_zero(grazeInner, 0.5)
    elseif (wave == 'p' && layer == 'k')
        grazeCore(r) = r / (cCore(r)) - p
        turningRadius = find_zero(grazeCore, 0.5)
    elseif (wave == 's' && layer == 'o')
        grazeOuterSecondary(r) = r / (cOuterSecondary(r)) - p
        turningRadius = find_zero(grazeOuterSecondary, 0.5)
    elseif (wave == 's' && layer == 'i')
        grazeInnerSecondary(r) = r / (cInnerSecondary(r)) - p
        turningRadius = find_zero(grazeInnerSecondary, 0.5)
    elseif (wave == 's' && layer == 'k')
        grazeCoreSecondary(r) = r / (cCoreSecondary(r)) - p
        turningRadius = find_zero(grazeCoreSecondary, 0.5)
    else
        println("ERROR: Invalid inputs.")
        print("Current Wave Type: ")
        println(wave)
        print("Current Layer: ")
        println(layer)
    end
    return turningRadius
end

function getAlpha(p, innerBoundary, outerRadius, wave, layer)
    minP = checkBoundary(p, innerBoundary, wave, layer)
    innerRadius = getRadius(minP, wave, layer) # lower bound of integral
    return midpoint(innerRadius, outerRadius, p, wave, layer)
end

function fillRadiusValues(p, minRadius, maxRadius)
    h = (maxRadius - minRadius) / 30
    radius_values = collect(minRadius:h:maxRadius)
    radius_path = vcat(reverse(radius_values), radius_values[2:31])
    return radius_path
end

function fillThetaValues(p, alpha, radius_values, minRadius, wave, layer, offset)
    # turningRadius = getRadius(p, true)
    theta_values = zeros(Float64, 31)
    for num in 1:30
        setindex!(theta_values, midpoint(minRadius, getindex(radius_values, num), p, wave, layer) + (pi / 2 - alpha), num)
    end
    setindex!(theta_values, (pi / 2 - alpha), 31) # pi/2
    
    theta_path = vcat((theta_values), zeros(Float64, 30))
    for num in 2:31
        setindex!(theta_path, -1 * getindex(reverse(theta_values), num) + 2 * (pi / 2 - alpha), num + 30)
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

function fillXHalf(radius_path, theta_path)
    x_values = zeros(Float64, 31)
    for num in 1:31
        setindex!(x_values, getindex(radius_path, num) * cos(getindex(theta_path, num)), num)
    end
    return x_values
end

function fillYHalf(radius_path, theta_path)
    y_values = zeros(Float64, 31)
    for num in 1:31
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

# -------------------------------------------------------------------------------------------------
function drawSegment(p, wave, layer, offset, parity) # currently unused
    if (layer == 'o') # outer layer
        if (wave == 'p' && p > outerGrazingP) # turning p-wave
            println()
        elseif (wave == 'p' && p < outerGrazingP) # reflecting p-wave
            if (parity % 2 == 0)
                println()
            else
                println()
            end
        elseif (wave == 's' && p > secondaryOuterGrazingP) # turning s-wave
            println()
        end
    elseif (layer == 'i') # inner layer
        if (wave == 'p' && p > innerGrazingP) # turning p-wave
            println()
        elseif (wave == 'p' && p < innerGrazingP) # reflecting p-wave
            println()
        elseif (wave == 's' && p > secondaryinnerGrazingP) # turning s-wave
            println()
        end
    elseif (layer == 'k')
        println("Not implemented yet.")
    end
end

function seismicPath(p, drawLayers, numTotalPaths)
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

    surface = 1
    outerInterface = firstInterface
    innerInterface = secondInterface
    fullPathRadiusVals = Float64[]
    fullPathThetaVals = Float64[]
    outerOffset = 0
    innerOffset = 0
    coreOffset = 0

    outerParity = 0
    innerParity = 0
    
    # TODO: fix offset calculations
    #= 
    if (cmp(drawLayers[1], 'o') == 0) # first path in outer layer
        outerOffset = 0
        outerParity = 0
    end
    elseif  (cmp(drawLayers[1], 'o') == 0) # first path in inner layer, measures for fixing exiting half of outer layer path
        outerOffset = -1 * getAlpha(p, outerInterface, surface, true)
        outerParity = 1
    else # first path in core layer
        outerOffset = 
        innerOffset = 
        outerParity = 
        innerParity = 
    end =#

    for currentLayer in drawLayers
        alpha = 0
        radius_values = Array{Float64}(undef, 61)
        theta_values = Array{Float64}(undef, 61)
        x_values = Array{Float64}(undef, 61)
        y_values = Array{Float64}(undef, 61)

        if (currentLayer == 'o')
            outerParity += 1
            println("\nDrawing in outer layer")
            
            # println("getting minRadius")
            pValue = checkBoundary(p, outerInterface, 'o')
            minRadius = getRadius(pValue, 'o')
            maxRadius = surface

            # println("getting alpha")
            alpha = getAlpha(p, outerInterface, maxRadius, 'o')
            radius_values = fillRadiusValues(p, minRadius, maxRadius)
            theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 'o', outerOffset)
        elseif (currentLayer == 'i') # inner layer
            innerParity += 1
            println("\nDrawing in inner layer")

            # println("getting minRadius")
            pValue = checkBoundary(p, innerInterface, 'i')
            minRadius = getRadius(pValue, 'i')
            # print("minRadius: ")
            # println(minRadius)
            maxRadius = outerInterface

            # println("getting alpha")
            alpha = getAlpha(p, minRadius, maxRadius, 'i')
            radius_values = fillRadiusValues(p, minRadius, maxRadius)
            theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 'i', innerOffset)
        elseif (currentLayer == 'k') # core layer
            println("\nDrawing in core layer")

            # println("getting minRadius")
            minRadius = getRadius(p, 'k') # 0
            # print("minRadius: ")
            # println(minRadius)
            maxRadius = innerInterface
            # print("innerInterface: ")
            # println(innerInterface)

            # println("getting alpha")
            alpha = getAlpha(p, 0, maxRadius, 'k')
            radius_values = fillRadiusValues(p, minRadius, maxRadius)
            theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 'k', coreOffset)
        end

        x_values = fillXValues(radius_values, theta_values)
        y_values = fillYValues(radius_values, theta_values)

        # path plotting
        if (currentLayer == 'o' && p > outerGrazingP)
            println("outer turning")
            fullPathRadiusVals = [fullPathRadiusVals; radius_values]
            fullPathThetaVals = [fullPathThetaVals; theta_values]
            plot!(x_values, y_values, label="outer", color=:brown)
            outerOffset += alpha * 2
            innerOffset += alpha * 2
            coreOffset  += alpha * 2
        end
        if (currentLayer == 'o' && p < outerGrazingP)
            println("outer reflecting")

            if (outerParity % 2 == 1) # only drawing first half of path
                println("outer, first half")
                fullPathRadiusVals = [fullPathRadiusVals; radius_values[1:31]]
                fullPathThetaVals = [fullPathThetaVals; theta_values[1:31]]
                plot!(x_values[1:31], y_values[1:31], label="outer, first half", color=:green)
                innerOffset += alpha
                coreOffset  += alpha
            elseif (outerParity % 2 == 0) # only drawing second half of path
                println("outer, second half")
                fullPathRadiusVals = [fullPathRadiusVals; radius_values[31:61]]
                fullPathThetaVals = [fullPathThetaVals; theta_values[31:61]]
                plot!(x_values[31:61], y_values[31:61], label="outer, second half", color=:red)
                outerOffset += alpha * 2
                innerOffset += alpha
                coreOffset  += alpha
            end
        end
        if (currentLayer == 'i' && p > innerGrazingP)
            println("inner turning")
            fullPathRadiusVals = [fullPathRadiusVals; radius_values]
            fullPathThetaVals = [fullPathThetaVals; theta_values]
            plot!(x_values, y_values, label="inner", color=:darkorange)
            outerOffset += alpha * 2
            innerOffset += alpha * 2
            coreOffset  += alpha * 2
        end
        if (currentLayer == 'i' && p < innerGrazingP)
            println("inner reflecting")

            if (innerParity % 2 == 1) # only drawing first half of path
                println("inner, first half")
                fullPathRadiusVals = [fullPathRadiusVals; radius_values[1:31]]
                fullPathThetaVals = [fullPathThetaVals; theta_values[1:31]]
                plot!(x_values[1:31], y_values[1:31], label="inner, first half", color=:cyan)
                outerOffset += alpha
                coreOffset  += alpha
            elseif (innerParity % 2 == 0) # only drawing second half of path
                println("inner, second half")
                fullPathRadiusVals = [fullPathRadiusVals; radius_values[31:61]]
                fullPathThetaVals = [fullPathThetaVals; theta_values[31:61]]
                plot!(x_values[31:61], y_values[31:61], label="inner, second half", color=:orange)
                outerOffset += alpha
                innerOffset += alpha * 2
                coreOffset  += alpha
            end
        end
        if (currentLayer == 'k')
            println("core path")
            fullPathRadiusVals = [fullPathRadiusVals; radius_values]
            fullPathThetaVals = [fullPathThetaVals; theta_values]
            plot!(x_values, y_values, label="core", color=:purple)
            outerOffset += alpha * 2
            innerOffset += alpha * 2
            coreOffset +=  alpha * 2
        end

    end # for loop

    reflectTotalPath(fullPathRadiusVals, fullPathThetaVals, getEpicentralDistance(p, drawLayers), numTotalPaths - 1)
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

function getEpicentralDistance(p, drawLayers)
    
    # determine total offset for one full path (epicentral distance)
    # m * d_p1 = 2pi * n
    # d_p1 = 2*a_o + 2*a_i*r
        # r: number of inner paths
    # m : number of full paths
    # n : winding number
    # use with findZeroes to get p value that results in a closing path

    numOuterPaths = 0
    numInnerPaths = 0
    numCorePaths = 0
    for layer in drawLayers
        if (layer == 'o')
            numOuterPaths += 1
        elseif (layer == 'i')
            numInnerPaths += 1
        elseif (layer == 'k')
            numCorePaths += 1
        end
    end

    outerInterface = firstInterface
    innerInterface = secondInterface
    surface = 1

    outerAlpha = 0
    innerAlpha = 0
    coreAlpha = 0

    distance = 0
    if (numOuterPaths > 0)
        if (p > outerGrazingP)
            println("smooth turning outer paths")
            outerAlpha = getAlpha(p, outerInterface, surface, 'o')
            distance += 2 * numOuterPaths * outerAlpha
        else
            println("reflecting outer paths")
            outerAlpha = getAlpha(p, outerInterface, surface, 'o')
            distance += numOuterPaths * outerAlpha
        end
    end
    if (numInnerPaths > 0)
        if (p > innerGrazingP)
            println("smooth turning inner paths")
            innerAlpha = getAlpha(p, innerInterface, outerInterface, 'i')
            distance += 2 * numInnerPaths * innerAlpha
        else
            println("reflecting inner paths")
            innerAlpha = getAlpha(p, innerInterface, outerInterface, 'i')
            distance += numInnerPaths * innerAlpha
        end
    end
    if (numCorePaths > 0)
        coreAlpha = getAlpha(p, 0, innerInterface, 'k')
        distance += 2 * numCorePaths * coreAlpha
    end
    return distance
end

function seismicWave(p, waveTypes, layerTypes, numTotalPaths)

    # p > primaryGrazingP, p > secondaryGrazingP, both turning
    # p > primaryGrazingP, p < secondaryGrazingP, P turning, S reflecting

    # p < primaryGrazingP, p < secondaryGrazingP, both reflecting
    # PPSS - two full reflections
    # PS - first half of P, second half of S
    # SP - first half of S, second half of P

    surface = 1
    outerInterface = firstInterface
    innerInterface = secondInterface
    fullPathRadiusVals = Float64[]
    fullPathThetaVals = Float64[]
    outerOffset = 0
    outerParity = 0
    innerOffset = 0
    innerParity = 0
    coreOffset = 0

    for num in 1:length(waveTypes)
        currentWaveType = getindex(waveTypes, num)
        currentLayer = getindex(layerTypes, num)

        alpha = 0
        radius_values = Array{Float64}(undef, 61)
        theta_values = Array{Float64}(undef, 61)
        x_values = Array{Float64}(undef, 61)
        y_values = Array{Float64}(undef, 61)

        if (currentWaveType == 'p')
            print("Drawing primary wave in ")
            if (currentLayer == 'o')
                println("outer layer")
                pValue = checkBoundary(p, outerInterface, 'p', 'o')
                minRadius = getRadius(pValue, 'p', 'o')
                maxRadius = surface

                alpha = getAlpha(p, outerInterface, maxRadius, 'p', 'o')
                radius_values = fillRadiusValues(p, minRadius, maxRadius)
                theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 'p', 'o', outerOffset)
            elseif (currentLayer == 'i')
                println("inner layer")
                pValue = checkBoundary(p, innerInterface, 'p', 'i')
                minRadius = getRadius(pValue, 'p', 'i')
                maxRadius = outerInterface

                alpha = getAlpha(p, innerInterface, maxRadius, 'p', 'i')
                radius_values = fillRadiusValues(p, minRadius, maxRadius)
                theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 'p', 'i', innerOffset)
            elseif (currentLayer == 'k')
                println("core layer")
                pValue = checkBoundary(p, 0, 'p', 'k')
                minRadius = getRadius(pValue, 'p', 'k') # 0
                maxRadius = innerInterface

                alpha = getAlpha(p, innerInterface, maxRadius, 'p', 'k')
                radius_values = fillRadiusValues(p, minRadius, maxRadius)
                theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 'p', 'k', coreOffset)
            end
        end
        if (currentWaveType == 's')
            print("\nDrawing secondary wave in ")
            #=
            pValue = checkBoundary(p, outerInterface, 's')
            minRadius = getRadius(pValue, 's')
            maxRadius = surface

            alpha = getAlpha(p, outerInterface, maxRadius, 's')
            radius_values = fillRadiusValues(p, minRadius, maxRadius)
            theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 's', outerOffset) =#
            if (currentLayer == 'o')
                println("outer layer")
                pValue = checkBoundary(p, outerInterface, 's', 'o')
                minRadius = getRadius(pValue, 's', 'o')
                maxRadius = surface

                alpha = getAlpha(p, outerInterface, maxRadius, 's', 'o')
                radius_values = fillRadiusValues(p, minRadius, maxRadius)
                theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 's', 'o', outerOffset)
            elseif (currentLayer == 'i')
                println("inner layer")
                pValue = checkBoundary(p, innerInterface, 's', 'i')
                minRadius = getRadius(pValue, 's', 'i')
                maxRadius = outerInterface

                alpha = getAlpha(p, innerInterface, maxRadius, 's', 'i')
                radius_values = fillRadiusValues(p, minRadius, maxRadius)
                theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 's', 'i', innerOffset)
            elseif (currentLayer == 'k')
                println("core layer")
                pValue = checkBoundary(p, 0, 's', 'k')
                minRadius = getRadius(pValue, 's', 'k') # 0
                maxRadius = innerInterface

                alpha = getAlpha(p, innerInterface, maxRadius, 'p', 'k')
                radius_values = fillRadiusValues(p, minRadius, maxRadius)
                theta_values = fillThetaValues(p, alpha, radius_values, minRadius, 's', 'k', coreOffset)
            else
                println("ERROR: Invalid layer.")
            end
        end
        print("currentLayer: ")
        println(currentLayer)
        # path plotting
        if (currentLayer == 'o')
            if (currentWaveType == 'p' && p >= outerGrazingP)
                println("primary outer turning")
                fullPathRadiusVals = [fullPathRadiusVals; radius_values]
                fullPathThetaVals = [fullPathThetaVals; theta_values]

                x_values = fillXValues(radius_values, theta_values)
                y_values = fillYValues(radius_values, theta_values)

                plot!(x_values, y_values, label="primary", color=:red)
                outerOffset += alpha * 2
                innerOffset += alpha * 2
                coreOffset += alpha * 2
            elseif (currentWaveType == 's' && p >= secondaryOuterGrazingP)
                println("secondary outer turning")
                fullPathRadiusVals = [fullPathRadiusVals; radius_values]
                fullPathThetaVals = [fullPathThetaVals; theta_values]

                x_values = fillXValues(radius_values, theta_values)
                y_values = fillYValues(radius_values, theta_values)

                plot!(x_values, y_values, label="secondary", color=:orange)
                outerOffset += alpha * 2
                innerOffset += alpha * 2
                coreOffset += alpha * 2
            else
                println("reflecting")
                if (outerParity % 2 == 0) # only drawing first half of path
                    println("outer, first half")
                    fullPathRadiusVals = [fullPathRadiusVals; radius_values[1:31]]
                    fullPathThetaVals = [fullPathThetaVals; theta_values[1:31]]

                    x_values = fillXHalf(radius_values[1:31], theta_values[1:31])
                    y_values = fillYHalf(radius_values[1:31], theta_values[1:31])

                    if (currentWaveType == 'p')
                        plot!(x_values, y_values, label="1st half (p)", color=:red)
                    else
                        plot!(x_values, y_values, label="1st half (s)", color=:orange)
                    end

                    outerOffset += alpha
                    innerOffset += alpha
                    coreOffset += alpha
                elseif (outerParity % 2 == 1) # only drawing second half of path
                    println("outer, second half")

                    secondHalfThetaVals = theta_values[31:61]
                    for num in 1:31
                        setindex!(secondHalfThetaVals, getindex(secondHalfThetaVals, num) + alpha, num)
                    end

                    fullPathRadiusVals = [fullPathRadiusVals; radius_values[31:61]]
                    fullPathThetaVals = [fullPathThetaVals; secondHalfThetaVals]

                    x_values = fillXHalf(radius_values[31:61], secondHalfThetaVals)
                    y_values = fillYHalf(radius_values[31:61], secondHalfThetaVals)

                    if (currentWaveType == 'p')
                        plot!(x_values, y_values, label="2nd half (p)", color=:red)
                    else
                        plot!(x_values, y_values, label="2nd half (s)", color=:orange)
                        end
                    outerOffset += alpha
                    innerOffset += alpha
                    coreOffset += alpha
                end
                outerParity += 1
            end
        end
        if (currentLayer == 'i')
            if (currentWaveType == 'p' && p >= innerGrazingP)
                println("primary inner turning")

                fullPathRadiusVals = [fullPathRadiusVals; radius_values]
                fullPathThetaVals = [fullPathThetaVals; theta_values]

                x_values = fillXValues(radius_values, theta_values)
                y_values = fillYValues(radius_values, theta_values)

                plot!(x_values, y_values, label="primary", color=:red)
                outerOffset += alpha * 2
                innerOffset += alpha * 2
            elseif (currentWaveType == 's' && p >= secondaryInnerGrazingP)
                println("secondary inner turning")

                fullPathRadiusVals = [fullPathRadiusVals; radius_values]
                fullPathThetaVals = [fullPathThetaVals; theta_values]

                x_values = fillXValues(radius_values, theta_values)
                y_values = fillYValues(radius_values, theta_values)

                plot!(x_values, y_values, label="secondary", color=:orange)
                outerOffset += alpha * 2
                innerOffset += alpha * 2
            else
                println("reflecting")
                if (innerParity % 2 == 0) # only drawing first half of path
                    println("inner, first half")

                    fullPathRadiusVals = [fullPathRadiusVals; radius_values[1:31]]
                    fullPathThetaVals = [fullPathThetaVals; theta_values[1:31]]

                    x_values = fillXHalf(radius_values[1:31], theta_values[1:31])
                    y_values = fillYHalf(radius_values[1:31], theta_values[1:31])

                    if (currentWaveType == 'p')
                        plot!(x_values, y_values, label="1st half (p)", color=:red)
                    else
                        plot!(x_values, y_values, label="1st half (s)", color=:orange)
                    end
                elseif (innerParity % 2 == 1) # only drawing second half of path
                    println("inner, second half")

                    secondHalfThetaVals = theta_values[31:61]
                    for num in 1:31
                        setindex!(secondHalfThetaVals, getindex(secondHalfThetaVals, num) + alpha, num)
                    end

                    fullPathRadiusVals = [fullPathRadiusVals; radius_values[31:61]]
                    fullPathThetaVals = [fullPathThetaVals; secondHalfThetaVals]

                    x_values = fillXHalf(radius_values[31:61], secondHalfThetaVals)
                    y_values = fillYHalf(radius_values[31:61], secondHalfThetaVals)

                    if (currentWaveType == 'p')
                        plot!(x_values, y_values, label="2nd half (p)", color=:red)
                    else
                        plot!(x_values, y_values, label="2nd half (s)", color=:orange)
                    end
                end
                outerOffset += alpha
                innerOffset += alpha
                innerParity += 1
            end
        end
        if (currentLayer == 'k')
            # TODO
        end
    # reflectTotalPath(fullPathRadiusVals, fullPathThetaVals, getWaveDistance(p, waveTypes), numTotalPaths - 1)
    end
end

function getWaveDistance(p, waveTypes)
    numPrimaryWaves = 0
    numSecondaryWaves = 0
    for wave in waveTypes
        if (wave == 'p')
            numPrimaryWaves += 1
        elseif (wave == 's')
            numSecondaryWaves += 1
        end
    end

    outerInterface = firstInterface
    innerInterface = secondInterface
    surface = 1

    primaryAlpha = 0
    secondaryAlpha = 0

    distance = 0
    if (numPrimaryWaves > 0)
        primaryAlpha = getAlpha(p, outerInterface, surface, 'p')
        if (p > outerGrazingP)
            println("smooth turning primary waves")
            distance += 2 * numPrimaryWaves * primaryAlpha
        else
            println("reflecting primary waves")
            distance += numPrimaryWaves * primaryAlpha
        end
    end
    if (numSecondaryWaves > 0)
        secondaryAlpha = getAlpha(p, outerInterface, surface, 's')
        if (p > secondaryGrazingP)
            println("smooth turning secondary waves")
            distance += 2 * numSecondaryWaves * secondaryAlpha
        else
            println("reflecting secondary waves")
            distance += numSecondaryWaves * secondaryAlpha
        end
    end
    return distance
end

# ------------------------------

using Roots
firstInterface = 0.667
secondInterface = 0.333

# maximum p-value is outerbound / c(outerbound)
outerMaxP     =               1 / cOuter(1)
outerGrazingP =  firstInterface / cOuter(firstInterface)
innerMaxP     =  firstInterface / cInner(firstInterface)
innerGrazingP = secondInterface / cInner(secondInterface)
coreMaxP      = secondInterface / cCore(secondInterface)

secondaryOuterMaxP     =               1 / cOuterSecondary(1)
secondaryOuterGrazingP =  firstInterface / cOuterSecondary(firstInterface)
secondaryInnerMaxP     =  firstInterface / cInnerSecondary(firstInterface)
secondaryInnerGrazingP = secondInterface / cInnerSecondary(secondInterface)
secondaryCoreMaxP      = secondInterface / cCoreSecondary(secondInterface)

findClosingPaths = false
drawingPaths = false
findClosingWaves = false
drawingWaves = true

if (findClosingPaths)
    windingNum = 2    # n
    numPatterns = 33   # m, total number of times to draw pattern
    # windingNum should be coprime with numPatterns * drawPath.length()

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

    #= 
    println("Enter a value for n (number of full circles made):")
    windingNum = readline()
    windingNum = parse(Float64, pValue)

    println("Enter a value for m (number of paths needed to make n circles):")
    numPatterns = readline()
    numPatterns = parse(Float64, pValue) =#

    # p-wave/s-wave mode conversion

    # grazingNumReflect = (pi * m) / grazingAlpha

    # println("Enter layers path goes through, in order:")
    # println("note: if reflecting paths within a single layer are desired, layers must be in pairs (""")
    # drawPath = lowercase(readline())
    drawPath = "oo"
    closingPath(p) = numPatterns / (2 * pi) * getEpicentralDistance(p, drawPath) - windingNum

    println("Enter initial guess: ")
    initialGuess = readline()
    initialGuess = parse(Float64, initialGuess)
    # println("calculated p-value: ")
    println(find_zero(closingPath, initialGuess))
end

if (drawingPaths)
    
    print(outerGrazingP)
    print(" < outerP < ")
    println(outerMaxP)
    print(innerGrazingP)
    print(" < innerP < ")
    println(innerMaxP)
    print("0 < coreP < ")
    println(coreMaxP)
    # wave may start in outer, inner, or core layer
    println("Enter layers to draw path in (outer \"O\", inner \"I\", or core \"K\"):\n")
    println("O's must not be adjacent to K's.")
    println("I's must be separated by an even number of O's.")
    println("K's must be separated by an even number of I's.")
    println("Acceptable examples: oooooo ooiioo oikkio ikiiki\n")
    drawLayers = lowercase(readline())

    println("Enter a value for p:")
    pValue = readline()
    pValue = parse(Float64, pValue)

    # oi p < outerGrazingP && p < innerMaxP
    # ik p < innerGrazingP && p < coreMaxP
    # if occursin("oi", drawLayers)

    using Plots
    plotly()
    # graphTitle = "p = " * string(pValue) * ", path " * drawLayers
    # plot(title=graphTitle, xlims=(-1, 1), ylims=(-1, 1), aspect_ratio=:equal, size=(1600, 900))
    plot(xlims=(-1, 1), ylims=(-1, 1), aspect_ratio=:equal, size=(1600, 900))
    m1(t) = cos(t)
    n1(t) = sin(t)
    plot!(m1, n1, 0, 2pi, label="Surface", color=:black)
    m2(t) = firstInterface * cos(t)
    n2(t) = firstInterface * sin(t)
    plot!(m2, n2, 0, 2pi, label="Outer Interface", color=:brown)
    m3(t) = secondInterface * cos(t)
    n3(t) = secondInterface * sin(t)
    plot!(m3, n3, 0, 2pi, label="Inner Interface", color=:tan)

    seismicPath(pValue, drawLayers, numPatterns - 1)

    gui()
end

if (findClosingWaves)
    windingNum = 1    # n
    numPatterns = 5   # m, total number of times to draw pattern
    # windingNum should be coprime with numPatterns * drawPath.length()

    # p-wave/s-wave mode conversion

    # grazingNumReflect = (pi * m) / grazingAlpha

    pathTypes = "ps"
    closingWave(p) = numPatterns / (2 * pi) * getWaveDistance(p, pathTypes) - windingNum

    println("Enter initial guess: ")
    initialGuess = readline()
    initialGuess = parse(Float64, initialGuess)
    # println("calculated p-value: ")
    println(find_zero(closingWave, initialGuess))
    println()
end

if (drawingWaves)
    
    print(outerGrazingP)
    print(" < primaryP < ")
    println(outerMaxP)
    print(secondaryOuterGrazingP)
    print(" < secondaryP < ")
    println(secondaryOuterMaxP)
    
    print(innerGrazingP)
    print(" < primaryP < ")
    println(innerMaxP)
    print(secondaryInnerGrazingP)
    print(" < secondaryP < ")
    println(secondaryInnerMaxP)
    
    println("Enter a value for p:")
    pValue = readline()
    pValue = parse(Float64, pValue)

    waveTypes = "sspps"
    layerTypes = "iiiii"


    println("Enter types of waves to draw (primary \"P\" or secondary \"S\"):\n")
    # waveTypes = lowercase(readline())

    println("Enter layers that each wave travels in (outer \"O\", inner \"I\", or core \"K\"):\n")
    print(waveTypes)
    println(" - the wave types that were previously entered. Make sure the number of layers match the number of wavetypes.")
    # waveTypes = lowercase(readline())

    # primary wave in outer layer transmitting to inner layer, converting to secondary
    # waveTypes ps, layerTypes oi
    
    using Plots
    plotly()
    # graphTitle = "p = " * string(pValue) * ", path " * drawLayers
    # plot(title=graphTitle, xlims=(-1, 1), ylims=(-1, 1), aspect_ratio=:equal, size=(1600, 900))
    plot(xlims=(-1, 1), ylims=(-1, 1), aspect_ratio=:equal, size=(1600, 900))
    m1(t) = cos(t)
    n1(t) = sin(t)
    plot!(m1, n1, 0, 2pi, label="Surface", color=:black)
    m2(t) = firstInterface * cos(t)
    n2(t) = firstInterface * sin(t)
    plot!(m2, n2, 0, 2pi, label="Outer Interface", color=:brown)
    m3(t) = secondInterface * cos(t)
    n3(t) = secondInterface * sin(t)
    plot!(m3, n3, 0, 2pi, label="Inner Interface", color=:gray)

    # seismicWave(pValue, waveTypes, layerTypes, numPatterns)
    seismicWave(pValue, waveTypes, layerTypes, 1)
    #= 
    goal: draw the first half of a reflecting primary wave in the outer layer, then a turning secondary wave in the inner layer
    input: "ps", "oi"
    what the program will do:
    -calculate minradius, maxradius, theta, and alpha for outer primary
    -use set of radius and theta values to calculate x- and y-values
    -plot x and y

    consequences of using drawSegment: 
    -parity must be updated in seismicWave outside the function's scope, after it is completed
    -can plot x and y within drawSegment with plot!(...)
    -cannot pass alpha, x-set and y-set back to seismicWave
    -alpha is determined outside the function's scope, before drawSegment is called =#
    gui()
end
