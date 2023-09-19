module LightControl

    include("SharedStructs.jl")
    using .SharedStructs

    include("MathFunctions.jl")
    using .MathFunctions

    using Interpolations
    
    function spaceFunc(r, r0, s)
        return 0.5 * (1 - tanh((r - r0) / s))
    end
    
    function rFunc(x, y, x0, y0)
        return sqrt((x - x0)^2 + (y - y0)^2)
    end
    
    sigmoid(x) = 1 / (1 + exp(-x))
    
    function smoothBump(x, a, xa, xb)
        return sigmoid((x - xa) / a) - sigmoid((x - xb) / a)
    end
    
    function waveFunc(t, cyc, len, del, tMax)
        s = 0.0
        for i in 0:floor(tMax / cyc)
            s += smoothBump(t, del, cyc * i, cyc * i + len)
        end
        return s
    end
    
    
    function getiFun(cyc, len, del, tMax)
        trange = 0:0.01:tMax
        wPoints = [waveFunc(t, cyc, len, del, tMax) for t in trange]
        data = hcat(trange, wPoints)
        return scale(interpolate(wPoints, BSpline(Cubic(Line(OnGrid())))), trange)
    end
    
    function gamFunc(x, y, xc, yc, r0, width, iF)
        return spaceFunc(rFunc(x, y, xc, yc), r0, width) * iF
    end

    function gamFuncCyl(x, y, xc, yc, r0, width, iF)
        return spaceFunc(x, r0, width) * iF
    end

    gammaSoAFunc2D(grid, r0, width, iF) = ScalarSoA2D([gamFunc(x, y, grid.Nx/2, grid.Nx/2, r0, width, iF) for x in 1:grid.Nx, y in 1:grid.Ny])

    gammaSoAFuncCyl(grid, r0, width, iF) = ScalarSoA2D([gamFuncCyl(x, y, grid.Nx/2, grid.Nx/2, r0, width, iF) for x in 1:grid.Nx, y in 1:grid.Ny])
    

    
end
