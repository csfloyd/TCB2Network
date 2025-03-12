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

    function thetaFunc(x, y, x0, y0)
        return atan((y-y0) + 1e-6, (x-x0) + 1e-6)
    end

    function Gaussian(x, mean, sig)
        #return (1/(sig * sqrt(2*pi))) * exp(- (x- mean)^2 / (2 * sig^2))
        return exp(- (x- mean)^2 / (2 * sig^2))
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

    function starPolarPattern(theta, scale)
        xmax = 15
        r0 = 0.1 * xmax 
        b = 0.1
        n = 5
        a = 0.81
        return scale * (r0 + sqrt( - log(2 * exp(-a^2) - exp(-(b * xmax * sin((theta - pi/2) * (n/2))) ^ 2) )))
    end

    function starxyPattern(x, y, x0, y0, scale, s)
        r0 = starPolarPattern(thetaFunc(x, y, x0, y0), scale)
        r = rFunc(x, y, x0, y0)
        return spaceFunc(r, r0, s)
    end

    function GetStarSoA(grid, scale, s)
        return ScalarSoA2D([LightControl.starxyPattern(x, y, grid.Nx/2, grid.Nx/2, scale, s) for x in 1:grid.Nx, y in 1:grid.Ny])
    end

    function movingCircle(x, y, Rcent, r0, s)
        r = rFunc(x, y, Rcent[1], Rcent[2])
        return spaceFunc(r, r0, s)
    end

    function circleTimeFunc(R, t, omega)
        return [R * cos(omega * t * 2*pi), -R * sin(omega * t * 2*pi)]
    end 

    function GetMovingCircleSoA(grid, R, r0, s, t, omega)
        Rcent = circleTimeFunc(R, t, omega) .+ [grid.Nx / 2, grid.Ny/2]
        return ScalarSoA2D([movingCircle(x, y, Rcent, r0, s) for x in 1:grid.Nx, y in 1:grid.Ny])
    end
    
    function getiFun(cyc, len, del, tMax)
        trange = 0:0.01:(2*tMax)
        wPoints = [waveFunc(t, cyc, len, del, 2*tMax) for t in trange]
        data = hcat(trange, wPoints)
        return scale(interpolate(wPoints, BSpline(Cubic(Line(OnGrid())))), trange)
    end
    
    function gamFunc(x, y, xc, yc, r0, width, iF)
        return spaceFunc(rFunc(x, y, xc, yc), r0, width) * iF
    end

    function gamFuncCyl(x, y, xc, yc, r0, width, iF)
        return spaceFunc(x, r0, width) * iF
    end

    function gamFuncrz(x, y, xc, yc, r0, width, zWidth, iF)
        return gamFuncCyl(x, y, xc, yc, r0, width, iF) * Gaussian(y, yc, zWidth)
    end

    gammaSoAFunc2D(grid, r0, width, iF) = ScalarSoA2D([gamFunc(x, y, grid.Nx/2, grid.Ny/2, r0, width, iF) for x in 1:grid.Nx, y in 1:grid.Ny])

    gammaSoAFuncCyl(grid, r0, width, iF) = ScalarSoA2D([gamFuncCyl(x, y, grid.Nx/2, grid.Ny/2, r0, width, iF) for x in 1:grid.Nx, y in 1:grid.Ny])
    
    gammaSoAFuncrz(grid, r0, width, zWidth, iF) = ScalarSoA2D([gamFuncrz(x, y, grid.Nx/2, grid.Ny/2, r0, width, zWidth, iF) for x in 1:grid.Nx, y in 1:grid.Ny])
   
    
    
end
