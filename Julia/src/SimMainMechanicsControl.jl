module SimMainMechanicsControl

    using Random
    using FileIO

    include("MathFunctions.jl")
    using .MathFunctions

    include("SharedStructs.jl")
    using .SharedStructs
    include("ReactAdvDiff.jl")
    include("Mechanics.jl")
    include("LightControl.jl")

    function CBProfile(grid, rDomain, CBScale, aOff, flat, n, r0Base, widthBase)
        gamma = LightControl.gammaSoAFuncCyl(grid, r0Base, widthBase, 1.0).Values
        return CBScale .* (1 .- flat) .* gamma .* (aOff .+ (1 .- aOff) .* (rDomain ./ r0Base) .^ n) .+ flat
    end 

    function InitializeAndRun(parameters)
                
        ########################################################################################
        #
        #                    Load the parameters into global variables
        #
        ########################################################################################

        println("Initializing...")

        # General params
        global domainType = parameters["domainType"]
        global Nx = parameters["Nx"] 
        if (domainType == "2D") || (domainType == "rz")
            global Ny = parameters["Ny"]
        else # case cyl
            global Ny = 1
        end
        global dx = parameters["dx"] # set to 1
        global dt = parameters["dt"] # set to 1
        global nSteps = parameters["nSteps"] # number of steps to run
        global timeStride = parameters["timeStride"] # save every timeStride frames
        global startCollecting = parameters["startCollecting"]
        global seed = parameters["seed"] # random seed - will use this if this it is not zero

        ### Boundary conditions
        # BE conditions
        global bcU_X = parameters["bcU_X"] # U condition in x direction
        global bcU_Y = parameters["bcU_Y"] # U condition in x direction

        ### mech params
        mP = Mechanics.MechParams(parameters["mu0"], parameters["lambda0"], parameters["gMin"])

        ########################################################################################
        #
        #               Initialize the various fields using the global variables
        #
        ########################################################################################

        ### Create the grid and the boundary inds
        global grid = Grid2D(Nx, Ny, dx)

        ### Set the random seed if using it
        if seed != 0
            Random.seed!(seed)
        end

        ### RAD init
        global CBISoA = ScalarSoA2D(grid, 0.0)
        global CBASoA = ScalarSoA2D(grid, 0.0)

        global rDomain = collect(1:grid.Nx) .* grid.dx

        if parameters["snapr0"]
            global r0C = parameters["r0"]
        else
            global r0C = parameters["r0Base"]
        end

	@show parameters["snapr0"]
	@show r0C
	
        global CBTot = CBProfile(grid, rDomain, parameters["CBScale"], parameters["aOff"], parameters["flat"], 
            parameters["n"], r0C, parameters["widthBase"])

        ### get light protocol
        iFun = LightControl.getiFun(parameters["cyc"], parameters["len"], parameters["del"], parameters["nSteps"] * parameters["dt"])
        gammaSoABase = LightControl.gammaSoAFuncCyl(grid, parameters["r0"], parameters["width"], 1) 

        gammaSoA = MultiplyScalarSoA2D(grid, gammaSoABase, iFun(0))
        
        CBASoA.Values .= gammaSoA.Values .* CBTot
        CBISoA.Values .= (1 .- gammaSoA.Values) .* CBTot

        ### mech init
        global uxSoA = ScalarSoA2D(grid)
        global uySoA = ScalarSoA2D(grid)

        dF = Mechanics.DispFields(uxSoA, uySoA)

        println("Done initializing.")

        ########################################################################################
        #
        #                            Run the main simulation loop
        #
        ########################################################################################

        # Create arrays to store the data
        global dFArray = []

        println("Beginning simulation...")
            
        # Beginning of loop
        @time for t in 1:nSteps

            # Copy down data to use for the updates of various fields, so that the order of updates doesn't matter
            dFN = deepcopy(dF)

	    if (t%100 == 0)
	    	    GC.gc()
            end

            # Push to the saved data if at multiple of timeStride
            if ((t == 1) || (t % timeStride == 0)) && (t > startCollecting)
	        #GC.gc()
                push!(dFArray, dFN)
            end

            # update gamma
            gammaSoA = MultiplyScalarSoA2D(grid, gammaSoABase, iFun(dt * (t-1)))
            CBASoA.Values .= gammaSoA.Values .* CBTot
            CBISoA.Values .= (1 .- gammaSoA.Values) .* CBTot
            
            Mechanics.PredictorCorrectorStepDispSoA!(grid, domainType, dF, CBASoA, CBISoA, mP, dt, bcU_X, bcU_Y, rDomain)

        end
        println("Done with simulation.")

        ### Return the saved trajectory data, to be selectivey exported from the calling scope
        return (dFArray)

    end # InitializeAndRun


end # Main
