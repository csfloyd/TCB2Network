module SimMain

    using Random
    using FileIO

    include("MathFunctions.jl")
    using .MathFunctions

    include("SharedStructs.jl")
    using .SharedStructs
    include("ReactAdvDiff.jl")
    include("Mechanics.jl")
    include("LightControl.jl")

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
        global chemOnly = parameters["chemOnly"]

        ### Boundary conditions
        # BE conditions
        global bcU_X = parameters["bcU_X"] # U condition in x direction
        global bcU_Y = parameters["bcU_Y"] # U condition in x direction
        # RAD conditions
        global bcRAD_X = parameters["bcRAD_X"] # RAD condition in x direction
        global bcRAD_Y = parameters["bcRAD_Y"] # RAD condition in y direction

        ### chem params
        cP = ReactAdvDiff.ChemParams(parameters["kAct"], parameters["kInact"], parameters["kTrap"], parameters["kRel"], 
            parameters["kIBind"], parameters["kIUnbind"], parameters["kABind"], parameters["kAUnbind"], 
            parameters["beta"], parameters["CSat"], parameters["DT"], parameters["DD"], parameters["DC"])

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
        global CDISoA = ScalarSoA2D(grid, Float64(parameters["CDI0"]))
        global CDASoA = ScalarSoA2D(grid, Float64(parameters["CDA0"]))
        global CBISoA = ScalarSoA2D(grid, Float64(parameters["CBI0"]))
        global CBASoA = ScalarSoA2D(grid, Float64(parameters["CBA0"]))
        global CCSoA = ScalarSoA2D(grid, Float64(parameters["CC0"]))
        global CDSoA = ScalarSoA2D(grid, Float64(parameters["CD0"]))
        global CDstSoA = ScalarSoA2D(grid, Float64(parameters["CDst0"]))
        global rDomain = collect(1:grid.Nx) .* grid.dx

        cF = ReactAdvDiff.ConcFields(CDISoA, CDASoA, CBISoA, CBASoA, CCSoA, CDSoA, CDstSoA)

        ### mech init
        global uxSoA = ScalarSoA2D(grid)
        global uySoA = ScalarSoA2D(grid)

        dF = Mechanics.DispFields(uxSoA, uySoA)

        ### get light protocol
	    lenFac = parameters["cyc"] / 30.0
	    len = parameters["len"] * lenFac
        #iFun = LightControl.getiFun(parameters["cyc"], parameters["len"], parameters["del"], parameters["nSteps"] * parameters["dt"])
	    iFun = LightControl.getiFun(parameters["cyc"] + len, len, parameters["del"], parameters["nSteps"] * parameters["dt"])
        
        if parameters["domainType"] == "cyl"
            gammaSoABase = LightControl.gammaSoAFuncCyl(grid, parameters["r0"], parameters["width"], 1)
        elseif parameters["domainType"] == "rz"
            gammaSoABase = LightControl.gammaSoAFuncrz(grid, parameters["r0"], parameters["width"], parameters["zWidth"], 1)
        else
            if parameters["starBool"]
                gammaSoABase = LightControl.GetStarSoA(grid, parameters["r0"], parameters["width"])
            else
                gammaSoABase = LightControl.gammaSoAFunc2D(grid, parameters["r0"], parameters["width"], 1)
            end
        end

        println("Done initializing.")

        ########################################################################################
        #
        #                            Run the main simulation loop
        #
        ########################################################################################

        # Create arrays to store the data
        global cFArray = [] 
        global dFArray = []

        println("Beginning simulation...")
            
        # Beginning of loop
        @time for t in 1:nSteps

            # Copy down data to use for the updates of various fields, so that the order of updates doesn't matter
            cFN = deepcopy(cF)
            dFN = deepcopy(dF)

	    if (t%100 == 0)
	    	    GC.gc()
            end

            # Push to the saved data if at multiple of timeStride
            if ((t == 1) || (t % timeStride == 0)) && (t > startCollecting)
	        #GC.gc()
                push!(cFArray, cFN)
                push!(dFArray, dFN)
            end

            # update gamma
            if parameters["movingCircle"][1]
                gammaSoA = LightControl.GetMovingCircleSoA(grid, parameters["movingCircle"][2], parameters["r0"], parameters["width"], dt*(t-1), parameters["movingCircle"][3])
            else 
                gammaSoA = MultiplyScalarSoA2D(grid, gammaSoABase, iFun(dt * (t-1)))
            end
    
            # Update concentrations
            ReactAdvDiff.PredictorCorrectorStepCasesRADSoA!(grid, domainType, cF.CDISoA, cFN, gammaSoA, ReactAdvDiff.RDI, cP, cP.DT, dt, bcRAD_X, bcRAD_Y, rDomain)
            ReactAdvDiff.PredictorCorrectorStepCasesRADSoA!(grid, domainType, cF.CDASoA, cFN, gammaSoA, ReactAdvDiff.RDA, cP, cP.DT, dt, bcRAD_X, bcRAD_Y, rDomain)
            ReactAdvDiff.PredictorCorrectorStepCasesRADSoA!(grid, domainType, cF.CBISoA, cFN, gammaSoA, ReactAdvDiff.RBI, cP, 0, dt, bcRAD_X, bcRAD_Y, rDomain)
            ReactAdvDiff.PredictorCorrectorStepCasesRADSoA!(grid, domainType, cF.CBASoA, cFN, gammaSoA, ReactAdvDiff.RBA, cP, 0, dt, bcRAD_X, bcRAD_Y, rDomain)
            ReactAdvDiff.PredictorCorrectorStepCasesRADSoA!(grid, domainType, cF.CCSoA, cFN, gammaSoA, ReactAdvDiff.RC, cP, cP.DC, dt, bcRAD_X, bcRAD_Y, rDomain)
            ReactAdvDiff.PredictorCorrectorStepCasesRADSoA!(grid, domainType, cF.CDSoA, cFN, gammaSoA, ReactAdvDiff.RD, cP, cP.DD, dt, bcRAD_X, bcRAD_Y, rDomain)
            ReactAdvDiff.PredictorCorrectorStepCasesRADSoA!(grid, domainType, cF.CDstSoA, cFN, gammaSoA, ReactAdvDiff.RDst, cP, cP.DD, dt, bcRAD_X, bcRAD_Y, rDomain)

            if !chemOnly
                # Update displacements
                Mechanics.PredictorCorrectorStepDispSoA!(grid, domainType, dF, cFN.CBASoA, cFN.CBISoA, mP, dt, bcU_X, bcU_Y, rDomain)
            end

        end
        println("Done with simulation.")

        ### Return the saved trajectory data, to be selectivey exported from the calling scope
        return (cFArray, dFArray)

    end # InitializeAndRun


end # Main
