module SimMainRL

    using Random
    using FileIO

    include("MathFunctions.jl")
    using .MathFunctions

    include("SharedStructs.jl")
    using .SharedStructs
    include("ReactAdvDiff.jl")
    include("Mechanics.jl")
    include("LightControl.jl")

    export SimParams,
    SimState,
    RewardParams,
    InitializeSimState,
    PredictedState,
    UpdateLightFromAction!,
    SimStep!

    struct SimParams
        grid::Grid2D
        domainType::String
        dt::Real
        ndt::Int # how many times to integrate in one RL step

        bcU_X::String
        bcU_Y::String 
        bcRAD_X::String
        bcRAD_Y::String 

        cP::Any # chemical params
        mP::Any # mechanical params
        cInits::Any # dict of initial concentrations

        rDomain::Any

        function SimParams(
            domainType, 
            Nx, 
            Ny, 
            dx, 
            dt, 
            ndt,
            bcU_X, 
            bcU_Y, 
            bcRAD_X, 
            bcRAD_Y, 
            cP, 
            mP, 
            cInits)

            return new(Grid2D(Nx, Ny, 1), domainType, dt, ndt, bcU_X, bcU_Y, bcRAD_X, bcRAD_Y, cP, mP, cInits, collect(1:Nx) .* dx)
        end
    end


    mutable struct SimState
        cF::Any # concFields
        dF::Any # dispFields
        gammaSoA::Any # light Field
        lastdF::Any

        function SimState(
            cF,
            dF, 
            gammaSoA)

            return new(cF, dF, gammaSoA, deepcopy(dF))
        end
    end

    struct RewardParams
        ust::Real 
        rst::Real
        ks::Real

        function RewardParams(ust, rst, ks)
            return new(ust, rst, ks)
        end
    end

    function InitializeSimState(sP, rng, fLP) # initialize to zero after 

        ### RAD init
        global CDISoA = ScalarSoA2D(sP.grid, Float64(sP.cInits["CDI0"]))
        global CDASoA = ScalarSoA2D(sP.grid, Float64(sP.cInits["CDA0"]))
        global CBISoA = ScalarSoA2D(sP.grid, Float64(sP.cInits["CBI0"]))
        global CBASoA = ScalarSoA2D(sP.grid, Float64(sP.cInits["CBA0"]))
        global CCSoA = ScalarSoA2D(sP.grid, Float64(sP.cInits["CC0"]))
        global CDSoA = ScalarSoA2D(sP.grid, Float64(sP.cInits["CD0"]))
        global CDstSoA = ScalarSoA2D(sP.grid, Float64(sP.cInits["CDst0"]))
        
        cF = ReactAdvDiff.ConcFields(CDISoA, CDASoA, CBISoA, CBASoA, CCSoA, CDSoA, CDstSoA)

        ### mech init
        global uxSoA = ScalarSoA2D(sP.grid)
        global uySoA = ScalarSoA2D(sP.grid)

        dF = Mechanics.DispFields(uxSoA, uySoA)

        ### get light protocol
        global gammaSoA = ScalarSoA2D(sP.grid, 0.0)

        return SimState(cF, dF, gammaSoA)
    end 

    function SimStep!(sS, sP) # integrates for ndt and updates the fields in sS

        sS.lastdF = deepcopy(sS.dF)
        
        for t in 1:sP.ndt

            # Copy down data to use for the updates of various fields, so that the order of updates doesn't matter
            cFN = deepcopy(sS.cF)
            dFN = deepcopy(sS.dF)
            
            # Update concentrations
            ReactAdvDiff.PredictorCorrectorStepRADCylSoA!(sP.grid, sP.domainType, sS.cF.CDISoA, cFN, sS.gammaSoA, ReactAdvDiff.RDI, sP.cP, sP.cP.DT, sP.dt, sP.bcRAD_X, sP.bcRAD_Y, sP.rDomain)
            ReactAdvDiff.PredictorCorrectorStepRADCylSoA!(sP.grid, sP.domainType, sS.cF.CDASoA, cFN, sS.gammaSoA, ReactAdvDiff.RDA, sP.cP, sP.cP.DT, sP.dt, sP.bcRAD_X, sP.bcRAD_Y, sP.rDomain)
            ReactAdvDiff.PredictorCorrectorStepRADCylSoA!(sP.grid, sP.domainType, sS.cF.CBISoA, cFN, sS.gammaSoA, ReactAdvDiff.RBI, sP.cP, 0, sP.dt, sP.bcRAD_X, sP.bcRAD_Y, sP.rDomain)
            ReactAdvDiff.PredictorCorrectorStepRADCylSoA!(sP.grid, sP.domainType, sS.cF.CBASoA, cFN, sS.gammaSoA, ReactAdvDiff.RBA, sP.cP, 0, sP.dt, sP.bcRAD_X, sP.bcRAD_Y, sP.rDomain)
            ReactAdvDiff.PredictorCorrectorStepRADCylSoA!(sP.grid, sP.domainType, sS.cF.CCSoA, cFN, sS.gammaSoA, ReactAdvDiff.RC, sP.cP, sP.cP.DC, sP.dt, sP.bcRAD_X, sP.bcRAD_Y, sP.rDomain)
            ReactAdvDiff.PredictorCorrectorStepRADCylSoA!(sP.grid, sP.domainType, sS.cF.CDSoA, cFN, sS.gammaSoA, ReactAdvDiff.RD, sP.cP, sP.cP.DD, sP.dt, sP.bcRAD_X, sP.bcRAD_Y, sP.rDomain)
            ReactAdvDiff.PredictorCorrectorStepRADCylSoA!(sP.grid, sP.domainType, sS.cF.CDstSoA, cFN, sS.gammaSoA, ReactAdvDiff.RDst, sP.cP, sP.cP.DD, sP.dt, sP.bcRAD_X, sP.bcRAD_Y, sP.rDomain)

            # Update displacements
            Mechanics.PredictorCorrectorStepDispSoA!(sP.grid, sP.domainType, sS.dF, cFN.CBASoA, cFN.CBISoA, sP.mP, sP.dt, sP.bcU_X, sP.bcU_Y, sP.rDomain)

        end 


    end


end # Main
