module ReactAdvDiff

    include("SharedStructs.jl")
    using .SharedStructs

    include("MathFunctions.jl")
    using .MathFunctions

    mutable struct ConcFields
        CDISoA::Any
        CDASoA::Any
        CBISoA::Any
        CBASoA::Any
        CCSoA::Any
        CDSoA::Any
        CDstSoA::Any
        function ConcFields(CDISoA, CDASoa, CBISoA, CBASoA, CCSoA, CDSoA, CDstSoA)
            return new(CDISoA, CDASoa, CBISoA, CBASoA, CCSoA, CDSoA, CDstSoA)
        end 
    end

    function UpdateConcFields(cF, cFN)
        for key in keys(cFN)
            SetScalarFromSoA2D!(cF[key], cFN[key])
        end 
    end

    struct ChemParams
        kAct::Real
        kInact::Real 
        kTrap::Real
        kRel::Real
        kIBind::Real 
        kIUnbind::Real 
        kABind::Real 
        kAUnbind::Real 
        beta::Real  
        CSat::Real 
        function ChemParams(kAct, kInact, kTrap, kRel, kIBind, kIUnbind, kABind, kAUnbind, beta, CSat)
            return new(kAct, kInact, kTrap, kRel, kIBind, kIUnbind, kABind, kAUnbind, beta, CSat)
        end
    end

    function BS(CBASoA, CBISoA, CSat)
        CBvals = CBASoA.Values .+ CBISoA.Values
        pAvals = CBASoA.Values ./ CBvals
        BSvals = (CSat^2) / 4 .- (CBvals .- CSat/2).^2
        return (ScalarSoA2D(pAvals), ScalarSoA2D(BSvals))
    end

    RDI(phi, cF, gSoA, cP) = -cP.kAct .* phi.Values .* cF.CCSoA.Values .+ 
        cP.kInact .* cF.CDASoA.Values .- cP.kIBind .* phi.Values .* BS(cF.CBASoA, cF.CBISoA, cP.CSat)[2].Values .+ 
        cP.kIUnbind .* cF.CBISoA.Values

    RDA(phi, cF, gSoA, cP) = cP.kAct .* cF.CDISoA.Values .* cF.CCSoA.Values .- 
        cP.kInact .* phi.Values .- cP.kABind .* phi.Values .* BS(cF.CBASoA, cF.CBISoA, cP.CSat)[2].Values .+ 
        cP.kAUnbind .* cF.CBASoA.Values

    RBI(phi, cF, gSoA, cP) = -cP.kAct .* phi.Values .* cF.CCSoA.Values .+ 
        cP.kInact .* cF.CBASoA.Values .+ cP.kIBind .* cF.CDISoA.Values .* BS(cF.CBASoA, phi, cP.CSat)[2].Values .- 
        cP.kIUnbind .* phi.Values
    
    RBA(phi, cF, gSoA, cP) = cP.kAct .* cF.CBISoA.Values .* cF.CCSoA.Values .- 
        cP.kInact .* phi.Values .+ cP.kABind .* cF.CDASoA.Values .* BS(phi, cF.CBISoA, cP.CSat)[2].Values .- 
        cP.kAUnbind .* phi.Values

    RC(phi, cF, gSoA, cP) = - cP.kAct .* (cF.CDISoA.Values .+ cF.CBISoA.Values) .* phi.Values .+
        cP.kInact .* (cF.CDASoA.Values .+ cF.CBASoA.Values) .- cP.kTrap .* cF.CDstSoA.Values .* phi.Values .+
        cP.kRel .* gSoA.Values .* cF.CDSoA.Values 
    
    RD(phi, cF, gSoA, cP) = cP.kTrap .* cF.CDstSoA.Values .* cF.CCSoA.Values .-
        cP.kRel .* gSoA.Values .* phi.Values 

    RDst(phi, cF, gSoA, cP) = - cP.kTrap .* phi.Values .* cF.CCSoA.Values .+
        cP.kRel .* gSoA.Values .* cF.CDSoA.Values .* cP.beta

    function DiffTerm2D(grid, phiSoA, D, bcDerivX, bcDerivY)
        return ScalarSoA2D(D .* ((bcDerivX(phiSoA.Values, FiniteSecondDifferenceX) .+ bcDerivY(phiSoA.Values, FiniteSecondDifferenceY)) ./ (grid.dx^2)))
    end

    function DiffTermCyl(grid, rDomain, phiSoA, D, bcDerivX, bcDerivY)
        return ScalarSoA2D(D .* (bcDerivX(phiSoA.Values, FiniteSecondDifferenceCyl) ./ (grid.dx^2) .+ bcDerivX(phiSoA.Values, FiniteDifferenceCyl) ./ (2.0 .* grid.dx .* rDomain)))
    end

    function PredictorCorrectorStepRADSoA!(grid, domainType, phiSoA, concFields, gammaSoA, reactionFunc, chemParams, D, dt, bcx, bcy, rDomain = [])
        if domainType == "2D"
            PredictorCorrectorStepRAD2DSoA!(grid, domainType, phiSoA, concFields, gammaSoA, reactionFunc, chemParams, D, dt, bcx, bcy)
        else 
            PredictorCorrectorStepRADCylSoA!(grid, domainType, phiSoA, concFields, gammaSoA, reactionFunc, chemParams, D, dt, bcx, bcy, rDomain)
        end
        
    end

    function PredictorCorrectorStepRAD2DSoA!(grid, domainType, phiSoA, concFields, gammaSoA, reactionFunc, chemParams, D, dt, bcx, bcy)

        bcDerivX = BCDerivDict[bcx]
        bcDerivY = BCDerivDict[bcy]
        b = GetBoundariesFromConditions(grid, bcx, bcy)

        predPhiSoA = deepcopy(phiSoA)

        rhs1 = ScalarSoA2D(reactionFunc(predPhiSoA, concFields, gammaSoA, chemParams)) 
        if D != 0
            diffSoA = DiffTerm2D(grid, predPhiSoA, D, bcDerivX, bcDerivY)
            AddScalarSoA2D!(rhs1, diffSoA)
        end
        predPhiSoA.Values[b[1]:b[2],b[3]:b[4]] .= predPhiSoA.Values[b[1]:b[2],b[3]:b[4]] .+ dt .* rhs1.Values[b[1]:b[2],b[3]:b[4]]

        rhs2 = ScalarSoA2D(reactionFunc(predPhiSoA, concFields, gammaSoA, chemParams)) 
        if D != 0
            diffSoA = DiffTerm2D(grid, predPhiSoA, D, bcDerivX, bcDerivY)
            AddScalarSoA2D!(rhs2, diffSoA)
        end

        phiSoA.Values[b[1]:b[2],b[3]:b[4]] .= phiSoA.Values[b[1]:b[2],b[3]:b[4]] .+ (0.5 .* dt .* (rhs1.Values[b[1]:b[2],b[3]:b[4]] .+ rhs2.Values[b[1]:b[2],b[3]:b[4]]))
    end

    function PredictorCorrectorStepRADCylSoA!(grid, domainType, phiSoA, concFields, gammaSoA, reactionFunc, chemParams, D, dt, bcx, bcy, rDomain)

        bcDerivX = BCDerivDict[bcx]
        bcDerivY = BCDerivDict[bcy]
        b = GetBoundariesFromConditions(grid, bcx, bcy)

        predPhiSoA = deepcopy(phiSoA)

        rhs1 = ScalarSoA2D(reactionFunc(predPhiSoA, concFields, gammaSoA, chemParams)) 
        if D != 0
            diffSoA = DiffTermCyl(grid, rDomain, predPhiSoA, D, bcDerivX, bcDerivY)
            AddScalarSoA2D!(rhs1, diffSoA)
        end
        predPhiSoA.Values[b[1]:b[2],1] .= predPhiSoA.Values[b[1]:b[2],1] .+ dt .* rhs1.Values[b[1]:b[2],1]

        rhs2 = ScalarSoA2D(reactionFunc(predPhiSoA, concFields, gammaSoA, chemParams)) 
        if D != 0
            diffSoA = DiffTermCyl(grid, rDomain, predPhiSoA, D, bcDerivX, bcDerivY)
            AddScalarSoA2D!(rhs2, diffSoA)
        end

        phiSoA.Values[b[1]:b[2],1] .= phiSoA.Values[b[1]:b[2],1] .+ (0.5 .* dt .* (rhs1.Values[b[1]:b[2],1] .+ rhs2.Values[b[1]:b[2],b[3]:b[4]]))
    end

end
