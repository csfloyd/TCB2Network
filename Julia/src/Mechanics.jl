module Mechanics

    include("SharedStructs.jl")
    using .SharedStructs

    include("MathFunctions.jl")
    using .MathFunctions

    mutable struct DispFields
        uxSoA::Any
        uySoA::Any
        function DispFields(uxSoA, uySoA)
            return new(uxSoA, uySoA)
        end 
    end

    struct MechParams
        lambda0::Real
        mu0::Real 
        gMin::Real
        function MechParams(lambda0, mu0, gMin)
            return new(lambda0, mu0, gMin)
        end
    end

    pAFunc(CBASoA, CBISoA) =  CBASoA.Values ./ (CBASoA.Values .+ CBISoA.Values)

    gFunc(pA, gMin) = -(1 - gMin) .* pA

    function gradgFunc2D(grid, CBASoA, CBISoA, gMin, bcx, bcy)
        bcDerivX = BCDerivDict[bcx]
        bcDerivY = BCDerivDict[bcy]

        pA = pAFunc(CBASoA, CBISoA)
        pref = (gMin - 1) ./ (CBASoA.Values .+ CBISoA.Values)
        dxCBA = bcDerivX(CBASoA.Values, FiniteDifferenceX) ./ (2.0 * grid.dx)
        dyCBA = bcDerivY(CBASoA.Values, FiniteDifferenceY) ./ (2.0 * grid.dx)
        dxCBI = bcDerivX(CBISoA.Values, FiniteDifferenceX) ./ (2.0 * grid.dx)
        dyCBI = bcDerivY(CBISoA.Values, FiniteDifferenceY) ./ (2.0 * grid.dx)

        dxg = pref .* ((1 .- pA) .* dxCBA .- pA .* dxCBI)
        dyg = pref .* ((1 .- pA) .* dyCBA .- pA .* dyCBI)

        return (dxg, dyg)
    end

    function gradgFuncCyl(grid, CBASoA, CBISoA, gMin, bcx, bcy)
        bcDerivX = BCDerivDict[bcx]

        pA = pAFunc(CBASoA, CBISoA)
        pref = (gMin - 1) ./ (CBASoA.Values .+ CBISoA.Values)
        dCBA = bcDerivX(CBASoA.Values, FiniteDifferenceCyl) ./ (2.0 * grid.dx)
        dCBI = bcDerivX(CBISoA.Values, FiniteDifferenceCyl) ./ (2.0 * grid.dx)

        dg = pref .* ((1 .- pA) .* dCBA .- pA .* dCBI)

        return dg
    end

    function StrainDerivFuncs(grid, dF, bcx, bcy)
        bcDerivX = BCDerivDict[bcx]
        bcDerivY = BCDerivDict[bcy]

        dxux = bcDerivX(dF.uxSoA.Values, FiniteDifferenceX) ./ (2.0 * grid.dx)
        dyux = bcDerivY(dF.uxSoA.Values, FiniteDifferenceY) ./ (2.0 * grid.dx)
        dxuy = bcDerivX(dF.uySoA.Values, FiniteDifferenceX) ./ (2.0 * grid.dx)
        dyuy = bcDerivY(dF.uySoA.Values, FiniteDifferenceY) ./ (2.0 * grid.dx)

        sxx = dxux
        sxy = 0.5 .* (dxuy .+ dyux)
        syy = dyuy

        divu = dxux .+ dyuy

        gradxTr = bcDerivX(divu, FiniteDifferenceX) ./ (2.0 * grid.dx)
        gradyTr = bcDerivY(divu, FiniteDifferenceY) ./ (2.0 * grid.dx)

        lapux = (bcDerivX(dF.uxSoA.Values, FiniteSecondDifferenceX) .+ bcDerivY(dF.uxSoA.Values, FiniteSecondDifferenceY)) ./ (grid.dx)^2
        lapuy = (bcDerivX(dF.uySoA.Values, FiniteSecondDifferenceX) .+ bcDerivY(dF.uySoA.Values, FiniteSecondDifferenceY)) ./ (grid.dx)^2
        
        return (sxx, sxy, syy, divu, gradxTr, gradyTr, lapux, lapuy)
    end

    function GetForces2D(grid, dispFields, CBASoA, CBISoA, mechParams, bcx, bcy)
        bcDerivX = BCDerivDict[bcx]
        bcDerivY = BCDerivDict[bcy]

        (sxx, sxy, syy, divu, gradxTr, gradyTr, lapux, lapuy) = StrainDerivFuncs(grid, dispFields, bcx, bcy)
        CB = CBASoA.Values .+ CBISoA.Values
        pA = pAFunc(CBASoA, CBISoA)
        g = gFunc(pA, mechParams.gMin)
        (dxg, dyg) = gradgFunc2D(grid, CBASoA, CBISoA, mechParams.gMin, bcx, bcy)

        dxCB = bcDerivX(CB, FiniteDifferenceX) ./ (2.0 * grid.dx)
        dyCB = bcDerivY(CB, FiniteDifferenceY) ./ (2.0 * grid.dx)

        forcex = 2 .* mechParams.mu0 .* (((dxCB .* sxx .+ dyCB .* sxy) .- 0.5 .* g .* dxCB) .+ CB .* (0.5 .* (gradxTr .+ lapux) .- 0.5 .* dxg)) .+
            mechParams.lambda0 .* (dxCB .* (divu .- g) .+ CB .* (gradxTr .- dxg));
        forcey = 2 .* mechParams.mu0 .* (((dxCB .* sxy .+ dyCB .* syy) .- 0.5 .* g .* dyCB) .+ CB .* (0.5 .* (gradyTr .+ lapuy) .- 0.5 .* dyg)) .+
            mechParams.lambda0 .* (dyCB .* (divu .- g) .+ CB .* (gradyTr .- dyg));

        return (forcex, forcey)

    end

    function GetForcesCyl(grid, dispFields, CBASoA, CBISoA, mechParams, bcx, bcy, rDomain)
        bcDerivX = BCDerivDict[bcx]
        bcDerivY = BCDerivDict[bcy]
        dF = dispFields

        CB = CBASoA.Values .+ CBISoA.Values
        pA = CBASoA.Values ./ CB
        g = gFunc(pA, mechParams.gMin)
        dg = gradgFuncCyl(grid, CBASoA, CBISoA, mechParams.gMin, bcx, bcy)

        dCB = bcDerivX(CB, FiniteDifferenceCyl) ./ (2.0 * grid.dx)
        du = bcDerivX(dF.uxSoA.Values, FiniteDifferenceCyl) ./ (2.0 * grid.dx)
        ubr = dF.uxSoA.Values ./ rDomain
        dubr = du ./ rDomain
        dusq = bcDerivX(dF.uxSoA.Values, FiniteSecondDifferenceCyl) ./ (grid.dx^2)
        ubrsq = dF.uxSoA.Values ./ (rDomain .^ 2)

        forcer = 2 .* mechParams.mu0 .* (dCB .* (du .- 0.5 .* g) .+ CB .* (dusq .+ dubr .- ubrsq .- 0.5 .* dg)) .+
            mechParams.lambda0 .* (dCB .* (ubr .+ du .- g) .+ CB .* (dusq .+ dubr .- ubrsq .- dg))

        return forcer

    end

    function PredictorCorrectorStepDispSoA!(grid, domainType, dispFields, CBASoA, CBISoA, mechParams, dt, bcx, bcy, rDomain = [])
        if domainType == "2D"
            PredictorCorrectorStepDisp2DSoA!(grid, dispFields, CBASoA, CBISoA, mechParams, dt, bcx, bcy)
        else 
            PredictorCorrectorStepDispCylSoA!(grid, dispFields, CBASoA, CBISoA, mechParams, dt, bcx, bcy, rDomain)
        end
    end

    function PredictorCorrectorStepDisp2DSoA!(grid, dispFields, CBASoA, CBISoA, mechParams, dt, bcx, bcy)

        b = GetBoundariesFromConditions(grid, bcx, bcy)

        predDispFields = deepcopy(dispFields)
        
        (rhs1x, rhs1y) = GetForces2D(grid, dispFields, CBASoA, CBISoA, mechParams, bcx, bcy)
        predDispFields.uxSoA.Values[b[1]:b[2],b[3]:b[4]] .= predDispFields.uxSoA.Values[b[1]:b[2],b[3]:b[4]] .+ dt .* rhs1x[b[1]:b[2],b[3]:b[4]]
        predDispFields.uySoA.Values[b[1]:b[2],b[3]:b[4]] .= predDispFields.uySoA.Values[b[1]:b[2],b[3]:b[4]] .+ dt .* rhs1y[b[1]:b[2],b[3]:b[4]]

        (rhs2x, rhs2y) = GetForces2D(grid, predDispFields, CBASoA, CBISoA, mechParams, bcx, bcy)
        dispFields.uxSoA.Values[b[1]:b[2],b[3]:b[4]] .= dispFields.uxSoA.Values[b[1]:b[2],b[3]:b[4]] .+ (0.5 .* dt .* (rhs1x[b[1]:b[2],b[3]:b[4]] .+ rhs2x[b[1]:b[2],b[3]:b[4]]))
        dispFields.uySoA.Values[b[1]:b[2],b[3]:b[4]] .= dispFields.uySoA.Values[b[1]:b[2],b[3]:b[4]] .+ (0.5 .* dt .* (rhs1y[b[1]:b[2],b[3]:b[4]] .+ rhs2y[b[1]:b[2],b[3]:b[4]]))
    end

    function PredictorCorrectorStepDispCylSoA!(grid, dispFields, CBASoA, CBISoA, mechParams, dt, bcx, bcy, rDomain)

        b = GetBoundariesFromConditions(grid, bcx, bcy)

        predDispFields = deepcopy(dispFields)
        
        rhs1 = GetForcesCyl(grid, dispFields, CBASoA, CBISoA, mechParams, bcx, bcy, rDomain)
        predDispFields.uxSoA.Values[b[1]:b[2],1] .= predDispFields.uxSoA.Values[b[1]:b[2],1] .+ dt .* rhs1[b[1]:b[2],1]
       
        rhs2 = GetForcesCyl(grid, predDispFields, CBASoA, CBISoA, mechParams, bcx, bcy, rDomain)
        dispFields.uxSoA.Values[b[1]:b[2],1] .= dispFields.uxSoA.Values[b[1]:b[2],1] .+ (0.5 .* dt .* (rhs1[b[1]:b[2],1] .+ rhs2[b[1]:b[2],1]))
        
    end
    

end
