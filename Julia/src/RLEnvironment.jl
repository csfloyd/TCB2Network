using ReinforcementLearning
using Random
using IntervalSets

include("SimMainRL.jl")
using .SimMainRL

include("LightControl.jl")

export TCB2Env

########################################################################################
#
#                    RL environment definitions
#
########################################################################################

    # model is https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl/blob/2e1de3e5b6b8224f50b3d11bba7e1d2d72c6ef7c/src/ReinforcementLearningEnvironments/src/environments/examples/SpeakerListenerEnv.jl

mutable struct TCB2Env <: AbstractEnv 
    sP::SimParams
    sS::SimState
    rP::RewardParams
    done::Bool
    t::Int 
    nStates::Int 
    nActions::Int
    bounds::Vector{Real} # range of increments to activity coefficients
    rewP::Real
    rewTm::Real
    rewTp::Real
    nSteps::Int # max number of steps in RL episode
    reward::Real
    rng::Any
    episodeCount::Int
end

function TCB2Env(sP, rP, rng; # overloaded constructor 
    bounds = ones(4), # used in a tanh clamp from the NN outputs
    rewP = 1.0,
    rewTm = 0, # not used yet
    rewTp = 0, # not used yet
    nSteps = 100
    )      

    nActions = length(bounds)
  
    sS = InitializeSimState(sP, rng, rP) 

    env = TCB2Env( # call default constructor  
        sP, 
        sS,
        rP,
        false,
        0, # t
        1, # nStates
        nActions,
        bounds,
        rewP, 
        rewTm, 
        rewTp,
        nSteps,
        0.0,
        rng,
        0
    )

    return env 
end

function RLBase.reset!(env::TCB2Env)

    env.t = 0
    env.reward = 0.0
    env.sS = InitializeSimState(env.sP, env.rng, env.rP)

end

function RLBase.is_terminated(env::TCB2Env) 
    if (env.t > env.nSteps)
        return true
    else 
        return false 
    end 
end 


function GetState(env; last = false) 

    rst = env.rP.rst
    if last
        u = env.sS.lastdF.uxSoA.Values[rst]
    else
        u = env.sS.dF.uxSoA.Values[rst]
    end

    du = u
    duS = u

    return (du, duS) # return raw and scaled version
end

function RLBase.state(env::TCB2Env)
    return [(GetState(env)[2] - env.rP.ust)] ./ env.rP.ust
end

RLBase.state_space(env::TCB2Env, ::Observation{Any}) = Space(vcat(
    Space([ClosedInterval(-Inf, Inf) for _ in 1:env.nStates])...)) # space of vectors of nStates numbers

RLBase.action_space(env::TCB2Env) = Space(vcat(
    Space([ClosedInterval(-env.bounds[a], env.bounds[a]) for a in 1:env.nActions])...)) # space of vectors of nStates numbers

function gaussian(r, r0, s)
    return exp(-0.5 * ((r - r0) / s)^2) / (s * sqrt(2*pi))
end

gaussianSoAFunc2D(grid, a, r0, s) = ScalarSoA2D([a * gaussian(r, r0, s) for r in 1:grid.Nx, y in 1:1])
tanhSoAFunc2D(grid, a, r0, s) = ScalarSoA2D([a * LightControl.spaceFunc(r, r0, s) for r in 1:grid.Nx, y in 1:1])

function UpdateLightFromAction!(env, action, bounds)


    r0 = env.rP.rst + bounds[1] * action[1]
    a = bounds[2] * 0.5 * (action[2] + 1.0) + 0.01
    s = 5 #bounds[3] * (0.5) * (action[3] + 1.0) + 5.0 # range from 5 to bounds

    #env.sS.gammaSoA = gaussianSoAFunc2D(env.sP.grid, a, r0, s)
    env.sS.gammaSoA = tanhSoAFunc2D(env.sP.grid, a, r0, s)
end

    
function (env::TCB2Env)(action::Vector) # call action vector, wraps function
    UpdateLightFromAction!(env, action, env.bounds)
    _step!(env)
end

function (env::TCB2Env)(action::Real) # call for single action, wraps function
    UpdateLightFromAction!(env, [action], env.bounds)
    _step!(env)
end

function _step!(env::TCB2Env) # wrap SimStep and check number of defects
    SimStep!(env.sS, env.sP)
    env.t += 1
end

function PredictedState(p, env)
    return p + env.sP.ndt * (- env.rP.ks * (p - env.rP.ust)) 
end 

function RLBase.reward(env::TCB2Env) 

    # u = GetState(env)[1]
    # env.reward = - env.rewP * (u - env.rP.ust)^2

    lastu = GetState(env; last = true)[1]
    u = GetState(env)[1]
    pred = PredictedState(lastu, env)
    env.reward = - env.rewP * (u - pred)^2
    
end 

