using ReinforcementLearning
using Random
using IntervalSets

include("SimMainRL.jl")
using .SimMainRL

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

function pmVec(agentHandler) # points from plus to minus
    return agentHandler.MinusDefects[1].Position .- agentHandler.PlusDefects[1].Position
end 

function GetState(env) 

    rst = env.rP.rst
    u = env.sS.dF.uxSoA.Values[rst]
    du = u - env.rP.ust
    duS = du / env.rP.ust 

    return (du, duS) # return raw and scaled version
end

function RLBase.state(env::TCB2Env)
    return GetState(env)[2]
end

RLBase.state_space(env::TCB2Env, ::Observation{Any}) = Space(vcat(
    Space([ClosedInterval(-Inf, Inf) for _ in 1:env.nStates])...)) # space of vectors of nStates numbers

RLBase.action_space(env::TCB2Env) = Space(vcat(
    Space([ClosedInterval(-env.bounds[a], env.bounds[a]) for a in 1:env.nActions])...)) # space of vectors of nStates numbers

function gaussian(r, r0, s)
    return exp(-0.5 * ((r - r0) / s)^2) / (s * sqrt(2*pi))
end

gaussianSoAFunc2D(grid, a, r0, s) = ScalarSoA2D([a * gaussian(r, r0, s) for r in 1:grid.Nx, y in 1:1])

function UpdateLightFromAction!(env, action, bounds)

    a = bounds[1] * 0.5 * (action[1] + 1.0)
    r0 = env.rP.rst + bounds[2] * 0.5 * (action[2] + 1.0)
    s = bounds[2] * 0.5 * (action[3] + 1.0)

    env.sS.gammaSoA = gaussianSoAFunc2D(env.sP.grid, a, r0, s)
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


function RLBase.reward(env::TCB2Env) 

    du = GetState(env)[1]

    env.reward = - env.rewP * du^2
end 

