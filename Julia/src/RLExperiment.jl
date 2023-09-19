using ReinforcementLearning
using StableRNGs
using Flux
using JLD2
using Flux.Losses
using Flux: params

include("ReactAdvDiff.jl")
include("Mechanics.jl")
include("SharedStructs.jl")
using .SharedStructs


# model is https://juliareinforcementlearning.org/docs/experiments/experiments/Policy%20Gradient/JuliaRL_MADDPG_SpeakerListener/#JuliaRL\\_MADDPG\\_SpeakerListener

# see DDPG code at https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl/blob/main/src/ReinforcementLearningZoo/src/algorithms/policy_gradient/ddpg.jl

# see MADDPG code at https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl/blob/main/src/ReinforcementLearningZoo/src/algorithms/policy_gradient/maddpg.jl


########################################################
#
#                    Define Hooks
#
#########################################################


## RLTCB2 Hook

Base.@kwdef mutable struct EpisodeInformation
    rewards::Vector{Real} = []
end

Base.@kwdef mutable struct RLTCB2Hook <: AbstractHook
    episodeList::Vector{EpisodeInformation} = []
    stateTrajList::Vector{Any} = []
    actorLossList::Vector{Real} = []
    criticLossList::Vector{Real} = []

    storingBool::Bool = true # track whether to store episodeInformation during current episode
    stateTrajBool::Bool = true # track whether to store stateTrajectory during current episode
    stepStride::Int = 1 # store every stepStride steps per episode
    stepCount::Int = 0
    episodeStride::Int = 5 # store every episodeStride steps per experiment
    stateTrajStride::Int = 50 # store every episodeStride steps per experiment
    episodeCount::Int = 0
end


function (hook::RLTCB2Hook)(::PostActStage, agent, env)
    if hook.storingBool && (hook.stepCount % hook.stepStride == 0)
        push!(hook.episodeList[end].rewards, reward(env))
    end 
    if hook.stateTrajBool && (hook.stepCount % hook.stepStride == 0)
        push!(hook.stateTrajList[end], deepcopy(env.sS))
    end 
    hook.stepCount += 1
end

function (hook::RLTCB2Hook)(::PreEpisodeStage, agent, env)

    GC.gc() # run garbage collection

    if hook.episodeCount % hook.episodeStride == 0
        push!(hook.episodeList, EpisodeInformation())
        hook.storingBool = true
    else 
        hook.storingBool = false 
    end
    if hook.episodeCount % hook.stateTrajStride == 0
        push!(hook.stateTrajList, [])
        hook.stateTrajBool = true
    else 
        hook.stateTrajBool = false 
    end
    hook.stepCount = 0
end

function (hook::RLTCB2Hook)(::PostEpisodeStage, agent, env)

    push!(hook.actorLossList, agent.policy.actor_loss)
    push!(hook.criticLossList, agent.policy.critic_loss)

    hook.episodeCount += 1
end



function CreateDDPGExperiment(env; # works for toy or the full RL env
    seed = 0,
    eps_or_hrs = "eps",
    nEpisodes = 10,
    batchSize = 128,
    updateFreq = 100,
    stepStride = 1,
    episodeStride = 5,
    stateTrajStride = 50,
    netLayers = 3,
    netWidth = 64,
    gamma = 0.95f0,
    rho = 0.99f0,
    act_limit = 1.0,
    act_noise = 1e-3,
    annealBool = false,
    annealTime = 100
    )
    
    rng = env.rng
    A = action_space(env)
    ns = length(state(env))
    na = length(A)

    init = glorot_uniform(rng)

    create_actor() = Chain( 
        Dense(ns, netWidth, relu; init = init),
        [Dense(netWidth, netWidth, relu; init = init) for _ in 1:netLayers]...,
        Dense(netWidth, na; init = init)
        ) 


    create_critic() = Chain(
        Dense(ns + na, netWidth, relu; init = init),
        [Dense(netWidth, netWidth, relu; init = init) for _ in 1:netLayers]...,
        Dense(netWidth, 1; init = init),
        ) 

    agent = Agent(
        policy = DDPGPolicy(
            behavior_actor = NeuralNetworkApproximator(
                model = create_actor(),
                #optimizer = ADAM(),
                optimizer = Flux.Optimise.Optimiser(ClipNorm(0.5), ADAM(5e-4)),
            ),
            behavior_critic = NeuralNetworkApproximator(
                model = create_critic(),
                #optimizer = ADAM(),
                optimizer = Flux.Optimise.Optimiser(ClipNorm(0.5), ADAM(5e-4)),
            ),
            target_actor = NeuralNetworkApproximator(
                model = create_actor(),
                #optimizer = ADAM(),
            ),
            target_critic = NeuralNetworkApproximator(
                model = create_critic(),
                #optimizer = ADAM(),
            ),

            γ = gamma, # discount factor to compute TD error
            ρ = rho, # controls the percentage of parameters that don't make it to the target: dest .= ρ .* dest .+ (1 - ρ) .* src
            na = na,
            batch_size = batchSize,

            start_steps = 75*batchSize, # wait until start_steps before learning
            start_policy = RandomPolicy(Space([-act_limit..act_limit for _ in 1:na]); rng = rng), # policy that is used until start_steps
            update_after = 75*batchSize, # wait until update_after before updating
            update_freq = updateFreq,

            act_limit = act_limit, # used to clamp the output of the actor NN, the action vector further passed to a tanh function of unit width with user-specified bounds for each element of the vector
            act_noise = act_noise, # how much noise to add to the actions
            rng = rng,
        ),
        trajectory = CircularArraySARTTrajectory(
            capacity = 1_000_000,
            state = Vector{Float64} => (ns,),
            action = Vector{Float64} => (na, ),
        ),
    )

    if eps_or_hrs == "eps"
        stop_condition = StopAfterEpisode(nEpisodes, is_show_progress=true)
        #stop_condition = StopAfterStep(nEpisodes * nSteps, is_show_progress=true)
    else 
        stop_condition = StopAfterNSeconds(nEpisodes * 3600.0)
    end


    hook = RLTCB2Hook(stepStride = stepStride, episodeStride = episodeStride, stateTrajStride = stateTrajStride)
    tagline = "# TCB2 RL with DDPG"

    Experiment(agent, env, stop_condition, hook, tagline)
end


function SaveResults(ex, parameters, pathToSave)
    save(pathToSave * "SavedData.jld2", 
    "ex", ex, 
    "parameters", parameters)
end

########################################################################################
#
#                    Initialize and run functions
#
########################################################################################

function InitializeAndRunExperiment(parameters)

    #######################################################################################
    #
    #                    Load the parameters into global variables
    #
    ########################################################################################

    println("Initializing...")

    # General params
    global domainType = parameters["domainType"]
    global Nx = parameters["Nx"] 
    if domainType == "2D"
        global Ny = parameters["Ny"]
    else
        global Ny = 1
    end
    global dx = parameters["dx"] # set to 1
    global dt = parameters["dt"]
    global ndt = parameters["ndt"]
    global timeStride = parameters["timeStride"] # save every timeStride frames
    global startCollecting = parameters["startCollecting"]
    global seed = parameters["seed"] # random seed - will use this if this it is not zero

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

    cInits = Dict([
        "CDI0" => parameters["CDI0"],
        "CDA0" => parameters["CDA0"],
        "CBI0" => parameters["CBI0"],
        "CBA0" => parameters["CBA0"],
        "CC0" => parameters["CC0"],
        "CD0" => parameters["CD0"],
        "CDst0" => parameters["CDst0"]
    ])

    ### mech params
    mP = Mechanics.MechParams(parameters["lambda0"], parameters["mu0"], parameters["gMin"])

    sP = SimMainRL.SimParams(domainType, Nx, Ny, dx, dt, ndt, bcU_X, bcU_Y, bcRAD_X, bcRAD_Y, cP, mP, cInits)

    ## parameters of the imposed force law 
    global rst = parameters["rst"] # stiffness times the drag
    global ust = parameters["ust"] # equilibrium separation

    rP = SimMainRL.RewardParams(ust, rst)

    ## parameters of the Environment object
    global nSteps = parameters["nSteps"] # number of steps in one episode
    global bounds = parameters["bounds"] # range of allowed increments of the activityCoefficients
    global rewP = parameters["rewP"] # strength of the penalty used to compute the reward
    rng = StableRNG(seed);

    env = TCB2Env(sP, rP, rng; nSteps = nSteps, bounds = bounds, rewP = rewP)

    ## parameters of the Experiment object
    global eps_or_hrs = parameters["eps_or_hrs"] # terminate after nEpisodes epsiodes ("eps") or nEpisodes hours ("hrs")
    global nEpisodes = parameters["nEpisodes"] # number of episodes / hours 
    global batchSize = parameters["batchSize"] # how many samples to include in replay buffer used to train the networks
    global updateFreq = parameters["updateFreq"] # how many steps to do before updating the network parameters 
    global stepStride = parameters["stepStride"] # clamp on NN action output
    global episodeStride = parameters["episodeStride"] # 
    global stateTrajStride = parameters["stateTrajStride"] #
    global netLayers = parameters["netLayers"] # number of hidden layers in the neural networks 
    global netWidth = parameters["netWidth"] # width of hidden layers in the neural networks
    global gamma = parameters["gamma"] # discount factor
    global rho = parameters["rho"] # soft update factor
    global act_limit = parameters["act_limit"] # clamp on NN action output
    global act_noise = parameters["act_noise"] # action noise scale
    global annealBool = parameters["annealBool"] # whether to anneal the noise in activity
    global annealTime = parameters["annealTime"] # time over which the noise exponentially decay (in units of steps)

    ex = CreateDDPGExperiment(env; seed = seed, eps_or_hrs = eps_or_hrs, nEpisodes = nEpisodes, 
        stepStride = stepStride, episodeStride = episodeStride, stateTrajStride = stateTrajStride, 
        batchSize = batchSize, updateFreq = updateFreq, netLayers = netLayers, netWidth = netWidth, 
        gamma = gamma, rho = rho, act_limit = act_limit, act_noise = act_noise, annealBool = annealBool, annealTime = annealTime)

    ###########################################
    ### Run the experiment and save results ###
    ###########################################

    run(ex)
    
    return ex
end


function RestartExperiment(pathToLoadRestartFile, nEpisodes, eps_or_hrs = "eps")

    ########################################################
    ### Load previous experiment, reset counter, and run ###
    ########################################################
    
    d = load(pathToLoadRestartFile)
    ex = d["ex"]

    if eps_or_hrs == "eps"
        ex.stop_condition = StopAfterEpisode(nEpisodes, is_show_progress=true)
    else 
        ex.stop_condition = StopAfterNSeconds(nEpisodes * 3600.0)
    end

    run(ex)
    
    return ex
end
