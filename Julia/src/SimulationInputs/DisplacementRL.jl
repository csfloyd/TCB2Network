srcPath = joinpath(@__DIR__, "../")
include(srcPath * "RLEnvironment.jl")
include(srcPath * "RLExperiment.jl")
include(srcPath * "Misc.jl")

args = Base.ARGS;

# set input/output
inputFile = string(args[1])
outputDir = string(args[2])

# set default parameters
defaultParams = Dict(
# grid params
"domainType" => "cyl",
"Nx" => 100,
"Ny" => 1,
"dx" => 1.0,
"dt" => 2e-4,
"ndt" => 2500, #2000,
"timeStride" => 1,
"startCollecting" => 0,
"seed" => 0,
    
# boundary conditions
"bcU_X" => "neuCyl",
"bcU_Y" => "neu",

"bcRAD_X" => "neuCyl",
"bcRAD_Y" => "neu",

# chemical parameters
"kFac" => 1.0,
"kAct" => 5e1,
"kInact" => 2e1, 
"kTrap" => 3e2,
"kRel" => 7e2,
"kIBind" => 2e-1, # 2
"kIUnbind" => 0.5e0,
"kABind" => 2e-1, # 2
"kAUnbind" => 0,
"beta" => 0.2,
"CSat" => 5,
    
"DT" => 10, # 8
"DD" => 100, # 120
"DC" => 300,
    
"CDI0" => 1.25, #2.5,
"CDA0" => 0,
"CBI0" => 0.05,
"CBA0" => 0,
"CC0" => 0,
"CD0" => 20,
"CDst0" => 20,

"lambda0" => 15.0,
"mu0" => 15.0,
"gMin" => 0.5,
    
## reward params
"ust" => -0.5,
"rst" => 50,
"ks" => 1e-2,

## parameters of the Environment object
"nSteps" => 50, # number of steps in one episode
"bounds" => [20, 0.15], #[20, 2.0], # range of allowed increments of the activityCoefficients
"rewP" => 1, # strength of the penalty used to compute the reward
    
## parameters of the Experiment object
"eps_or_hrs" => "hrs", # terminate after nEpisodes epsiodes ("eps") or nEpisodes hours ("hrs")
"nEpisodes" => 24, # number of episodes / hours
"stepStride" => 1,
"episodeStride" => 5,
"stateTrajStride" => 5,
"batchSize" => 128, # how many samples to include in replay buffer used to train the networks
"updateFreq" => 10, # how many steps to do before updating the network parameters 
"netLayers" => 1,
"netWidth" => 32,
"gamma" => 0.99f0, #0.99f0, 
"rho" => 0.995f0,
"adamRate" => 5e-4,
    
"act_limit" => 1.0,
"act_noise" => 5e-2,
"annealBool" => false,
"annealTime" => 1000,

## restart parameters
"restartBool" => false, # whether to load experiment from a previous SavedData file or not 
"parentDirectoryName" => "Dirs_seed_ks", # the name of the parent folder containing all the param trials
"restartLabel" => "_R" # the label which has been appended to parentDirectoryName directory manually by the user
);

# set parameters
inputParams = Misc.ParseInput(inputFile);
parameters = deepcopy(defaultParams)
for p in keys(inputParams)
    parameters[p] = inputParams[p]
end

println("Starting experiment.")
if parameters["restartBool"]

    # the next 2 lines convert, e.g., ".../Dirs/Dirs_ks/ks_1/" to ".../Dirs/Dirs_ks_R/ks_1/SavedData.jld2"
    labeledOutputDir = replace(outputDir, parameters["parentDirectoryName"] => parameters["parentDirectoryName"] * parameters["restartLabel"]) 
    pathToLoadRestartFile = labeledOutputDir * "SavedData.jld2" 

    @time ex = RestartExperiment(pathToLoadRestartFile, parameters["nEpisodes"], parameters["eps_or_hrs"]);
else
    @time ex = InitializeAndRunExperiment(parameters);
end

println(outputDir);
SaveResults(ex, parameters, outputDir);

println("Done saving results.")
