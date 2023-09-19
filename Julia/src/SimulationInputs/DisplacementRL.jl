srcPath = joinpath(@__DIR__, "../")
include(srcPath * "RLEnvironment.jl")
include(srcPath * "RLExperiment.jl")
include(srcPath * "Misc.jl")

args = Base.ARGS;

# set input/output
inputFile = string(args[1])
outputDir = string(args[2])

# set default parameters
parameters = Dict(
# grid params
"domainType" => "cyl",
"Nx" => 150,
"Ny" => 1,
"dx" => 1.0,
"dt" => 2e-4,
"ndt" => 1000,
"timeStride" => 1,
"startCollecting" => 0,
"seed" => 0,
    
# boundary conditions
"bcU_X" => "neuCyl",
"bcU_Y" => "neu",

"bcRAD_X" => "neuCyl",
"bcRAD_Y" => "neu",

# chemical parameters
"kAct" => 5e1,
"kInact" => 5e0,
"kTrap" => 3e2,
"kRel" => 7e2,
"kIBind" => 1e-1,
"kIUnbind" => 1e0,
"kABind" => 1e-1,
"kAUnbind" => 0,
"beta" => 1.0,
"CSat" => 5,
    
"DT" => 1,
"DD" => 25,
"DC" => 500,
    
"CDI0" => 2.5,
"CDA0" => 0,
"CBI0" => 0.1,
"CBA0" => 0,
"CC0" => 0,
"CD0" => 20,
"CDst0" => 20,

"lambda0" => 3.0,
"mu0" => 3.0,
"gMin" => 0.5,
    
## reward params
"ust" => -0.25,
"rst" => 75,

## parameters of the Environment object
"nSteps" => 100, # number of steps in one episode
"bounds" => [1.0, 15.0, 5.0], # range of allowed increments of the activityCoefficients
"rewP" => 1, # strength of the penalty used to compute the reward
    
## parameters of the Experiment object
"eps_or_hrs" => "eps", # terminate after nEpisodes epsiodes ("eps") or nEpisodes hours ("hrs")
"nEpisodes" => 1000, # number of episodes / hours
"stepStride" => 1,
"episodeStride" => 5,
"stateTrajStride" => 5,
"batchSize" => 128, # how many samples to include in replay buffer used to train the networks
"updateFreq" => 10, # how many steps to do before updating the network parameters 
"netLayers" => 1,
"netWidth" => 32,
"gamma" => 0.99f0, 
"rho" => 0.9995f0,
    
"act_limit" => 1.0,
"act_noise" => 0e-1,
"annealBool" => false,
"annealTime" => 250,

## restart parameters
"restartBool" => false, # whether to load experiment from a previous SavedData file or not 
"parentDirectoryName" => "Dirs_ks_cutoff", # the name of the parent folder containing all the param trials
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
