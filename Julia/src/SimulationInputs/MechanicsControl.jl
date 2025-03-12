srcPath = joinpath(@__DIR__, "../")
include(srcPath * "SimMainMechanicsControl.jl")
include(srcPath * "Misc.jl")
using FileIO

args = Base.ARGS;

# set input/output
inputFile = string(args[1])
outputDir = string(args[2])

# set default parameters
defaultParams = Dict(
# grid params
"domainType" => "cyl",
"Nx" => 250, # 400
"Ny" => 1,
"dx" => 1.0,
"dt" => 1.0e-4,
"nSteps" => 100e4,
"timeStride" => 2e3,
"startCollecting" => 0,
"seed" => 0,

"CBScale" => 2.0,
"aOff" => 0.5,
"flat" => 0.2,
"n" => 2,
"r0Base" => 37.5,
"widthBase" => 1.0,
"snapr0" => false,

"r0" => 37.5,
"width" => 4.0,

"len" => 1,
"cyc" => 30,
"del" => 0.5,

# boundary conditions
"bcU_X" => "neuCyl",
"bcU_Y" => "neu",

"bcRAD_X" => "neuCyl",
"bcRAD_Y" => "neu",

"lambda0" => 3.0,
"mu0" => 3.0,
"gMin" => 0.5

);

# set parameters
inputParams = Misc.ParseInput(inputFile);
parameters = deepcopy(defaultParams)
for p in keys(inputParams)
    parameters[p] = inputParams[p]
end

(dFArray) = SimMainMechanicsControl.InitializeAndRun(parameters);

# save results
pathName = outputDir *  "SavedData" * ".jld2"

@time save(pathName,
"parameters", parameters,
"dFArray", dFArray
)

println("Done saving results.")
