srcPath = joinpath(@__DIR__, "../")
include(srcPath * "SimMain.jl")
include(srcPath * "Misc.jl")
using FileIO

args = Base.ARGS;

# set input/output
inputFile = string(args[1])
outputDir = string(args[2])

# set default parameters
defaultParams = Dict(
# grid params
"domainType" => "2D",
"Nx" => 600, # 400
"Ny" => 600,
"dx" => 1.0,
"dt" => 1.0e-4,
"nSteps" => 5e4,
"timeStride" => 1e4,
"startCollecting" => 0,
"seed" => 0,
    
# boundary conditions
"bcU_X" => "neu",
"bcU_Y" => "neu",

"bcRAD_X" => "neu",
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

"chemOnly" => false,
    
# light control
"width" => 4,
"r0" => 50,
"starBool" => true,
    
#"width" => 4,
#"r0" => 50, #75,
#"starBool" => false,

"movingCircle" => (false, 50, 2.5),


"len" => 100,
"cyc" => 30e6,
"del" => 0.5,

#"len" => 1,
#"cyc" => 30,
#"del" => 0.5,

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

if false
#   parameters["kABind"] = parameters["kFac"] * parameters["kABind"]
#   parameters["kIBind"] = parameters["kFac"] * parameters["kIBind"]
#   parameters["kAUnbind"] = parameters["kFac"] * parameters["kAUnbind"]
#   parameters["kIUnbind"] = parameters["kFac"] * parameters["kIUnbind"]
#   parameters["kAct"] = parameters["kFac"] * parameters["kAct"]
#   parameters["kInact"] = parameters["kFac"] * parameters["kInact"]
#   parameters["kTrap"] = parameters["kFac"] * parameters["kTrap"]
#   parameters["kRel"] = parameters["kFac"] * parameters["kRel"]
#    parameters["DT"] = parameters["kFac"] * parameters["DT"]
#    parameters["DD"] = parameters["kFac"] * parameters["DD"]
#    parameters["DC"] = parameters["kFac"] * parameters["DC"]
 end

(cFArray, dFArray) = SimMain.InitializeAndRun(parameters);

# save results
pathName = outputDir *  "SavedData" * ".jld2"

@time save(pathName,
"parameters", parameters,
"cFArray", cFArray,
"dFArray", dFArray
)

println("Done saving results.")
