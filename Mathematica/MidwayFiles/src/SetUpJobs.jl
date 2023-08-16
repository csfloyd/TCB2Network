function CreateLists2(v1Min, v1Step, v1Num, v2Min, v2Step, v2Num)
    retList = [];
    for i in 1:v1Num, j in 1:v2Num
        v1 = v1Min + v1Step * (i-1)
        v2 = v2Min + v2Step * (j-1)
        push!(retList, "List[" * string(v1) * "," * string(v2) * "]")
    end
    return retList
end

DCVec = [20, 40, 60, 80, 100];
r0Vec = [25, 50, 75, 100, 125, 150];
facVec = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0];
facVec = [0, 1, 2, 3, 4, 5]

facVec = CreateLists2(0, 1, 6, 0, 1, 6)

trialParamDict = Dict(
#"DC" => DCVec
#"r0" => r0Vec
"fac" => facVec
)

global baseDir = "/project/svaikunt/csfloyd/TCB2/"
global Dirs = baseDir * "Dirs/Dirs_lambda_mu_cyl/"

function ReplaceBatchRunDelete(batchFile, val, outputDir)
    f = open(batchFile)
    (tmppath, tmpio) = mktemp()
    try
        lines = readlines(f)
        for l in lines
            sl = split(l)
            if (length(sl) > 0) && (sl[1] == "math")
                ns = replace(l, "tempVal" => val)
                ns = replace(ns, "outputDir" => outputDir)
                write(tmpio, ns)
                write(tmpio, "\n")

            else
                write(tmpio, l)
                write(tmpio, "\n")

            end
        end
    finally
        close(f)
        close(tmpio)
    end
    newBatch = outputDir * "sTCB2.sh"
    mv(tmppath, newBatch, force = true)
    oldDir = pwd()
    cd(outputDir)
    try
        run(`sbatch sTCB2.sh`)
    catch
        println("Failed to submit the batch job.")
    end
    cd(oldDir)
    rm(newBatch, force = true)
end

function AddLineToFile!(inputFile, newLine)
    f = open(inputFile, "a")
    write(f, newLine)
    write(f, "\n")
    close(f)
end


function MakeDirectoriesAndRun(currParamDict, currDir)

    global baseDir

    currParam = collect(keys(currParamDict))[1]

    if length(collect(keys(currParamDict))) == 1 # reached the bottom
        for val in currParamDict[currParam]
            # create the new directory
            newDir = currDir * "$currParam"*"_"*"$val/"
            mkdir(newDir)
            cp(baseDir * "TCB2.sh", newDir * "TCB2.sh", force = true)
            ReplaceBatchRunDelete(newDir * "TCB2.sh", string(val), newDir)
            rm(newDir * "TCB2.sh")
        end
        return
    end

    newParamDict = deepcopy(currParamDict)
    delete!(newParamDict, currParam)

    for val in currParamDict[currParam]
        # create the new directory
        newDir = currDir * "$currParam"*"_"*"$val/"
        mkdir(newDir)
        # do the recursive call
        MakeDirectoriesAndRun(newParamDict, newDir)
    end
end


# make sure the directory is empty
try
    rm(Dirs, recursive = true)
catch
end
mkdir(Dirs)

# do it all
MakeDirectoriesAndRun(trialParamDict, Dirs)
