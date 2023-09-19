ustVec = [-0.4, -0.2, 0.0, 0.2, 0.4]


trialParamDict = Dict(
"ust" => ustVec
)

global baseDir = "/project/svaikunt/csfloyd/RLTCB2/"
global Dirs = baseDir * "Dirs/Dirs_ust/"

function ReplaceBatchRunDelete(batchFile, inputFile, outputDir)
    f = open(batchFile)
    (tmppath, tmpio) = mktemp()
    try
        lines = readlines(f)
        for l in lines
            sl = split(l)
            if (length(sl) > 0) && (sl[1] == "julia")
                ns = replace(l, "inputFile" => inputFile)
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
    newBatch = outputDir * "slatboltz.sh"
    mv(tmppath, newBatch, force = true)
    oldDir = pwd()
    cd(outputDir)
    try
        run(`sbatch slatboltz.sh`)
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


function MakeDirectoriesAndRun(currParamDict, currDir, currInputFile)

    global baseDir

    currParam = collect(keys(currParamDict))[1]

    if length(collect(keys(currParamDict))) == 1 # reached the bottom
        for val in currParamDict[currParam]
            # create the new directory
            newDir = currDir * "$currParam"*"_"*"$val/"
            mkdir(newDir)
            # update the input file
            newInputFile = newDir * "inputFile.txt"
            newLine = string(currParam) * "     " * string(val)
            cp(currInputFile, newInputFile)
            AddLineToFile!(newInputFile, newLine)
            cp(baseDir * "latboltz.sh", newDir * "latboltz.sh", force = true)
            ReplaceBatchRunDelete(newDir * "latboltz.sh", newInputFile, newDir)
            rm(newDir * "latboltz.sh")
        end
        return
    end

    newParamDict = deepcopy(currParamDict)
    delete!(newParamDict, currParam)

    for val in currParamDict[currParam]
        # create the new directory
        newDir = currDir * "$currParam"*"_"*"$val/"
        mkdir(newDir)
        # update the input file
        newInputFile = newDir * "tempInputFile.txt"
        newLine = string(currParam) * "     " * string(val)
        cp(currInputFile, newInputFile)
        AddLineToFile!(newInputFile, newLine)
        # do the recursive call
        MakeDirectoriesAndRun(newParamDict, newDir, newInputFile)
        # delete the temporary input file
        rm(newDir * "tempInputFile.txt")
    end
end


# make sure the directory is empty
try
    rm(Dirs, recursive = true)
catch
end
mkdir(Dirs)

# do it all
MakeDirectoriesAndRun(trialParamDict, Dirs, baseDir * "/baseInput.txt")
