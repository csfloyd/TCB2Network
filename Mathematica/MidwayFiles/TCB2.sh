#!/bin/bash
#SBATCH --job-name=TCB2
#SBATCH --output=TCB2_batch.out
#SBATCH --error=TCB2_batch.err
#SBATCH --time=32:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-svaikunt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64000

module unload mathematica
module load mathematica
math -script /project/svaikunt/csfloyd/TCB2/src/TCB2Solve.m tempVal outputDir

