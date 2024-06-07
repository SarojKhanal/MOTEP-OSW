#!/bin/bash -l
#SBATCH --job-name=EIProject     # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=48       # cpu-cores per task (>1 if multi-threaded tasks) (= 192 [GB] / (ntasks * mem-per-cpu [GB]))
#SBATCH --time=71:59:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=skhanal8@jhu.edu

module purge
module load julia/1.7.3
module load gurobi/9.5.2

julia TEP.jl
