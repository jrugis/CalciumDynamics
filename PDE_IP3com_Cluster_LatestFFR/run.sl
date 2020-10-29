#!/bin/bash -e
#SBATCH --account=nesi99999
#SBATCH --job-name=psim-mat    #Name to appear in squeue
#SBATCH --time=10:00:00      #Max walltime
#SBATCH --mem=24G          #Max memory
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --partition=large,bigmem

module load MATLAB/2018b

#Job run 
#srun matlab -nodisplay -r "solver('Fc0Fip0t50nc70c02b20437d537da4def222acbffbf90.mat', 0)"
#srun matlab -nodisplay -r "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); solver('Fc0Fip0t50nc70c02b20437d537da4def222acbffbf90.mat', 0)"
srun matlab -nodisplay -r "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); solver('Fc0Fip0t150nc70c02b20437d537da4def222acbffbf90.mat', 0)"
