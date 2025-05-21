#!/bin/bash -l
#SBATCH --ntasks=8
##SBATCH --ntasks=50
##SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=4500m
#SBATCH --time=02:00:00
##SBATCH --time=06:00:00
#SBATCH --partition=batch

export release_epoch="EPOCH00"
export neighbours_dirname="neighbours"

module load pkgsrc/2022Q1
export PATH=$PATH:/usr/pkg/bin
export SINGULARITYENV_APPEND_PATH=/usr/pkg/bin
export SINGULARITYENV_APPEND_LD_LIBRARY_PATH=/usr/pkg/lib
export containerImage=/raid-17/LS/kaplan/users/obrienan/singularity/vast-combined_mortimer


srun --mpi=pmi2 --ntasks=$SLURM_NTASKS \
    singularity exec -B /usr/pkg/bin -B /raid-17/LS/kaplan/users/obrienan/racs-wip/racs-wip-low-1 $containerImage \
        python /app/swarp.py /raid-17/LS/kaplan/users/obrienan/racs-wip/racs-wip-low-1/${neighbours_dirname}/${release_epoch} --mpi --racs
