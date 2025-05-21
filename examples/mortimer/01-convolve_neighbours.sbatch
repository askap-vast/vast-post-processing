#!/bin/bash -l
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=18000m
#SBATCH --time=07:00:00
#SBATCH --partition=768g

# one convolution requires approx 16.21 GB memory
# for all 903 RACS fields: --ntasks=24 --mem=390g --time=07:00:00 --partition=768g 

release_epoch="EPOCH00"

module load pkgsrc/2022Q1
export PATH=$PATH:/usr/pkg/bin
export SINGULARITYENV_APPEND_PATH=/usr/pkg/bin
export SINGULARITYENV_APPEND_LD_LIBRARY_PATH=/usr/pkg/lib
export containerImage=/raid-17/LS/kaplan/users/obrienan/singularity/vast-combined_mortimer

srun --mpi=pmi2 --ntasks=$SLURM_NTASKS \
    singularity exec -B /usr/pkg/bin:/usr/pkg/bin -B /raid-17/LS/kaplan/users/obrienan/racs-wip/racs-wip-low-1/neighbours:/data $containerImage \
        python /app/convolve_neighbours.py /data/${release_epoch} --mpi --racs
        
        # or, if only a subset of fields is required ...
        # python /app/convolve_neighbours.py /data/${release_epoch} --mpi --racs \
        #     --field RACS_1824+00A \
        #     --field RACS_0703+31A \
        #     --field RACS_1435-31A \
        #     --field RACS_2157-37A \
        #     --field RACS_1430-68A
