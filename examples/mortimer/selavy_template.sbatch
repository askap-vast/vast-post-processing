#!/bin/bash -l
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --ntasks=37
#SBATCH --mem-per-cpu=6500M
#SBATCH --job-name={job_name}
#SBATCH --chdir={working_dir_path}

module load pkgsrc/2022Q1
export PATH=$PATH:/usr/pkg/bin
export SINGULARITYENV_APPEND_PATH=/usr/pkg/bin
export SINGULARITYENV_APPEND_LD_LIBRARY_PATH=/usr/pkg/lib
export containerImage=/raid-17/LS/kaplan/users/obrienan/singularity/yandasoft_1.6.0-mpich

parset={parset_path}
log={log_path}
working_dir={working_dir_path}

srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --job-name={job_name} \
    singularity exec \
        --bind /usr/pkg/bin \
        --bind /raid-17/LS/kaplan/users/obrienan/racs-wip/racs-wip-low-1/neighbours \
        --home ${{working_dir}} \
        $containerImage selavy -c "$parset" >> "$log"
