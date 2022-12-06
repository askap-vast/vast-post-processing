#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --partition=batch

release_epoch="EPOCH00"
neighbours_dirname="neighbours-swarp"

module load pkgsrc/2022Q1
export PATH=$PATH:/usr/pkg/bin
export SINGULARITYENV_APPEND_PATH=/usr/pkg/bin
export SINGULARITYENV_APPEND_LD_LIBRARY_PATH=/usr/pkg/lib
export containerImage=/raid-17/LS/kaplan/users/obrienan/singularity/vast-combined_mortimer
export data_dir=/raid-17/LS/kaplan/users/obrienan/racs-wip/racs-wip-low-1

srun --job-name selavy-config-${release_epoch} \
    singularity exec \
        --bind /usr/pkg/bin \
        --bind ${data_dir} \
        $containerImage python /app/selavy_combined.py ${data_dir}/${neighbours_dirname}/${release_epoch} ${data_dir}/${neighbours_dirname}/selavy_template.in ${data_dir}/${neighbours_dirname}/selavy_template.sbatch --stokes I --racs
