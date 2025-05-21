#!/bin/bash -l
#SBATCH --partition=batch
#SBATCH --time=01:30:00
#SBATCH --ntasks=37
#SBATCH --mem-per-cpu=6500M
#SBATCH --job-name=selavy-RACS-EPOCH00
#SBATCH --array=0-4%5

module load pkgsrc/2022Q1
export PATH=$PATH:/usr/pkg/bin
export SINGULARITYENV_APPEND_PATH=/usr/pkg/bin
export SINGULARITYENV_APPEND_LD_LIBRARY_PATH=/usr/pkg/lib
export containerImage=/raid-17/LS/kaplan/users/obrienan/singularity/yandasoft_1.6.0-mpich

field_list_file=${1:-"field_list.txt"}
release_epoch="EPOCH00"

mapfile -t field_list < ${field_list_file}
field=${field_list[${SLURM_ARRAY_TASK_ID}]}
working_dir=/raid-17/LS/kaplan/users/obrienan/racs-wip/racs-wip-low-1/neighbours/${release_epoch}/${field}
parset=selavy.${field}.${release_epoch}.I.conv.in
log="${working_dir}/selavy.${field}.${release_epoch}.I.conv.log"
slurm_log="${working_dir}/slurm-%A_%a.out"

# echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}, field: ${field}, parset: ${parset}, log: ${log}, working_dir: ${working_dir}"
# echo "Cleaning up /tmp ..."
# find /tmp -maxdepth 1 -user $(whoami) -name "rootfs-*" -exec rm -r {} \;
# echo "Done"

srun \
    --chdir=${working_dir} \
    --mpi=pmi2 \
    --ntasks=${SLURM_NTASKS} \
    --job-name=selavy-${field}-${release_epoch} \
    --output=${slurm_log} \
    --error=${slurm_log} \
        singularity exec \
            --bind /usr/pkg/bin \
            --bind /raid-17/LS/kaplan/users/obrienan/racs-wip/racs-wip-low-1/neighbours \
            --home ${working_dir} \
            $containerImage selavy -c "$parset" >> "$log"
