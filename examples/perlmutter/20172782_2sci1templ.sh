#!/bin/bash
#SBATCH --account m4385
#SBATCH --constraint gpu
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
###SBATCH --mem 60G
#SBATCH --gpus-per-task 1
#SBATCH --qos shared
#SBATCH --time 00:10:00
#SBATCH --output 20172782_2sci1temp.out

echo 'phrosty slurm job starting at' `date`
echo 'phrosty slurm job id ' $SLURM_JOB_ID
echo 'phrosty slurm job running on ' $SLURM_JOB_NODELIST
echo

podman-hpc run --gpu \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/phrosty,target=/phrosty \
    --mount type=bind,source=$PWD/sn_info_dir,target=/sn_info_dir \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/lc_out_dir,target=/lc_out_dir \
    --mount type=bind,source=$SCRATCH/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data,target=/ou2024 \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/PQ+HDF5_ROMAN+LSST_LARGE,target=/ou2024_snana \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/ROMAN+LSST_LARGE_SNIa-normal,target=/ou2024_snana_lc_dir \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/sims_sed_library,target=/ou2024_sims_sed_library \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env PYTHONPATH=/phrosty \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    -w /home \
    -it \
    registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev \
    python phrosty/phrosty/pipeline.py \
      -c phrosty/examples/perlmutter/phrosty_config.yaml \
      --oid 20172782 \
      -r 7.551093401915147 \
      -d -44.80718106491529 \
      -b R062 \
      -t phrosty/examples/perlmutter/20172782_instances_templates_1.csv \
      -s phrosty/examples/perlmutter/20172782_instances_science.csv \
      -p 15 \
      -w 5
