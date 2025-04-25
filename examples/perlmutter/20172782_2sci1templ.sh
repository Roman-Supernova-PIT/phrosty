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
    --mount type=bind,source=$PWD/sfft,target=/sfft \
    --mount type=bind,source=$PWD/sn_info_dir,target=/sn_info_dir \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/lc_out_dir,target=/lc_out_dir \
    --mount type=bind,source=$SCRATCH/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data,target=/sims_dir \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/PQ+HDF5_ROMAN+LSST_LARGE,target=/snana_pq_dir \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/ROMAN+LSST_LARGE_SNIa-normal,target=/snid_lc_dir \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env PYTHONPATH=/roman_imsim:/phrosty:/sfft \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env SIMS_DIR=/sims_dir \
    --env SNANA_PQ_DIR=/snana_pq_dir \
    --env SN_INFO_DIR=sn_info_dir \
    --env DIA_OUT_DIR=dia_out_dir \
    --env TERM=xterm \
    -w /home \
    -it \
    registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev \
    python phrosty/phrosty/pipeline.py \
      --oid 20172782 \
      -r 7.551093401915147 \
      -d -44.80718106491529 \
      -b R062 \
      -t 20172782_instances_templates_1.csv \
      -s 20172782_instances_science_two_images.csv \
      -p 15 \
      -w 5
