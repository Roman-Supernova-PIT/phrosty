#!/bin/bash
#SBATCH --account m4385
#SBATCH --constraint gpu
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 40G
#SBATCH --gpus-per-task 1
#SBATCH --qos shared
#SBATCH --time 00:30:00
#SBATCH --job-name 20172782_slurm_demo
#SBATCH --output 20172782_slurm_demo.out

echo 'phrosty slurm job starting at' `date`
echo 'phrosty slurm job id ' $SLURM_JOB_ID
echo 'phrosty slurm job running on ' $SLURM_JOB_NODELIST
echo

podman-hpc images
echo

podman-hpc run --gpu \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/phrosty,target=/phrosty \
    --mount type=bind,source=$SCRATCH/phrosty_intermediate,target=/scratch \
    --mount type=bind,source=$PWD/sn_info_dir,target=/sn_info_dir \
    --mount type=bind,source=$SCRATCH/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$SCRATCH/lc_out_dir,target=/lc_out_dir \
    --mount type=bind,source=$SCRATCH/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=/global/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data,target=/sims_dir \
    --mount type=bind,source=/global/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/PQ+HDF5_ROMAN+LSST_LARGE,target=/snana_pq_dir \
    --mount type=bind,source=/global/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/ROMAN+LSST_LARGE_SNIa-normal,target=/snid_lc_dir \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env PYTHONPATH=/roman_imsim:/phrosty:/sfft \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env SIMS_DIR=/sims_dir \
    --env SNANA_PQ_DIR=/snana_pq_dir \
    --env TERM=xterm \
    -w /home \
    -it \
    registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev \
    python phrosty/phrosty/pipeline.py \
      -c phrosty/examples/perlmutter/phrosty_config.yaml \
      -oc ou2024 \
      -ic ou2024 \
      --oid 20172782 \
      -r 11.728845489889757 \
      -d -42.702185894608725 \
      -b R062 \
      -t phrosty/examples/perlmutter/20172782_instances_templates_1.csv \
      -s phrosty/examples/perlmutter/20172782_instances_science.csv \
      -p 9 \
      -w 5 \
      -v
