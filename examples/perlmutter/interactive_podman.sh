#!/bin/sh

podman-hpc run --gpu \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/phrosty,target=/phrosty \
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
    --annotation run.oci.keep_original_groups=1 \
    -it \
    registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev \
    /bin/bash
