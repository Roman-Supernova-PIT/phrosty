podman-hpc run --gpu \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$HOME/secrets,target=/secrets \
    --mount type=bind,source=$PSCRATCH/snpit_temp,target=/snpit_temp \
    --mount type=bind,source=$PSCRATCH,target=/scratch \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/m4385/env,target=/snpit_env \
    --mount type=bind,source=/pscratch/sd/m/masao/roman_snpit,target=/roman_snpit_masao_scratch \
    --mount type=bind,source=/pscratch/sd/m/masao/roman_snpit/database_dirs,target=/data \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data,target=/ou2024 \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/PQ+HDF5_ROMAN+LSST_LARGE,target=/ou2024_snana \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/ROMAN+LSST_LARGE_SNIa-normal,target=/ou2024_snana_lc_dir \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/sims_sed_library,target=/ou2024_sims_sed_library \
    --mount type=bind,source=/dvs_ro/cfs/cdirs/m4385/calib_data/A25ePSF,target=/a25epsf \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env PYTHONPATH=/roman_imsim \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    --env SNPIT_DEFAULT_CONFIG=/snpit_env/configs/nov2025_container_config.yaml \
    --env SNPIT_CONFIG=/snpit_env/configs/nov2025_container_config.yaml \
    --annotation run.oci.keep_original_groups=1 \
    -it \
    registry.nersc.gov/m4385/roman-snpit-env:cuda-dev-0.1.41 \
    /bin/bash