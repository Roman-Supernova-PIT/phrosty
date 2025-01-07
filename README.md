# phrosty

**PHotometry for ROman Simulations protoTYpe.**

Aldoroty, L.,  *et al.*, 2025, in prep

This is lightcurve software being developed for the Roman Supernova Project Infrastructure Team (SNPIT).  It is designed to take a sequence of Roman images, divided into "template" and "science" images.  It subtracts each template image from each science image.  It then measures the flux in a point source at a specified RA/Dec on each difference image.  A combination of all of the fluxes measured for a given science image provide an estimate of the flux at the MJD of that science image.  All of these fluxes together provide the lightcurve.

Currently, *phrosty* only handles the OpenUniverse sims [ref needed].  It reads images from those sims in FITS format, using the standard FITS WCSes found in the headers of those files for image alignment.  The galsim package and a particular version of the roman imsim package (https://github.com/matroxel/roman_imsim.git) provide the PSF used both for PSF photometry and for the PSF matching in the SFFT image subtraction software (Hu *et al.*, 2022, ApJ, 936, 157).

## Preqrequisites and Environment

Although SFFT works in both CPU and GPU environments, currently some of the code in *phrosty* requires a CUDA-based GPU.  It requires a machine with an NVidia GPU that has at least 20MB (or more) of GPU RAM.  The Perlmutter GPU nodes at NERSC (with 40GB of GPU RAM) meet this requirement, but consumer graphics cards with only 12GB of RAM aren't sufficient.

A number of standard python libraries are required, all of which may be installed from pip, including `numpy`, `scipy`, `astropy`, `scikit-image`, `cupy-cuda12x`, `pandas`, `requests`, `fitsio`, `pyarrow`, `fastparquet`, and (maybe?) `nvmath-python[cu12]`.  *Phrosty* also requires `galsim`, which may be installed from pip.

Phrosty requires commit `74a9053` from `roman_imsim`, which may be git cloned from https://github.com/matroxel/roman_imsim.git.  (This is nominally `v2.0` of that archive, but as of this writing, there is both a branch and a tag `v2.0`, which confuses the `git archive` command.  Thus, we give you the specific commit.)

The fledgling/proposed Roman SNPIT docker environment at https://github.com/Roman-Supernova-PIT/snpit_docker_env provides a Dockerized environment that *phrosty* should be able to run in.  As of this writing, we're using version `v0.0.1` of that environment.  *phrosty* currently requires the (more bloated) dev environment, but we hope in the future to make it run with the regular environment.  NOTE: this environment will be deprecated in the near future, in favor of a Docker environment created under https://github.com/Roman-Supernova-PIT/environment.  Once that is available, this README will be updated.  You can pull this docker image from the following locations:

* `docker.io/rknop/roman-snpit-dev:v0.0.1`
* `registry.nersc.gov/m4385/rknop/roman-snpit-dev:v0.0.1`

### Manual installs

Not included in the docker image describewd above is SFFT— because, as of this writing, the verison used with phrosty was under development.  SFFT must be separately cloned from git@github.com:Roman-Supernova-PIT/sfft.git, the directory to which you clone that archive must by in your `PYTHONPATH`.  As of this writing, *phrosty* requires the `fixes_20241022` branch, though hopefully we will get that branch merged back to the `master` branch.  See "Examples" below for samples of this in action.

### Necessary directories data files

Currently, *phrosty* depends on the following environment variables to find data:

* `SIMS_DIR` : OpenUniverse roman data (on nersc: `/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data`)
* `SN_INFO_DIR` : A directory that holds information about the OpenUniverse sims, and transients from the OpenUniverse sims.  It must have a file `tds.yaml` that is properly configured for whatever environment you're running.
* `SNANA_PQ_DIR` : A directory with the transient parquet files from the OpenUniverse sims (on nersc: /dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/PQ+HDF5_ROMAN+LSST_LARGE)
* `DIA_OUT_DIR` : A directory where difference images (and psfs and other associated images) will be written.


## Running the code

TODO.

For now, see the one extant example below.

## Examples


### NERSC/Perlmutter

An example running *phrosty* on an interactive node, and on a slurm-allocated node, on Perlmutter at NERSC may be found in the examples/perlmutter directory.

## Running tests

The tests need to be run within an environment where phrosty will run properly.  At the moment, the only environment we've succesfully used this with is the Docker environment described in the NERSC/Perlmutter environment above.  Follow that example through running `bash interactive_podman.sh`.  Then, instead of doing what the example says, run
```
cd /home/phrosty/tests
export DIA_OUT_DIR=../../dia_out_dir
export SN_INFO_DIR=../../sn_info_dir
pytest -v
```

TODO: we need to clean this up and make the test environment more self-contained.

## Plotting lightcurves

TODO.

For now, you can use your favorite plotting program and just read the `.csv` files produced by the pipeline.  However, we need to document how you read in the OpenUniverse truth files for plotting lightcurves against the truth files.
