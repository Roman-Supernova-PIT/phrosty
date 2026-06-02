Difference imaging forced photometry pipeline
=============================================

phrosty: **PHotometry for ROman with SFFT for tYpe Ia supernovae**

Aldoroty, L.,  *et al.*, 2026, in prep

This is lightcurve software being developed for the Roman Supernova Project Infrastructure Team (SNPIT).  It is designed to take a sequence of Roman images, divided into "template" and "science" images.  It subtracts each template image from each science image.  It then measures the flux in a point source at a specified RA/Dec on each difference image.  A combination of all of the fluxes measured for a given science image provide an estimate of the flux at the MJD of that science image.  All of these fluxes together provide the lightcurve.

Statement of Need
-----------------
Difference imaging analysis (DIA) is a technique in transient analysis that allows both transient detection and photometric measurements. Broadly, this method involves subtracting two images, one with a transient and a "reference" image without a transient, to isolate the object of interest. 

The Nancy Grace Roman Space Telescope (_Roman_) is NASA's first survey-oriented flagship mission. Its design is optimized for cosmological studies, which includes extraction and precision analysis of Type Ia Supernova (SN Ia) light curves. Analysis of these data pose several challenges: (1) the enormous quantity of data that will be returned on a regular basis, and (2) the spatially-dependent nature of the image properties as a function of location on the detector, which directly affects ease of (3) obtaining the required &lt;1\% flux precision on SN Ia measurements.

_Roman_'s High Latitude Time Domain Survey (HLTDS) is specifically designed for cosmological studies with SNe Ia. The core component of this survey will last two years, and have a cumulative 158 days of exposure time. The survey will yield an average of 241 Level 2 (rate files) from imaging mode observations per day, which totals to approximately 0.4 TB of data. Other DIA and forced photometry survey pipelines are insufficient to process this quantity of data. For example, the Dark Energy Survey (DES) pipeline, which uses the HOTPANTS algorithm, takes 10 minutes to create a difference image from an image from a single detector, with 8.3 million pixels per detector. Using this pipeline, it would take 80 hours to make difference images from all 241 _Roman_ HLTDS images obtained within one 24 hour period, which have nearly 17 million pixels per detector. In reality, the amount of time spent processing these data will be much more than 80 hours because (1) there are many steps associated with SN Ia light curve analysis in addition to creating the difference images, (2) data from other concurrent surveys will also be processed for transient studies, and (3) this amount of time also excludes re-processing data when pipelines or calibration information are updated.

`phrosty` is a DIA and forced photometry pipeline that addresses these challenges. We integrate the saccadic fast fourier transform (SFFT) (Hu *et al.*, 2022, ApJ, 936, 157) method for difference imaging analysis (DIA), which uses a flexible delta basis function that can accommodate Roman’s spatially-varying point spread function (PSF). `phrosty` uses both GPU and CPU computation; this architecture has been carefully chosen to address the challenge of processing the large quantity of survey data within a reasonable amount of time. 

Because Roman software is currently under-development, there is no similar publicly available all-in-one pipeline that takes _Roman_ images as input and outputs light curves. The architecture choices made in development have substantially accelerated processing time such that `phrosty` is among the fastest astronomical DIA pipelines available. 

Prerequisites and Environment
-----------------------------

Although SFFT works in both CPU and GPU environments, currently some of the code in *phrosty* requires a CUDA-based GPU.  It requires a machine with an NVidia GPU that has at least 20MB (or more) of GPU RAM.  The Perlmutter GPU nodes at NERSC (with 40GB of GPU RAM) meet this requirement, but consumer graphics cards with only 12GB of RAM aren't sufficient.

A number of standard python libraries are required, all of which may be installed from pip, including ``numpy``, ``scipy``, ``astropy``, ``scikit-image``, ``cupy-cuda12x``, ``pandas``, ``requests``, ``fitsio``, ``pyarrow``, ``fastparquet``, and (maybe?) ``nvmath-python[cu12]``.  If you are using OpenUniverse images, *phrosty* also requires ``galsim``, which may be installed from pip; and ``roman_imsim``, which may be git cloned from ``https://github.com/matroxel/roman_imsim.git``.

*phrosty* should run in the `roman snpit environment <https://github.com/Roman-Supernova-PIT/environment>`_.  In particular, it should run in the docker image created from that environment. As of this writing, you can pull the latest version of the docker image from either of:

* ``rknop/roman-snpit-env:cuda-dev``
* ``registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev``

(For an example of actually running phrosty in this dockerized environment, see ``examples/perlmutter/README.md`` and/or [the documentation](https://roman-supernova-pit.github.io/phrosty/).)

Running the code
----------------

Please see the [documentation](https://roman-supernova-pit.github.io/phrosty/).

Examples
--------

See [the documentation](https://roman-supernova-pit.github.io/phrosty/). Members of the PIT, if you have any issues running phrosty, ping Rob Knop and Lauren Aldoroty in the ``#photometry`` channel on the SN PIT Slack.

Running tests
-------------

The tests need to be run within an environment where *phrosty* will run properly.  At the moment, the only environment we've successfully used this with is the Docker environment described in the NERSC/Perlmutter environment above.  Follow that example through running ``bash interactive_podman.sh``.  Then run::

  cd /home/phrosty/tests
  export DIA_OUT_DIR=../../dia_out_dir
  export SN_INFO_DIR=../../sn_info_dir
  pytest -v

License
=======

This project is Copyright (c) Lauren Aldoroty and licensed under
the terms of the BSD-3Clause license. This package is based upon
the `Roman Supernova PIT packaging guide <https://github.com/Roman-Supernova-PIT/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.

