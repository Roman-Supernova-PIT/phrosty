Difference imaging forced photometry pipeline
=============================================

phrosty: **PHotometry for ROman with SFFT for tYpe Ia supernovae**

Aldoroty, L.,  *et al.*, 2026, in prep

This is lightcurve software being developed for the Roman Supernova Project Infrastructure Team (SNPIT).  It is designed to take a sequence of Roman images, divided into "template" and "science" images.  It subtracts each template image from each science image.  It then measures the flux in a point source at a specified RA/Dec on each difference image.  A combination of all of the fluxes measured for a given science image provide an estimate of the flux at the MJD of that science image.  All of these fluxes together provide the lightcurve.

Prerequisites and Environment
-----------------------------

Although SFFT works in both CPU and GPU environments, currently some of the code in *phrosty* requires a CUDA-based GPU.  It requires a machine with an NVidia GPU that has at least 20MB (or more) of GPU RAM.  The Perlmutter GPU nodes at NERSC (with 40GB of GPU RAM) meet this requirement, but consumer graphics cards with only 12GB of RAM aren't sufficient.

A number of standard python libraries are required, all of which may be installed from pip, including ``numpy``, ``scipy``, ``astropy``, ``scikit-image``, ``cupy-cuda12x``, ``pandas``, ``requests``, ``fitsio``, ``pyarrow``, ``fastparquet``, and (maybe?) ``nvmath-python[cu12]``.  If you are using OpenUniverse images, *phrosty* also requires ``galsim``, which may be installed from pip; and ``roman_imsim``, which may be git cloned from ``https://github.com/matroxel/roman_imsim.git``. *However*, if you use the provided containerized environments, the prerequisite packages have been dealt with for you.

*phrosty* should run in the `roman snpit environment <https://github.com/Roman-Supernova-PIT/environment>`_.  In particular, it should run in the docker image created from that environment. As of this writing, you can pull the latest version of the docker image from either of:

* ``rknop/roman-snpit-env:cuda-dev``
* ``registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev``

(For an example of actually running phrosty in this dockerized environment, see ``examples/perlmutter/README.md`` and/or `the documentation <https://roman-supernova-pit.github.io/phrosty/>`_.)

Running the code
----------------

Please see `the documentation <https://roman-supernova-pit.github.io/phrosty/>`_.

Examples
--------

See `the documentation <https://roman-supernova-pit.github.io/phrosty/>`_. Members of the PIT, if you have any issues running phrosty, ping Rob Knop and Lauren Aldoroty in the ``#photometry`` channel on the SN PIT Slack.

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

