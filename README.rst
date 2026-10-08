Difference imaging forced photometry pipeline
=============================================

phrosty: **PHotometry for ROman with SFFT for tYpe Ia supernovae**

Aldoroty, L.,  *et al.*, 2026, in prep

This is lightcurve software being developed for the Roman Supernova Project Infrastructure Team (SNPIT).  It is designed to take a sequence of Roman images, divided into "template" and "science" images.  It subtracts each template image from each science image.  It then measures the flux in a point source at a specified RA/Dec on each difference image.  A combination of all of the fluxes measured for a given science image provide an estimate of the flux at the MJD of that science image.  All of these fluxes together provide the lightcurve.

Prerequisites and Environment
-----------------------------

SFFT works in both CPU and GPU environments. In order to use the GPU backend, ``phrosty`` requires a machine with an NVIDIA GPU that has at least 30GB (or more) of GPU RAM.  The GPU nodes at NERSC/Perlmutter (with 40GB of GPU RAM) and SMDC meet this requirement, but consumer graphics cards with only 12GB of RAM aren't sufficient.

Please refer to `the docs <https://roman-supernova-pit.github.io/phrosty/installation.html>`_ for a more detailed description of the prerequisites and setup. 

Running the code
----------------

Please see `the documentation <https://roman-supernova-pit.github.io/phrosty/usage.html>`_ for instructions and usage examples.

License
=======

This project is Copyright (c) Lauren Aldoroty and licensed under
the terms of the BSD-3Clause license. This package is based upon
the `Roman Supernova PIT packaging guide <https://github.com/Roman-Supernova-PIT/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.
