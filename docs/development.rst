.. highlight:: shell

===========
Development
===========

If you're one of the phrosty developers, or are otherwise interested in contributing, here is some useful information.

.. _running-tests:
Running Tests
-------------

**Warning**: all docker images are currently built only for the ``amd64`` architecture (also sometimes known as ``x86_64``).  If you're on ``arm64``, then things may not work, and if they do work, they may be horribly slow as you're emulating a different architecture.  NERSC's Perlmutter cluster, and any Linux system running on hardware with an Intel or AMD CPU, are on the ``amd64`` architecture.  Current Macs, and any other systems based on an ARM chip, are on the ``arm64`` architecture.

Running all of the tests requires an NVIDIA GPU with enough memory.  We are able to run them on 40GB NVIDIA GPUs, a GPU with only 12GB is not enough.  (TODO: figure out the actual cutoff.)

To run the tests, you will need to run the SNPIT container as briefly described in :ref:`running-snpit-container`.  However, you must also make sure to mount the ``photometry_test_data`` archive in your container.

First, if you haven't already, get a copy of phrosty::
  git clone https://github.com/Roman-Supernova-PIT/phrosty.git

Second, in the same directory, get a copy of the photometry test datay::
  git clone https://github.com/Roman-Supernova-PIT/photometry_test_data.git

Next, the latest ``cuda-dev`` image to your local system.  If you're using docker, do this with one of::
  docker pull registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev
  docker pull docker.io/rknop/roman-snpit-env:cuda-dev

Run the container with::
  docker run -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    rknop/roman-snpit-env:cuda-dev \
    /bin/bash

If you have access to the ``registry.nersc.gov`` container repository, use that image instead of ``rknop/roman-snpit-env:cuda-dev``.  If you're on Perlmutter, replace ``docker`` with ``podman-hpc``.

Once inside the container, cd into ``/home/phrosty/phrosty/tests`` and run::
  SNPIT_CONFIG=phrosty_test_config.yaml pytest -v

Manually running a test lightcurve
----------------------------------

Currently, we do not have tests written for ``phrosty/pipeline.py``; this will be rectified soon, and hopefully this documentation will be updated to reflect that.

If you want to build a lightcurve using the test data, follow the instructions in the previous section for getting the ``phrosty`` and ``photometry_test_data`` archives, and for pulling the ``cuda-dev`` version of the roman SNPIT docker image.

You need a few of extra directories when running your container, to store temporary and output data.  You can put these where you want, but for this example we are going to assume you make three directories underneath the same place you checked out the two git archvies: ``phrosty_temp``, ``dia_out_dir``, and ``lc_out_dir``.

With these directories in place, run a container with::
  docker run -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    --mount type=bind,source=$PWD/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/lc_out_dir,target=/lc_out_dir \
    rknop/roman-snpit-env:cuda-dev \
    /bin/bash

If you placed any of the new directories anywhere other than underneath your current working directory, modify the ``source=...`` parts of the command above to reflect that.  If you have access to the ``registry.nersc.gov`` repostitory, replace ``rknop/roman-snpit-env:cuda-dev`` with ``registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev``.  On Perlmutter, replace ``docker`` with ``podman-hpc``.
    
Inside the container, cd into ``/home/phrosty`` and try running:

