.. highlight:: shell

===========
Development
===========

If you're one of the phrosty developers, or are otherwise interested in contributing, here is some useful information.

Note that fully running phrosty requires an NVIDIA GPU with enough memory.  Empirically, as 12GB GPU is not enough.  We have succesfully run phrosty on NVIDIA GPUs with 40GB of memory.

.. _running-tests:

Running Tests
-------------

**Warning**: all docker images are currently built only for the ``amd64`` architecture (also sometimes known as ``x86_64``).  If you're on ``arm64``, then things may not work, and if they do work, they may be horribly slow as you're emulating a different architecture.  NERSC's Perlmutter cluster, and any Linux system running on hardware with an Intel or AMD CPU, are on the ``amd64`` architecture.  Current Macs, and any other systems based on an ARM chip, are on the ``arm64`` architecture.

Running all of the tests requires an NVIDIA GPU with enough memory.  We are able to run them on 40GB NVIDIA GPUs, a GPU with only 12GB is not enough.  (TODO: figure out the actual cutoff.)

To run  the tests, forst make sure you've set up your environment and pulled down the necessary docker images as described in :ref:`phrosty installation prerequisites<phrosty-installation-prerequisites>`.  

If you haven't already, get a copy of phrosty::

  git clone https://github.com/Roman-Supernova-PIT/phrosty.git

Second, in the same directory, get a copy of the photometry test datay::

  git clone https://github.com/Roman-Supernova-PIT/photometry_test_data.git

Make a couple of necessary directories::

  mkdir dia_out_dir
  mkdir phrosty_temp
  
Run the container with::

  docker run --gpus=all -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/phrosty_temp,target=/phrosty_temp \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env PYTHONPATH=/roman_imsim \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    rknop/roman-snpit-env:cuda-dev \
    /bin/bash

**On NERSC Perlmutter**, run the container with::

  podman-hpc run --gpu -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/phrosty_temp,target=/phrosty_temp \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env PYTHONPATH=/roman_imsim \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    --annotation run.oci.keep_original_groups=1 \
    registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev \
    /bin/bash

Once inside the container, cd into ``/home/phrosty/phrosty/tests`` and run::

  SNPIT_CONFIG=phrosty_test_config.yaml pytest -v


Manually running a test lightcurve
------------------------------------

Currently, we do not have tests written for ``phrosty/pipeline.py``; this will be rectified soon, and hopefully this documentation will be updated to reflect that.

If you want to build a lightcurve using the test data, follow the instructions in the previous section for getting the ``phrosty`` and ``photometry_test_data`` archives, and for pulling the ``cuda-dev`` version of the roman SNPIT docker image.

You need a few of extra directories when running your container, to store temporary and output data.  You can put these where you want, but for this example we are going to assume you make three directories underneath the same place you checked out the two git archvies: ``phrosty_temp``, ``dia_out_dir``, and ``lc_out_dir``.

With these directories in place, run a container with::

  docker run --gpus=all -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    --mount type=bind,source=$PWD/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/lc_out_dir,target=/lc_out_dir \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env PYTHONPATH=/roman_imsim \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    rknop/roman-snpit-env:cuda-dev \
    /bin/bash

**on NERSC Perlmutter**, the command would be::

  podman-hpc run -gpu -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    --mount type=bind,source=$PWD/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/lc_out_dir,target=/lc_out_dir \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env PYTHONPATH=/roman_imsim \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    --annotation run.oci.keep_original_groups=1 \
    registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev \
    /bin/bash

If you placed any of the new directories anywhere other than underneath your current working directory, modify the ``source=...`` parts of the command above to reflect that.
    
Inside the container, cd into ``/home/phrosty`` and try running::

  nvidia-smi

If you don't get errors, it should list the nvidia GPUs you have available.  If it doesn't list GPUs, then the rest of this won't work.

Next, try running::

  cd /home/phrosty
  pip install -e .
  SNPIT_CONFIG=phrosty/tests/phrosty_test_config.yaml python phrosty/pipeline.py --help | less

You should see all the options you can pass to phrosty.  There are a lot, because there are (verbose) options for everything that's in the config file.  Press ``q`` to get out of ``less``.

Try running::

  SNPIT_CONFIG=phrosty/tests/phrosty_test_config.yaml python phrosty/pipeline.py \
    --oid 20172782 \
    --ra 7.551093401915147 \
    --dec -44.80718106491529 \
    -b Y106 \
    -t phrosty/tests/20172782_instances_templates_1.csv \
    -s phrosty/tests/20172782_instances_science_2.csv \
    -p 3 -w 3 \
    -v

If all is well, after it's done running the output will end with something like::

  [2025-08-13 17:35:24 - INFO] - Results saved to /lc_out_dir/data/20172782/20172782_Y106_all.csv

On your host system (as well as inside the container), you should see new files in ``lc_out_dir``, ``dia_out_dir``, and ``phrosty_temp``.  (Inside the container, these are at ``/lc_out_dir``, ``/dia_out_dir``, and ``/phrosty_temp``.)
