.. highlight:: shell

.. _phrosty-installation:

============
Installation
============

.. contents::

.. _system-requirements:

System Requirements
-------------------

``phrosty`` can run using either a ``cupy`` (CUDA 12.4, requires an NVIDIA GPU) or ``numpy`` backend (CPU). Empirically, you will need at least 36 GB GPU memory or 56 GB CPU memory to run these backends, respectively, for a standard 4088 x 4088 px *Roman* image.

To properly set up ``phrosty``, you need to follow **one** of the sections here, followed by `Install from sources<install-from-sources>`.

.. _phrosty-local:

Running locally
---------------

.. _phrosty-local-common-setup:

Setting up directories and getting the code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You will want to make yourself and environment to work in, such as a conda environment or python venv, as you will need to pip install some packages.  Alternatively, you can work in a docker container.

If you're using conda or a venv
"""""""""""""""""""""""""""""""

Go into your environment and install the ``snappl`` and ``sfft`` packages that ``phrosty`` needs:

.. _code-block: console

  pip install roman-snpit-snappl sfft-romansnpit crds

If you're using a docker container
""""""""""""""""""""""""""""""""""

Depending on whether you're going to run the CPU or GPU version of ``phrosty`` (read below) you will want to pull one of two docker images:

  * ``docker pull docker.io/rknop/roman-snpit-env:cpu-dev``
  * ``docker pull docker.io/rknop/roman-snpit-env:cuda-dev``

Make directories
""""""""""""""""

Pick a directory to work in; we shall call that ``$RUNDIR``; make sure that's your current directory.

You will need some directories for phrosty to work in:

 * ``mkdir -p packages`` : you will do all git cloning in this subdirectory
 * ``mkdir -p temp_dir``: This is where temporary files will be written.
 * ``mkdir -p dev_storage``: A general directory for storing files that you might want to keep around a while
 * ``mkdir -p dev_storage/dia_out_dir``: Output image files are written here.
 * ``mkdir -p dev_storage/ltcv_dir``: Output lightcurves are written here.
 * ``mkdir -p dev_storage/intermediate_dir``: Intermediate files are written here.

In ``$RUNDIR/packages``, check out the ``phrosty`` archive.  (Eventually, ``phrosty`` will be on pip. As of writing this, it is not.)  If you're going to run tests and/or the examples in :ref:`usage`, also check out the ``photometry_test_data`` archive::

.. _code-block: console

  cd packages
  git clone https://github.com/Roman-Supernova-PIT/phrosty.git
  git clone https://github.com/Roman-Supernova-PIT/photometry_test_data.git
  cd ..

(TODO: figure out if there's a ``git-lfs`` thing people have to do.)

You will need a standard default config file, which is referenced by the ``phrosty`` config files used in the examples below.  Assuming you are still in ``$RUNDIR``, run:

.. _code_block: console

  curl -L https://raw.githubusercontent.com/Roman-Supernova-PIT/environment/refs/heads/main/local_nodb.yaml -O
  curl -L https://raw.githubusercontent.com/Roman-Supernova-PIT/environment/refs/heads/main/container_nodb.yaml -O

This will copy down a standard Roman SNPIT/snappl config file, which you use when you don't want to connect to any database.


.. _phrosty-local-cpu:

If you are going to run in a virtual/conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, make sure you are in ``$RUNDIR``, and run:

.. _code-block: console

  cd packages/phrosty
  pip install .

(Note: if you're developing ``phrosty``, you might want that last line to be ``pip install -e .``.  If you want to run the ``phrosty`` tests, do ``pip install -e .[test]``.)

NOTE: Eventually, ``phrosty`` will be on pip. As of writing this, it is not.  When that happens, you can replace the last three lines with just ``pip install roman-snpit-phrosty``.

You need to set two environment variables to tell ``phrosty`` where to find the config files:

.. _code_block: console

   export SNPIT_DEFAULT_CONFIG=${PWD}/local_nodb.yaml
   export SNPIT_CONFIG=${PWD}/packages/phrosty/phrosty_config_default.yaml

Finally, also set the following environment variables:

.. _code_block: console

  export CRDS_SERVER_URL=https://roman-crds.stsci.edu
  export CRDS_PATH=${HOME}/crds_cache

You can make ``CRDS_PATH`` exist in ``$RUNDIR`` instead if you want.

You should be good to go now.

If you have a 40GB NVIDIA GPU
"""""""""""""""""""""""""""""

(...and if you want to use it rather than running on the cpu...)

You will also need to install ``cupy``.  This can be challenging, and there may be issues of getting versions of ``cupy`` that are consistent with the NVIDIA drivers and CUDA version installed on your system.  You can try:

.. _code-block: console

  pip install cupy-cuda12x

but you may find that you need a different CUDA version. For more help with matching your CUDA version to your ``cupy`` version, see `the cupy documentation<https://docs.cupy.dev/en/stable/install.html>`_.`


If you are going to run in a docker container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You go into the container by running:

.. _code-block: console

  docker run -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/packages,target=/packages \
    --mount type=bind,source=$PWD/temp_dir,target=/temp_dir \
    --mount type=bind,source=$PWD/dev_storage,target=/dev_storage \
    --mount type=bind,source=$PWD/packages/photometry_test_data,target=/photometry_test_data \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    --env CRDS_SERVER_URL=https://roman-crds.stsci.edu \
    --env CRDS_PATH=/home/crds_cache \
    --env SNPIT_DEFAULT_CONFIG=/home/container_nodb.yaml \
    --env SNPIT_CONFIG=/packages/phrosty/phrosty_default_config.yaml \
    --annotation run.oci.keep_original_groups=1 \
    rknop/roman-snpit-env:cpu-dev \
    /bin/bash

If all is well, this will put you in a docker container.  You can tell you're in the container because your prompt will change to something like ``root@47394bd41fbe:/#`` (where the string of hexidecimal numbers will be different every time you start a container).  Your ``$RUNDIR`` is mounted at ``/home`` inside the container.

Next, you will want to get ``phrosty`` installed inside your environment.  (You will need to do this every time you restart the container.)

.. _code-block: console

  cd /packages/phrosty
  pip install .
  cd /home

If you are developing ``phrosty``, you might want to do ``pip install -e .``, and if you think you might want to run the tests, you might want to do ``pip install -e .[test]``.

At this point, you should be good to go.

When you're done with your container, you can just ``exit`` to get out of it.  You may also want to do ``docker ps`` followed by ``docker rm <container-id>`` to clean up cruft left behind on your system.  (You can find ``<container-id>`` by looking at the output of ``docker ps``.)

If you have a 40GB NVIDIA GPU and want to use it
"""""""""""""""""""""""""""""""""""""""""""""""""

Add ``--gpus=all`` between ``docker run`` and ``-it``.  Also, replace ``rknop/roman-snpit-env:cpu-dev`` with ``rknop-snpit-env:cuda-dev``.  This can be fraught; the versions of the NVIDIA drivers you have on your system have to be compatible with what's inside the container.  Once you're inside the container, verify that you can see your GPU with:

.. _code-block: console

   nvidia-smi


.. _phrosty-smdc:

Installing on SMDC
^^^^^^^^^^^^^^^^^^

**This section will work for SN PIT members.**

General instructions for accessing SMDC can be found `in the wiki <https://github.com/Roman-Supernova-PIT/Roman-Supernova-PIT/wiki/NASA-SMDC-%28AWS%29>`_.

There is some information at `"Working with PIT Images" here <https://github.com/Roman-Supernova-PIT/Roman-Supernova-PIT/wiki/SMCE-Containers>`_, which is largely applicable if you want to use the Singularity containers; however, if all you want to do is go into a standard SNPIT container, the instructions linked below have everything you need.

First, ``salloc`` a node. If you want a GPU node, do::

  salloc -p gpu-int --time=04:00:00

If you want to run on a CPU node, do::

  salloc -p mem-lg --time=04:00:00

Sometimes your correct set of groups won't be correctly populated on a compute node due to a race condition between populating the container and correctly configuring the active directory lookup. You will see a message about this when you get your node that says that the groups weren't loaded properly. Also, if you type ``groups`` on the login node, you'll see ``[your username] spack cluster_users snpit``. If you type ``groups`` on the GPU node, you'll see ``[your username] nogroup``. You will also hit a permissions issue running ``phrosty`` when it tries to write files outside your home directory. (In particular, you need to be in the ``snpit`` group.)  To start a new terminal that will have the groups loaded correctly, do::

  ssh localhost

Then, follow the instructions in `the snappl documentation about running on SMDC <https://roman-supernova-pit.github.io/snappl/environment.html#running-on-smdc>`_.

.. _phrosty-nersc-perlmutter:

Installing on NERSC Perlmutter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**This section will work for SN PIT members, and maybe anyone else with access to NERSC Perlmutter, which ostensibly could be you.**

First, follow the instructions in `the snappl docs<https://roman-supernova-pit.github.io/snappl/environment.html#running-on-nersc>`_. If you aren't in the SN PIT, skip the stuff about a password and "secrets" folder to use the database.

If you are on NERSC Perlmutter, you have access to NVIDIA GPUs with 40 GB GPU memory. Run ``module list``.  Make sure that ``cudatoolkit/12.4`` shows up in your list of modules.  If not, you may need to adjust the modules you have loaded.

(Ideally, because ``phrosty`` runs inside a container, the specific version of CUDA on the host system wouldn't matter.  However, containers can be touchy about hooking up GPUs inside containers.  If you're having trouble seeing the GPU inside the container, check for mismatches between the versions of CUDA inside and outside of the container.)

If you're on NERSC Perlmutter, use ``podman-hpc`` in place of ``docker``.  Pull the image with::

  podman-hpc pull registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev

**Note:** do *not* run ``podman-hpc image pull ...``.  That will superficially seem to work, but will skip a step that gets run when you just do ``podman-hpc pull ...``.

If you get a permission error trying to do this, try::

  podman-hpc login registry.nersc.gov

Give it your usual NERSC username and password (without any OTP).  Once that's done, try the ``podman-hpc pull`` command again.  If you don't seem to have access to the registry, then you can just pull ``docker.io/rknop/roman-snpit-env:cuda-dev`` instead.

After you've pulled, run ``podman-hpc images``.  You should see output something like::

  REPOSITORY                                          TAG                 IMAGE ID      CREATED         SIZE        R/O
  registry.nersc.gov/m4385/rknop/roman-snpit-env      cuda-dev            6b39a47ffc5b  25 minutes ago  8.6 GB      false
  registry.nersc.gov/m4385/rknop/roman-snpit-env      cuda-dev            6b39a47ffc5b  25 minutes ago  8.6 GB      true

In particular, notice that there is both a R/O "true" and R/O "false" version.  The R/O "true" version is the one you need to be able to run on nodes other than the node you pulled the image (such as compute nodes).  (The Image ID will probably not match what you see above, if we're released version of the snpit docker image since this documentation was written.  The "CREATED" time is also very likely to be longer ago than what you see here.)

If you've pulled images before, and you're now working on a new login node, you will only see the ``R/O=true`` image.  That's the only image you really need, so in that case there's no need to pull the image again.  (You will only see the ``R/O=true`` image on compute nodes.)

**If you have trouble with podman**: Refer to `NERSC's documentation on podman-hpc <https://docs.nersc.gov/development/containers/podman-hpc/overview/>`_.  In particular, if you want to clean the slate and start over, try running::

  podman-hpc system reset

to delete all of your podman images and contexts.  Then try pulling the image again.

Assuming you're in the directory above your ``phrosty`` and ``photometry_test_data`` checkouts, you can run the container with ``bash /global/cfs/cdirs/m4385/env/interactive-podman-nov2025.sh``. At this time, both of these files are the same, but you have the ability to modify the one in ``examples/perlmutter`` and not the one in ``m4385/env``.

If you absolutely must make your own container for some reason, see `the interactive podman scripts in our environment directory <https://github.com/Roman-Supernova-PIT/environment/blob/main/interactive-podman-nov2025.sh>`_ for reference.

If you're inside the container, your prompt will be something like ``root@f24c2ad04d6d:/#`` (though with a different string of hexidecimal digits (hexits?)).  If you do ``ls -F /``, you will see the various specific directories that are mounted in the above command.

Verify that you have access to GPUs by running::

  nvidia-smi

.. _photometry-test-data:

Installing the photometry test data (recommended but optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to run tests, and some of the examples, then you will also need to pull the photometry test data into ``$RUNDIR``::

  git clone https://github.com/Roman-Supernova-PIT/photometry_test_data.git

.. _Github repo: https://github.com/Roman-Supernova-PIT/phrosty
.. _tarball: https://github.com/Roman-Supernova-PIT/phrosty/tarball/master
