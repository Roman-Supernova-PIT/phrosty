.. highlight:: shell

.. _phrosty-installation:

============
Installation
============

.. _system-requirements:

System Requirements
-------------------

Phrosty is built using CUDA 12.4 and requires an NVIDIA GPU to run.  Your host system needs to have an NVIDIA GPU, and needs to have a version of CUDA compatible with 12.4.  (It's *possible* that phrosty will run with other versions of CUDA on the host system, but we haven't verified this.)  To run all of phrosty, you need a GPU with enough memory.  Empirically, 12GB is not enough GPU memory to fully run phrosty on OpenUniverse simulations of Roman data (though you will be able to run some of the tests).  We have succesfully run phrosty on an NVIDIA GPU with 40GB of memory.

If you are on NERSC Perlmutter, run ``module list``.  Make sure that ``cudatoolkit/12.4`` shows up in your list of modules.  If not, you may need to adjust the modules you have loaded.

(Ideally, because phrosty runs inside a container, the specific version of CUDA on the host system wouldn't matter.  However, containers can be touchy about hooking up GPUs inside containers.  If you're having trouble seeing the GPU inside the container, check for mismatches between the versions of cuda inside and outside of the container.)

.. _phrosty-installation-prerequisites:

Prerequisites
-------------

Currently, phrosty is designed inside a container built from the the `Roman Supernova PIT environment <https://github.com/Roman-Supernova-PIT/environment>`_.

You can pull the necessary container image from one of the following two sources:

* ``registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev``
* ``docker.io/rknop/roman-snpit-env:cuda-dev``

The latter image location should work for everybody.  For members of the Roman SN PIT, the first image location should be more efficient when you're running on NERSC Perlmutter; see below.the image, 

Because phrosty (and other libraries it depends on, such as snappl) is under heavy development, it's possible that the latest container will not work properly with phrosty at any given moment.

On NERSC Perlmutter
^^^^^^^^^^^^^^^^^^^

If you're on NERSC Perlmutter, use ``podman-hpc`` in place of ``docker``.  Pull the image with::
  podman-hpc pull registry.nersc.gov/m4385/rknpo/roman-snpit-env:cuda-dev

**Note:** do _not_ run ``podman-hpc image pull ...``.  That will superficially seem to work, but will skip a step that gets run when you just do ``podman-hpc pull ...``.
  
If you get a permission error trying to do this, try::
  podman-hpc login registry.nersc.gov

Give it your usual NERSC username and password (without any OTP).  Once that's done, try the ``podman-hpc pull`` command again.  If you don't seem to have access to the registry, then you can just pull ``docker.io/rknop/roman-snpit-env:cuda-dev`` instead.

After you've pulled, run ``podman-hpc images``.  You should see output something like::
  REPOSITORY                                          TAG                 IMAGE ID      CREATED         SIZE        R/O
  registry.nersc.gov/m4385/rknop/roman-snpit-env      cuda-dev            6b39a47ffc5b  25 minutes ago  8.6 GB      false
  registry.nersc.gov/m4385/rknop/roman-snpit-env      cuda-dev            6b39a47ffc5b  25 minutes ago  8.6 GB      true

In particular, notice that there is both a R/O "true" and R/O "false" version.  The R/O "true" version is the one you need to be able to run on nodes other than the node you pulled the image (such as compute nodes).  (The Image ID will probably not match what you see above, if we're released version of the snpit docker image since this documentation was written.  The "CREATED" time is also very likely to be longer ago than what you see here.)

If you've pulled images before, and you're now working on a new login node, you will only see the ``R/O=true`` image.  That's the only image you really need, so in that case there's no need to pull the image again.  (You will only see the ``R/O=true`` image on compute nodes.)

**If you have trouble with podman**: Podman is, as of this writing, still a beta feature on Nersc, and has some warts.  Refer to `NERSC's documentation on podman-hpc <https://docs.nersc.gov/development/containers/podman-hpc/overview/>`_.  In particular, if you want to clean the slate and start over, try running::
  podman-hpc system reset
  
to delete all of your podman images and contexts.  Then try pulling the image again.

On other HPC Systems
^^^^^^^^^^^^^^^^^^^^

If your HPC system doesn't support containers, then you're out of luck.  Most HPC systems support containers using ``apptainer`` (previously known as ``singularity``).  This system is not a drop-in replacment for docker, as the way you obtain images, and the semantics for running containers, are different.  (There are also differences in terms of how isolated the environment is; singularity is less of a "real" isolated container than docker, usually.)

TODO : document use of singularity.


Stable release
--------------

phrosty is under heavy development and does not currently have stable releases.  Eventually it will be released on pypi either under the name ``phrosty`` or ``roman-snpit-phrosty`` (depending on what's available).

.. _install-from-sources:

From sources
------------

Currently, the only way to install phrosty is to download it from the `github repo <https://github.com/Roman-Supernova-PIT/phrosty>`_.  Clone it with::

    git clone https://github.com/Roman-Supernova-PIT/phrosty.git

(you can also clone it via the ``git@`` code link if you know what you're doing.)

Installing the photometry test data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to :ref:`run tests<running-tests>`, and some of the examples, then you will also need to pull the photometry test data::
  git clone https://github.com/Roman-Supernova-PIT/photometry_test_data.git

.. _running-snpit-container:

Running the SNPIT container
---------------------------

To use phrosty inside the container, you will need to run it with ``docker`` or ``podman``, and bind-mount the directory where you've cloned phrosty.  Phrosty requires a handful of additional directories:
* ``lc_out_dir`` : a place to write output lightcurves
* ``dia_out_dir`` : a place to write output difference images
* ``phrosty_temp`` : a place to write temp files; you want this on a fast filesystem
* ``phrosty_intermediate`` : a place to write intermediate data products for diagnostic purposes; you want this on a fast filesystem

You configure these directories with the phrosty config ``.yaml`` file.  For the config file we use for tests, inside the container these directories must show up at ``/lc_out_dir``, ``/dia_out_dir``, and ``/phrosty_temp``.  (The test environment unifies ``phrosty_intermediate`` and ``phrosty_temp``.)  You can make all of these diretories as subdirectories of your current directory::
  mkdir lc_out_dir
  mkdir dia_out_dir
  mkdir phrosty_temp

If you put them somewhere else, then make sure to modify the docker command below appropriately.
   
Assuming you're currently in the directory which is the parent of your ``phrosty`` and ``photometry_test_data`` checkouts, you can run a docker container suitable for running tests by running the following::
  docker run --gpus=all -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    --mount type=bind,source=$PWD/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/lc_out_dir,target=/lc_out_dir \
    --env PYTHONPATH=/roman_imsim \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    rknop/roman-snpit-env:cuda-dev \
    /bin/bash

(Substitute ``registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev`` for ``rknop/roman-snpit-env:cuda-dev`` if you pulled the docker image from there.)

If all is well, this will put you in a docker container.  You can tell you're in the container because your prompt will change to something like ``root@47394bd41fbe:/#`` (where the string of hexidecimal numbers will be different every time you start a container).  Verify that you've got access to the GPUs by running, inside the container::
  nvidia-smi

If you get an error message, or don't see at least one NVIDIA GPU listed, then phrosty will not work.

On NERSC Perlmutter
^^^^^^^^^^^^^^^^^^^

The procedure above is mostly right.  However, we **strongly** recommend you put your output directories on  on the Perlmutter scratch disk, at least for testing and development::
  mkdir $SCRATCH/phrosty_lc_out_dir
  mkdir $SCRATCH/phrosty_dia_out_dir
  mkdir $SCRATCH/phrosty_temp

(In fact, it's probably a good idea to put the other directories on ``$SCRATCH`` as well.)
  
Then, assuming you're in the directory above your ``phrosty`` and ``photometry_test_data`` checkouts, and assuming you've made the other two necessary directories, you can run the container with::
  podman-hpc run --gpu -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    --mount type=bind,source=$SCRATCH/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=$SCRATCH/phrosty_dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$SCRATCH/phrosty_lc_out_dir,target=/lc_out_dir \
    --env PYTHONPATH=/roman_imsim \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    --annotation run.oci.keep_original_groups=1 \
    registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev \
    /bin/bash  

If you're inside the container, your prompt will be something like ``root@f24c2ad04d6d:/#`` (though with a different string of hexidecimal digits (hexits?)).  If you do ``ls -F /``, you will see the various specific directories you mounted, such as ``/dia_out_dir`` and ``/photometry_test_data``.

Verify that you have access to GPUs by running::
  nvidia-smi

.. _Github repo: https://github.com/Roman-Supernova-PIT/phrosty
.. _tarball: https://github.com/Roman-Supernova-PIT/phrosty/tarball/master
