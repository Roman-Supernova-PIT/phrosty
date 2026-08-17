.. highlight:: shell

.. _phrosty-installation:

============
Installation
============

.. contents::

.. _system-requirements:

System Requirements
-------------------

``phrosty`` can run using either a ``cupy`` (CUDA 12.4, requires an NVIDIA GPU) or ``numpy`` backend (CPU). Empirically, you will need at least 36 GB GPU memory or 56 GB CPU memory to run these backends, respectively. 

To properly set up ``phrosty``, you need to follow **one** of the sections in `Environment set-up<phrosty-environment-setup>`, followed by `Install from sources<install-from-sources>`.

.. _phrosty-environment-setup:

Environment set-up
------------------

There have been a number of ways to run ``phrosty`` as the package has evolved. Here is my best attempt at preserving all of what is still relevant. If you aren't on the SNPIT, then you probably want "`I do not have a 40 GB-memory NVIDIA GPU<phrosty-local>`". If you are on the SNPIT, then you know which section you need.

.. _phrosty-local:

I do not have 40 GB-memory NVIDIA GPU
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You are most people. 

Do::

  pip install r

.. _phrosty-general-docker:

I have a 40 GB-memory NVIDIA GPU
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have a docker container with all the prerequisites set up for you. Pull the image by doing::
  
  docker pull docker.io/rknop/roman-snpit-env:cuda-dev

**Note that if you choose to run natively, i.e. without the Docker container, you will need to install ``cupy``**::

  pip install cupy-cuda12x

To use phrosty inside the container, you will need to run it with ``docker`` or ``podman``, and bind-mount the directory where you've cloned phrosty.  Phrosty requires a handful of additional directories:

* ``lc_out_dir`` : a place to write output lightcurves
* ``dia_out_dir`` : a place to write output difference images
* ``phrosty_temp`` : a place to write temporary files; you want this on a fast filesystem
* ``scratch`` : another place to write temporary files. Yes, we need to consolidate these, but it is a low-priority task at this time.

You configure these directories with the phrosty config ``.yaml`` file.  For the config file we use for tests, inside the container these directories must show up at ``/lc_out_dir``, ``/dia_out_dir``, ``/phrosty_temp``, and ``/scratch``. You can make all of these diretories as subdirectories of your current directory::

  mkdir lc_out_dir
  mkdir dia_out_dir
  mkdir phrosty_temp
  mkdir scratch

If you put them somewhere else, then make sure to modify the docker command below appropriately.

Assuming you're currently in the directory which is the parent of your ``phrosty`` and ``photometry_test_data`` checkouts:

  docker run --gpus=all -it \
    --mount type=bind,source=$PWD,target=/home \
    --mount type=bind,source=$PSCRATCH,target=/scratch \
    --mount type=bind,source=$PWD/photometry_test_data,target=/photometry_test_data \
    --mount type=bind,source=$PWD/phrosty_temp,target=/phrosty_temp \
    --mount type=bind,source=$PWD/dia_out_dir,target=/dia_out_dir \
    --mount type=bind,source=$PWD/lc_out_dir,target=/lc_out_dir \
    --env LD_LIBRARY_PATH=/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:/usr/local/cuda/lib64/stubs \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env NUMEXPR_NUM_THREADS=1 \
    --env OMP_NUM_THREADS=1 \
    --env VECLIB_MAXIMUM_THREADS=1 \
    --env TERM=xterm \
    --annotation run.oci.keep_original_groups=1 \
    rknop/roman-snpit-env:cuda-dev-0.1.41 \
    /bin/bash

**You may need to modify these paths.** Note that 0.1.41 will increment over time.

If all is well, this will put you in a docker container.  You can tell you're in the container because your prompt will change to something like ``root@47394bd41fbe:/#`` (where the string of hexidecimal numbers will be different every time you start a container).  Verify that you've got access to the GPUs by running, inside the container::

  nvidia-smi

If you get an error message, or don't see at least one NVIDIA GPU listed, then this will not work.

.. _phrosty-smdc:

Installing on SMDC
^^^^^^^^^^^^^^^^^^

**This section will work for SN PIT members.** 

General instructions for accessing SMDC can be found `in the wiki <https://github.com/Roman-Supernova-PIT/Roman-Supernova-PIT/wiki/NASA-SMDC-%28AWS%29>`_.

There is some information at `"Working with PIT Images" here <https://github.com/Roman-Supernova-PIT/Roman-Supernova-PIT/wiki/SMCE-Containers>`_, which is largely applicable if you want to use the Singularity containers.

First, ``salloc`` a node. If you want a GPU node, do::

  salloc -p gpu-int --time=04:00:00

If you want to run on a CPU node, do::

  salloc -p mem-lg --time=04:00:00

Sometimes your correct set of groups won't be correctly populated on a compute node due to a race condition between populating the container and correctly configuring the active directory lookup. You will see a message about this when you get your node that says that the groups weren't loaded properly. Also, if you type `groups` on the login node, you'll see `[your username] spack cluster_users snpit`. If you type `groups` on the GPU node, you'll see `[your username] nogroup`. You will also hit a permissions issue running `phrosty` when it tries to write files outside your home directory. To start a new terminal that will have the groups loaded correctly, do::

  ssh localhost

I want to use the Singularity/Apptainer container
"""""""""""""""""""""""""""""""""""""""""""""""""
PENDING: This will link to snappl documentation when snappl PR #214 is merged. 


I want to use the shared virtual environment
""""""""""""""""""""""""""""""""""""""""""""

PENDING: This will link to snappl documentation also.

I want my own development virtual environment that I can change
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**UPDATE THIS WHEN environment/#23 PR IS MERGED.**

You need to install ``cupy``. Do::

  pip install cupy-cuda12x

.. _phrosty-nersc-perlmutter:

Installing on NERSC Perlmutter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**This section will work for SN PIT members, and maybe anyone else with access to NERSC Perlmutter, which ostensibly could be you.** 

If you are on NERSC Perlmutter, you have access to NVIDIA GPUs with 40 GB GPU memory. Run ``module list``.  Make sure that ``cudatoolkit/12.4`` shows up in your list of modules.  If not, you may need to adjust the modules you have loaded.

(Ideally, because phrosty runs inside a container, the specific version of CUDA on the host system wouldn't matter.  However, containers can be touchy about hooking up GPUs inside containers.  If you're having trouble seeing the GPU inside the container, check for mismatches between the versions of cuda inside and outside of the container.)

Currently, phrosty is designed inside a container built from the the `Roman Supernova PIT environment <https://github.com/Roman-Supernova-PIT/environment>`_.

Because phrosty (and other libraries it depends on, such as snappl) is under heavy development, it's possible that the latest container will not work properly with phrosty at any given moment.

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

Assuming you're in the directory above your ``phrosty`` and ``photometry_test_data`` checkouts, you can run the container with ``bash /global/cfs/cdirs/m4385/env/interactive-podman-nov2025.sh``. At this time, both of these files are the same, but you have the ability to modify the one in `examples/perlmutter` and not the one in `m4385/env`. 

If you absolutely must make your own container for some reason, see `the interactive podman scripts in our environment directory <https://github.com/Roman-Supernova-PIT/environment/blob/main/interactive-podman-nov2025.sh>`_ for reference.

If you're inside the container, your prompt will be something like ``root@f24c2ad04d6d:/#`` (though with a different string of hexidecimal digits (hexits?)).  If you do ``ls -F /``, you will see the various specific directories that are mounted in the above command.

Verify that you have access to GPUs by running::

  nvidia-smi

.. _install-from-sources:

Installing from sources
-----------------------

You will need the SNPIT's photometry package ``snappl``, as well as our version of SFFT, in order to run ``phrosty``. The latest stable versions are on ``pip``:

  pip install roman-snpit-snappl sfft-romansnpit

Currently, the only way to install ``phrosty`` is to download it from the `github repo <https://github.com/Roman-Supernova-PIT/phrosty>`_.  Clone it with::

    git clone https://github.com/Roman-Supernova-PIT/phrosty.git

(you can also clone it via the ``git@`` code link if you know what you're doing.)

Then, ``cd`` to the ``phrosty`` folder and ``pip install .``. Or ``pip install -e .`` if you expect to do development.

For SNPIT development only: If you need a more recent version of SFFT than what's in the docker image, use the Roman SNPIT SFFT fork::

    https://github.com/Roman-Supernova-PIT/sfft.git

Make sure you check these out to the same parent folder. 

Installing the photometry test data (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to run tests, and some of the examples, then you will also need to pull the photometry test data::

  git clone https://github.com/Roman-Supernova-PIT/photometry_test_data.git

.. _Github repo: https://github.com/Roman-Supernova-PIT/phrosty
.. _tarball: https://github.com/Roman-Supernova-PIT/phrosty/tarball/master

.. Pulling the container image on other HPC Systems
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. If your HPC system doesn't support containers, then you're out of luck.  Most HPC systems support containers using ``apptainer`` (previously known as ``singularity``).  This system is not a drop-in replacment for docker, as the way you obtain images, and the semantics for running containers, are different.  (There are also differences in terms of how isolated the environment is; singularity is less of a "real" isolated container than docker, usually.)

.. TODO : document use of singularity.