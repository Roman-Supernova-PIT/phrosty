*****
Usage
*****

Phrosty may be run from the command line by running ``python phrosty/pipeline.py`` (assuming you are in the top level of a github checkout).  If you're in :ref:`the necessary environment to run phrosty<phrosty-installation-prerequisites>`, then try running::

  python phrosty/pipeline.py --help

.. _example-usage:
Example Usage
=============

In addition to the examples below, see :ref:`running-tests`.

Running on Perlmutter
---------------------

This example is primarily intended for members of the Roman SN PIT, as it will require having an account on the NERSC Perlmutter cluster, and will require reading files that may not be accessible to people who aren't in the right unix groups.

Setting up the environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

You must first install the :ref:`docker image with the Roman SNPIT environment<running-snpit-container>`.  On Perlmutter, that means using ``podman-hpc``.  Pull the image with::

  podman-hpc pull registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev

**Note:** do _not_ run ``podman-hpc image pull ...``.  That will superficially seem to work, but will skip a step that gets run when you just do ``podman-hpc pull ...``.

(If the pull fails with a permission error, try running::

  podman-hpc login registry.nersc.gov

with your usual NERSC username and password, and the redoing the pull.  Once you've logged in once on perlmutter, it will probably remember it for a long time, having written some credential file somewhere in your home directory, so this login step will not usually be necessary.  It may be that only Roman SN PIT members have access to the m4385 repository.  If you're in the PIT and this doesn't work, contact Rob Knop.  If you're not in the PIT and this doesn't work, try using the image ``docker.io/rknop/roman-snpit-env:cuda-dev``.)

Verify that you got it with the command `podman-hpc images`, which should show (at least)::
  REPOSITORY                                      TAG         IMAGE ID      CREATED            SIZE        R/O
  registry.nersc.gov/m4385/rknop/roman-snpit-env  cuda-dev    07db212294db  About an hour ago  10.7 GB     false
  registry.nersc.gov/m4385/rknop/roman-snpit-env  cuda-dev    07db212294db  About an hour ago  10.7 GB     true

(The Image ID may not match exactly, and of course the Created time will be larger.)  The ``R/O false`` version is the image you pulled directly, and will only be available on the login node where you ran the pull.  More important is the ``R/O true`` (the readonly, or "squashed", image).  This one should be available on all login nodes and all compute nodes.  (You can verify this by running ``podman-hpc images`` on another login node; you should see only the ``R/O true`` image there.)

If you've pulled images before, and you're now working on a new login node, you will only see the ``R/O=true`` image.  That's the only image you really need, so in that case there's no need to pull the image again.

**If you have trouble with podman**: Podman is, as of this writing, still a beta feature on Nersc, and has some warts.  Refer to `NERSC's documentation on podman-hpc <https://docs.nersc.gov/development/containers/podman-hpc/overview/>`_.  In particular, if you want to clean the slate and start over, try running::
  podman-hpc system reset
to delete all of your podman images and contexts.  Then try pulling the image again.

Pick a place to work
^^^^^^^^^^^^^^^^^^^^

Work in one of two places.  You make yourself a subdirectory underneath ``/pscratch/sd/<u>/<username>``, where ``<username>`` is your NERSC username and `<u>` is the first letter of your username.  (You can get to this directory with ``cd $SCRATCH``; this is your top-level scratch directory, and NERSC sets the ``SCRATCH`` environment variable to point to it.)  Alternatively, you can create yourself a subdirectory somewhere underneath ``/global/cfs/cdirs/m4385/users``.  This is the shared SNPIT space on the NERSC community file system, so if you're going to work there, be aware that you're using up our shared file allocation.  At the moment, that's not a worry.

I'm going to call the place you've picked to work your "parent" directory.

Get phrosty
^^^^^^^^^^^

In your parent directory, :ref:`clone the phrosty repository<install-from-sources>`.

Locate existing directories
^^^^^^^^^^^^^^^^^^^^^^^^^^^

phrosty currently reads data from the OpenUniverse sims.  On NERSC, you can find the necessary information at the following directories.  These directories will be bind-mounted to the locations in parentheses (see below re: bind mounting).

* ``/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data`` (``/ou2024``)
* ``/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/PQ+HDF5_ROMAN+LSST_LARGE`` (``/ou2024_snana``)
* ``/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/ROMAN+LSST_LARGE_SNIa-normal`` (``/ou2024_snana_lc_dir``)
* ``/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/sims_sed_library`` (``/ou2024_sims_sed_library``)

Create needed directories
^^^^^^^^^^^^^^^^^^^^^^^^^

You need to make the following directories.  (They don't have to have exactly these names.  However, for purposes of the example, create these directories with these names as subdirectories under your parent directory.)

* ``dia_out_dir``
* ``lc_out_dir``

In addition, create a directory ``phrosty_temp`` somewhere underneath ``$SCRATCH``, e.g.::
  mkdir $SCRATCH/phrosty_temp
This directory will be mounted to ``/phrosty_temp`` inside the container.  (The further examples below will assume that this is where you made it.)

Secure lists of images for your supernova
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pick a supernova to run on.  TODO: more information.

For this example, we're going to run on the object with id 20172782.  In the ``examples/perlmutter`` directory under your ``phrosty`` checkout), you can find three ``.csv`` files that have information about the template and/or science images we're going to use:
* ``20172782_instances_templates_1.csv`` — a single R-band template image
* ``20172782_instances_templates_10.csv`` — 10 R-band template images
* ``20172782_instances_science.csv`` — 54 science images
* ``20172782_instances_science_2.csv`` — 2 science images

(Template images where chosen based on their simulated date relative to when the simulated supernova was active.)

For this example, you don't have to do anything, you will just use the files that are there.  However, if you are pushing this further, you will need to know how to find files, and how to construct your own ``.csv`` files.

If you look at these ``.csv`` files, there are give pieces of information on each line:
* The filename of the OpenUniverse image, relative to ``/ou2024/RomanTDS/images`` inside the container (see below).  On Perlmutter outside the container, these are relative to ``/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data/RomanTDS/images``.
* The pointing of the image
* The SCA on which the supernova is present for this pointing
* The MJD of the pointing
* The band (filter) of the exposure

Running interactively
^^^^^^^^^^^^^^^^^^^^^

The easiest way to just run something is to do it on an interactive node on Perlmutter.  (See :ref:`below<perlmutter-running-slurm>` for running it with slurm.)

First, get yourself a session on an interactive GPU node with::
  salloc -t 04:00:00 -A m4385 --constraint=gpu -q interactive

after a minute or so, that should log you into one of the nodes with a session that will last 4 hours.  (This is overkill; if you know it won't be that long, shorten the time after the ``-t`` flag.)  You can verify that you're on a compute node by running ``nvidia-smi``; you should see four different GPUs listed each with either 40MB or 80GB of memory, but no GPU processes running.

cd into your "parent" directory (if you're not there already).

Look at the file ``phrosty/examples/perlmutter/interactive_podman.sh``.  (There's no need to edit it; this is so you can see what's going on.)  You'll see number of ``--mount`` parameters.  Each of these takes a directory on the host machine (the ``source``) and maps it to a directory inside the podman container (the ``target``); this is "bind mounting".  For example, you will see your phrosty checkout goes to ``/phrosty`` inside the container.  In addition, several environment variables are set, and an "annotation" that is needed for ``podman-hpc`` to be able to handle accessing directories that are group-readable, but not world-readable.

Now do::
  bash phrosty/examples/perlmutter/interactive_podman.sh

This will create a container from the ``roman-snpit-env`` image, and put in a bash shell inside the container.  This will put you inside the container.  Your prompt will change to something like ``root@56356f1a4b9b:/usr/src#`` (where the hex barf will be different every time).  At any time, run ``ls -F /``; if you see directories ``phrosty``, ``phrosty_temp``, ``dia_out_dir``, and the others that were mounted by ``interactive_podman.sh``, then you know you're working inside the container, rather than on the host machine.  Verify that the GPUs are visible inside the container with ``nvidia-smi``.

Go to the ``/home`` directory, which is where your parent directory should be mounted::
  cd /home

The main Python executable for running the pipeline is ``phrosty/phrosty/pipeline.py``.  Run::
  python phrosty/phrosty/pipeline.py --help

to see how it works, and to see what the various parameters you can specify are.

Run this on your example lightcurve with::
  python phrosty/phrosty/pipeline.py \
    -c phrosty/examples/perlmutter/phrosty_config.yaml \
    --oid 20172782 \
    -r 7.551093401915147 \
    -d -44.80718106491529 \
    -b R062 \
    -t phrosty/examples/perlmutter/20172782_instances_templates_1.csv \
    -s phrosty/examples/perlmutter/20172782_instances_science_2.csv \
    -p 3 \
    -w 3

(If you run with ``.csv`` files that have larger number of images, you probably want to pass a larger number to `-p`; this is a number of parallel CPU processes that will run at once, and is limited by how many CPUs and how much memory you have available.  The code will only run one GPU process at once.  You can also try increasing `-w`, but this is more limited by filesystem performance than the number of CPUs and the amount of memory you have available.  We've set these both to 3 right now because there are only 3 files being processed (one template and two science images).  Empirically, on Perlmutter nodes, you can go up to something like `-p 15`; while there are (many) more CPUs than that, memory is the limiting factor.  Also, empirically, on Perlmutter, you can go up to something like `-w 5` before you reach the point of diminishing returns.  This is more variable, because whereas you have the node's CPUs to yourself, you're sharing the filesystem with the rest of the users of the system.)

If all is well, you should see a final line that looks something like::
  [2025-01-07 18:30:05 - phrosty - INFO] Results saved to /lc_out_dir/data/20172782/20172782_R062_all.csv

Outside the container (i.e. on Perlmutter), you should be able to find the file ``data/20172782/20172782_R062_all.csv`` underneath the ``lc_out_dir`` subdirectory of your parent directory.  Congratulations, this has the lightcurve!  (TODO: document the columns of this ``.csv`` file, but you can approximately guess what they are based on the column headers.)

You will also find new files in the ``dia_out_dir`` subdirectory, including several large ``.fits`` files.


Running with the NSight Profiler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO



.. _perlmutter-running-slurm:
Running a SLURM batch job
^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Phrosty Functionality
=====================

<<ALSO DOCUMENT FUNCTIONALITY &>>

