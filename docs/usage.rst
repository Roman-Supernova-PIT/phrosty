*****
Usage
*****

Phrosty may be run from the command line by running ``python phrosty/pipeline.py`` (assuming you are in the top level of a github checkout).  If you're in :ref:`the necessary environment to run phrosty<phrosty-installation-prerequisites>`, then try running::

  python phrosty/pipeline.py -c phrosty/tests/phrosty_test_config.yaml --help

phrosty's behavior, and where it looks to find various images and other files it needs, are defined by a yaml config file.  You can find two examples of these files in:
* ``examples/perlmutter/phrosty_config.yaml``
* ``phrosty/tests/phrosty_test_config.yaml``
  
.. _example-usage:

Example Usage
=============

In addition to the examples below, see :ref:`running-tests`.

.. _manual-test-lightcurve:

Manually running a test lightcurve
----------------------------------

In this example, you will use data packaged with the photometry test archive to build a two-point lightcurve.

First, make sure your system meets the :ref:`system-requirements` and that you've downloaded the roman-snpit docker image as described in the :ref:`phrosty installation preqreuisties<phrosty-installation-prerequisites>`.

Next, make sure you've pulled down the ``phrosty`` archive as described in :ref:`installing phrosty from sources<install-from-sources>`.  Make sure also to install the photomery test data, as described there.

Finally, follow the instructions under :ref:`running-snpit-container`.

If all has gone well, you are now sitting inside a container that's ready to run phrosty.  Verify that you're in the container with ``ls -F /``.  Make sure that you see ``dia_out_dir/``, ``lc_out_dir/``, ``photometry_test_data/``, and ``phrosty_temp/`` in the list of directories.  Next, run ``nvidia-smi``, and make sure it shows you a GPU with 40MB.  Part of that output will look something like this::

  \|=========================================+========================+======================|
  |   0  NVIDIA A100-SXM4-40GB          Off |   00000000:03:00.0 Off |                    0 |
  | N/A   30C    P0             52W /  400W |       1MiB /  40960MiB |      0%      Default |
  |                                         |                        |             Disabled |

This shows a NVIDIA A100 GPU with 40GB of memory.  A different system might show::

  \|=========================================+========================+======================|
  |   0  NVIDIA GeForce RTX 3080 Ti     Off |   00000000:09:00.0  On |                  N/A |
  |  0%   60C    P8             35W /  350W |     774MiB /  12288MiB |      1%      Default |

This system has a NVIDIA RTX 3080 Ti consumer graphics card.  Notice that it only has 12288MiB of memory; 12GB of memory is not enough to run the example.

Inside the container, run::
  cd /home/phrosty
  pip install -e .

That will install the checked out version of phrosty in your currently running environment.  Note that if you exit the container and run a new container, you will have to ``pip install -e`` phrosty again, as the changes you make to a container only persist as long as that same container is still running.

Next, try running:
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

On your host system (as well as inside the container), you should see new files in wherever you put ``lc_out_dir``, ``dia_out_dir``, and ``phrosty_temp``.  (Inside the container, these are at ``/lc_out_dir``, ``/dia_out_dir``, and ``/phrosty_temp``.)


.. _perlmutter-example:

Running on Perlmutter
---------------------

While the previous example should have worked on Perlmutter, this is a somewhat more realistic example.  It doesn't use the photometry test data, but rather points to the full set of OpenUniverse2024 data available on Perlmutter.  This example is primarily intended for members of the Roman SN PIT, as it will require having an account on the NERSC Perlmutter cluster, and will require reading files that may not be accessible to people who aren't in the right unix groups.

This example will *probably* not work on a login node.  It might.  However, there is only a single GPU on each login node, and because of how some very annoying python libraries are written (jax, I'm looking at you), often people have allocated a large fraction of the available GPU memory even if they don't really need it.  You will probably need to run on a compute node.  This is described below.

Setting up the environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

Get your environment set up as described under the :ref:`phrosty installation prerequisites<phrosty-installation-prerequisites>`.

Pick a place to work
^^^^^^^^^^^^^^^^^^^^

Work in one of two places.  You make yourself a subdirectory underneath ``/pscratch/sd/<u>/<username>``, where ``<username>`` is your NERSC username and `<u>` is the first letter of your username.  (You can get to this directory with ``cd $SCRATCH``; this is your top-level scratch directory, and NERSC sets the ``SCRATCH`` environment variable to point to it.)  Alternatively, you can create yourself a subdirectory somewhere underneath ``/global/cfs/cdirs/m4385/users``.  This is the shared SNPIT space on the NERSC community file system, so if you're going to work there, be aware that you're using up our shared file allocation.  At the moment, that's not a worry.

I'm going to call the place you've picked to work your "parent" directory.

Get phrosty
^^^^^^^^^^^

In your parent directory, :ref:`clone the phrosty repository<install-from-sources>`.  For this example, you do not need to install the photometry test data.

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
* ``20172782_instances_science.csv`` — 53 science images
* ``20172782_instances_science_2.csv`` — 2 science images

(Template images where chosen based on their simulated date relative to when the simulated supernova was active.)

For this example, you don't have to do anything, you will just use the files that are there.  However, if you are pushing this further, you will need to know how to find files, and how to construct your own ``.csv`` files.

If you look at these ``.csv`` files, there are give pieces of information on each line:
* The filename of the OpenUniverse image, relative to ``/ou2024/RomanTDS/images`` inside the container (see below).  On Perlmutter outside the container, these are relative to ``/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data/RomanTDS/images``.
* The pointing of the image
* The SCA on which the supernova is present for this pointing
* The MJD of the pointing
* The band (filter) of the exposure

.. _perlmutter-interactive:

Running interactively
^^^^^^^^^^^^^^^^^^^^^

The easiest way to just run something is to do it on an interactive node on Perlmutter.  (See :ref:`below<perlmutter-running-slurm>` for running it with slurm.)

First, get yourself a session on an interactive GPU node with::
  salloc -t 04:00:00 -A m4385 --constraint=gpu -q interactive

after a minute or so, that should log you into one of the nodes with a session that will last 4 hours.  (This is overkill; if you know it won't be that long, shorten the time after the ``-t`` flag.)  You can verify that you're on a compute node by running ``nvidia-smi``; you should see four different GPUs listed each with either 40MB or 80GB of memory, but no GPU processes running.

cd into your "parent" directory (if you're not there already).

Look at the file ``phrosty/examples/perlmutter/interactive_podman.sh``.  (There's no need to edit it; this is so you can see what's going on.  If you read all of the :ref:`installation instructions<phrosty-installation>`, you will recognize a lot of what's there.)  You'll see number of ``--mount`` parameters.  Each of these takes a directory on the host machine (the ``source``) and maps it to a directory inside the podman container (the ``target``); this is "bind mounting".  For example, you will see your phrosty checkout goes to ``/phrosty`` inside the container.  In addition, several environment variables are set, and an "annotation" that is needed for ``podman-hpc`` to be able to handle accessing directories that are group-readable, but not world-readable.

Now do::
  bash phrosty/examples/perlmutter/interactive_podman.sh

This will create a container from the ``roman-snpit-env`` image, and put in a bash shell inside the container.  This will put you inside the container.  Your prompt will change to something like ``root@56356f1a4b9b:/usr/src#`` (where the hex barf will be different every time).  At any time, run ``ls -F /``; if you see directories ``phrosty``, ``phrosty_temp``, ``dia_out_dir``, and the others that were mounted by ``interactive_podman.sh``, then you know you're working inside the container, rather than on the host machine.  Verify that the GPUs are visible inside the container with ``nvidia-smi``.

Go to the ``/home`` directory, which is where your parent directory should be mounted::
  cd /home

The main Python executable for running the pipeline is ``phrosty/phrosty/pipeline.py``.  Run::
  SNPIT_CONFIG=phrosty/examples/perlmutter/phrosty_config.yaml python phrosty/phrosty/pipeline.py --help

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

**WARNING**: this section has not been tested recently so may be out of date.  TODO: try this again and update the docs after so doing.

When developing/debugging the pipeline, it's useful to run with a profiler, so you can see where the code is spending most of its time.  The huge ``roman-snpit-env:cuda-dev`` Docker image includes the NVIDIA NSight Systems profiler, and (at least as of this writing) the *phrosty* code includes hooks to flag parts of the code to the nsight profiler.  You can generate a profile for your code by doing everything described in :ref:`perlmutter-interactive` above, only replacing the final ``python`` command with::

  nsys profile \
    --trace-fork-before-exec=true \
    --python-backtrace=cuda \
    --python-sampling=true \
    --trace=cuda,nvtx,cublas,cusparse,cudnn,cudla,cusolver,opengl,openacc,openmp,osrt,mpi,nvvideo,vulkan,python-gil \
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

*Ideally*, this would create a file ``report1.nsys-rep``, but something about that is broken; I'm not sure what.  It does create a file ``report1.qdstrm``, which you can then somehow somewhere else convert to a ``.nsys-rep`` file.  On a Linux system, if you've installed the ``nsight-compute`` and ``nsight-systems`` packages (see `Nvidia's Nsight Systems installation guide <https://docs.nvidia.com/nsight-systems/InstallationGuide/index.html)>`_), you can download the ``.qdstrm`` file to your system and run::

  /opt/nvidia/nsight-systems/2024.4.2/host-linux-x64/QdstrmImporter -i <name>.qdstrm

where ``<name>.qstrm`` is the file you downloaded.  (Note that the directory may have something other than ``2024.4.2`` in it, depending on what version you've installed.  For best comptibility with the version of Nsight in the current (as of this writing) snpit docker image, I recommend trying to install something close to ``nsight-compute-2024.3.1`` and  ``nsight-systems-2024.4.2``; exactly what is avialable seems to vary with time.)  This should produce a file ``<name>.nsys-rep``.  Then, on your local desktop, run::
  nsys-ui <name>.nsys-rep

to look at the profile.


.. _perlmutter-running-slurm:

Running a SLURM batch job
^^^^^^^^^^^^^^^^^^^^^^^^^

For reference, see `the NERSC documentation on running jobs on Perlmutter <https://docs.nersc.gov/systems/perlmutter/running-jobs/>`_.  You need to set up your environment and run all of the steps above *before* "Running interactively".

**Create a job script**: to submit a job to a batch queue, you need to write a slrum script, which is just a shell script with some directives in the comments at the top.  An example script may be found in the file ``examples/perlmutter/20172782_slurm_demo.sh`` in your phrosty checkout.  If you look at this script, you will see that it contains mostly a combination of the ``podman-hpc`` and ``python phrosty/phrosty/pipeline.py`` commands above under "running interactively".  Instead of starting a shell with ``/bin/bash``, the ``podman-hpc`` command just runs the job directly.  It also adds a ``-w /home`` flag so it will be working in the right location.

At the top are the directives that control how the job is submitted.  Many of these you can leave as is.  (If you're morbidly curious, see `full documentation on the sbatch command <https://slurm.schedmd.com/sbatch.html>`_.  The ones you are most likely to want to change are

* ``#SBATCH --output <filename>`` : this is the filename that will hold all of the output written to the console for your job.  It will be written in the directory where you run ``slurm``.
* ``#SBATCH --qos shared`` : this tells slurm which queue to submit to.  See `NERSC's information on Perlmutter queues <https://docs.nersc.gov/jobs/policy/>`_.  By default, you want to submit to the ``shared`` queue.  Phrosty only currently uses a single GPU.  Each Perlmutter node has 4 GPUs, so if you submit to a queue that gives you an entire node, you're wasting it.  The shared queue has the advantage that *usually* jobs will start faster than they will on node-exclusive queues.  (You can sometimes wait days for a job on the regular queue to start!)  Additionally, our NERSC allocation will only be charged for the fraction of the node that we used.  However, when you're first testing, and you're only running a very small number of images, you might want to submit to the ``debug`` queue.  That allocates an entire node for the job, but _might_ start faster than jobs on the shared queue start.  (Try the shared queue first, though, because the job may well start within a few minutes.)
* ``#SBATCH --time 00:20:00`` : This is how long the job will run before the queue manager kills it.  The default, 20 minutes, *should* be plenty of time for the sample job that has 1 template image and 53 science images.  (However, if the NERSC filesystems are behaving poorly, it may not be enough time.)  If you're running a bigger job, then you need to specify more time.  Right now, assume something like ~1-2 minutes per image (which you can divide by the number of processes you run with ``-p``; see below), plus a few minutes of overhead.  Because phrosty is under heavy development and things are changing, this number will be highly variable.

You can probably leave the rest of the flags as is.  The ``--cpus-per-task`` and ``--gpus-per-task`` flags are set so that it will only ask for a quarter of a node.  (The queue manager is very particular about numbers passed to GPU nodes on the shared queue.  It needs you to ask for exactly 32 CPU cores for each GPU, and it needs you to ask for _exactly_ the right amount of memory.  The extra comment marks on the ``####SBATCH --mem`` line tell slurm to ignore it, as it seems to get the default right, and it's not worth fiddling with it to figure out what you should ask for.  A simple calculation would suggest that 64GB per GPU is what you should ask for, but when you do that, slurm thinks you're asking for 36 CPUs worth of memory, not 32 CPUs worth of memory.  The actual number is something like 56.12GB, but again, since the default seems to do the right thing, it's not worth fiddling with this.)

If look look at the bottom of the script, you will see that the number of parallel worker jobs that phrosty uses is set to 15 (``-p 15`` as a flag to ``python phrosty/phrosty/pipeline.py``).  The total number of processes that the python program runs at once is this, plus the number of FITS writer threads (given by ``-w``), plus one for the master process that launches all of the others.   You will notice that this total is less than the 32 CPUs that we nominally have.  To be safe, assume that each of the ``-p`` processes will use ~6GB of memory.  By limiting ourselves to 9 processes, we should safely fit within the amount of CPU memory allocated to the job (allowing for some overhead for the driver process and the FITS writer processes).  (TODO: we really want to get this memory usage down.)   Based on performance, you might want to play with the number of FITS writing threads (the number after ``-w``); assume that each FITS writer process will use ~1GB of memory.  (TODO: investigate how much they really use.)

**Make sure expected directories exists**: If you look at the batch script, you'll see a number of ``--mount`` flags that bind-mount directories inside the container.  From the location where you submit your job, all of the ``source=`` part of those ``--mount`` directives must be available.  For the demo, you will need to create the following directories underneath where you plan to submit the script::
  mkdir lc_out_dir
  mkdir dia_out_dir
  mkdir $SCRATCH/phrosty_temp

**Submitting your job**: Once you've are satisfied with your job script, submit it with::
  sbatch phrosty/examples/perlmutter/20172782_slurm_demo.sh

(Assuming you are running from the parent directory of your phrosty checkout.)

(replacing the argument with the actual name of your script).  (This example assumes that your current working directory is the parent directory of your phrosty checkout.)  If all is well, you should see an output something like::
  Submitted batch job 35680404

That number is the job id of your job.  If you see other things, they are probably error messages, and you need to fix what went wrong.

**Monitoring your job**: If you run::
  squeue --me

you will see all jobs you have submitted that are either pending or still running.  In the ``ST`` (short for "state") column, if you see ``PD``, it means your job is still pending, and hasn't started yet.  If you see ``R``, it means your job is running; in this case, the ``TIME`` column will tell you how long your job has been running.  If you see ``CG``, it means your job has recently finished (either succesfully, or with an error), and the system is currently cleaning it up.  If you see nothing, it means either that your job failed to submit (in which case you should have gotten an error message after your ``sbatch`` command above), or that it's finished.  ("Finished" may mean "exited right away with an error".)  Look at your output file to see what happened.

While the job is running, you can look at the output file to see how far it's gone and how it's doing.  (This is the file you specified on the ``#SBATCH --output`` line of your slurm script)

If you want to see the status of jobs that have completed, there are a few jobs you can run; try each of::
  scontrol show job <jobid>
  sacct -j <jobid>
  sacct -j <jobid> -o jobid,jobname,maxvmsize,reqmem,cputime --units=G
  seff <jobid>

(For more things you can pass to ``sacct``, see `its documentation <https://slurm.schedmd.com/sacct.html>`_.)  For all of those, ``<jobid>`` is the ID of your job on the slurm system.  While the job is still running you can see that job id in the left column of the output of ``squeue --me``.  After your job is over, you can look at the output file.  Assuming you used the example slurm script from this directory, you should see the jobid near the top of the output file.

**Checking job results:** Look at your output file.  The last line should be something like::
  [2025-02-10 15:43:32 - phrosty - INFO] Results saved to /lc_out_dir/data/20172782/20172782_R062_all.csv

and, ideally, there should be no lines anywhere in the file with ``ERROR`` near the beginning of the log message.

Note that ``/lc_out_dir/...`` is the absolute path _inside_ the container; it maps to ``lc_out_dir/...`` underneath your working directory where you ran ``sbatch``.  You will find the lightcurve in that ``.csv`` file.  There will also be a number of files written to the ``dia_out_dir`` directory.


Phrosty Functionality
=====================

<<ALSO DOCUMENT FUNCTIONALITY &>>

