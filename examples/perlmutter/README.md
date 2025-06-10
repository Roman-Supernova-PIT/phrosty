# Running phrosty on Perlmutter

Currently, this is primarily for internal development use of members of the Roman SN PIT.  You are welcome to try it if you aren't in the PIT, but at the moment we can't support you.  If you are in the PIT and have trouble with any of this, ping Rob Knop and Lauren Aldoroty in the `#photometry` channel on the Roman SNPIT Slack.

<a name="settingupenvironment"><a>

## Set up environment.

### podman image

This example will use `podman-hpc` to run phrosty on Perlmutter.  You need to pull the Docker image with:
```
podman-hpc pull registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev
```
Note: do _not_ run `podman-hpc image pull ...`.  That will superficially seem to work, but will skip a step that gets run when you just do `podman-hpc pull ...`.

(If the pull fails with a permission error, try running
```
podman-hpc login registry.nersc.gov
```
with your usual NERSC username and password, and the redoing the pull.  Once you've logged in once on perlmutter, it will probably remember it for a long time, having written some credential file somewhere in your home directory, so this login step will not usually be necessary.  It may be that only Roman SN PIT members have access to the m4385 repository.  If you're in the PIT and this doesn't work, contact Rob Knop.  If you're not in the PIT and this doesn't work, try using the image `docker.io/rknop/roman-snpit-env:cuda-dev`.)

Verify that you got it with the command `podman-hpc images`, which should show (at least):
```
REPOSITORY                                      TAG         IMAGE ID      CREATED            SIZE        R/O
registry.nersc.gov/m4385/rknop/roman-snpit-env  cuda-dev    07db212294db  About an hour ago  10.7 GB     false
registry.nersc.gov/m4385/rknop/roman-snpit-env  cuda-dev    07db212294db  About an hour ago  10.7 GB     true
```

(The Image ID may not match exactly, and of course the Created time will be larger.)  The `R/O false` version is the image you pulled directly, and will only be available on the login node where you ran the pull.  More important is the `R/O true` (the readonly, or "squashed", image).  This one should be available on all login nodes and all compute nodes.  (You can verify this by running `podman-hpc images` on another login node; you should see only the `R/O true` image there.)

If you've pulled images before, and you're now working on a new login node, you will only see the `R/O=true` image.  That's the only image you really need, so in that case there's no need to pull the image again.

### If you have trouble with podman

Podman is, as of this writing, still a beta feature on Nersc, and has some warts.  Refer to [NERSC's documentation on podman-hpc](https://docs.nersc.gov/development/containers/podman-hpc/overview/).  In particular, if you want to clean the slate and start over, try running
```
podman-hpc system reset
```
to delete all of your podman images and contexts.  Then try pulling the image again.


<a name="parentdirectory"></a>

### Pick a place to work

Work in one of two places.  You make yourself a subdirectory underneath `/pscratch/sd/<u>/<username>`, where `<username>` is your NERSC username and `<u>` is the first letter of your username.  (You can get to this directory with `cd $SCRATCH`; this is your top-level scratch directory, and NERSC sets the `SCRATCH` environment variable to point to it.)  Alternatively, you can create yourself a subdirectory somewhere underneath `/global/cfs/cdirs/m4385/users`.  This is the shared SNPIT space on the NERSC community file system, so if you're going to work there, be aware that you're using up our shared file allocation.  At the moment, that's not a worry.

I'm going to call the place you've picked to work your "parent" directory.

### Get phrosty

```
git clone https://github.com/Roman-Supernova-PIT/phrosty.git
cd phrosty
git checkout main
cd ..
```

This will pull the actual phrosty code, including the pipeline you'll run.  This README file you're reading is within that repository (in `examples/perlmutter/README.md`).  The `git checkout main` is probably redundant, becuase it's likely to check that branch out by default.  (Again, you may prefer to clone `git@github.com:Roman-Supernova-PIT/phrosty.git`, but you will already know if you want to do that.)  The last `cd ..` puts you back in your parent directory.


### Locate existing directories

*phrosty* currently reads data from the OpenUniverse sims.  On NERSC, you can find the necessary information at the following directories.  (The names in parentheses afterwards tell you which environment variables correspond to these directories.  You won't actually set these environment variables directly; that's handled in the shell scripts you'll run below.)

* `/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data` ("SIMS_DIR")
* `/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/PQ+HDF5_ROMAN+LSST_LARGE` ("SNANA_PQ_DIR")
* `/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/ROMAN+LSST_LARGE_SNIa-normal` ("SNID_LC_DIR")

### Create needed directories

You need to make the following directories.  (They don't have to have exactly these names.  However, for purposes of the example, create these directories with these names as subdirectories under your parent directory.)

* `sn_info_dir`
* `dia_out_dir`
* `lc_out_dir`

In addition, create a directories `phrosty_temp` and `phrosty_intermediate` somewhere underneath `$SCRATCH`, e.g.:
```mkdir $SCRATCH/phrosty_temp```
```mkdir $SCRATCH/phrosty_intermediate```

(The further examples below will assume that this is where you made it.)

### Populate your sn_info_dir

This is where you put information that *phrosty* needs to find information about the OpenUniverse sims, and about any supernova from that sim you want to run it on.  The first file you need is `tds.yaml`; copy that file from this directory (i.e. `phrosty/examples/perlmutter/tds.yaml`) into your `sn_info_dir`; when in the `sn_info_dir`, run:)
```
cp -p phrosty/examples/perlmutter/tds.yaml ./
```
(This file is a modified version of the standard OpenUniverse sims file `/global/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data/RomanTDS`, fixing some paths for our own purposes.)

<a name="securelistsofimagesforyoursupernova"></a>

### Secure lists of images for your supernova.

Pick a supernova to run on.  TODO: more information.

For this example, we're going to run on the object with id 20172782.  In the directory with this README file (i.e. `examples/perlmutter` under your `phrosty` checkout), you can find three `.csv` files that have information about the template and/or science images we're going to use:
* `20172782_instances_templates_1.csv` — a single R-band template image
* `20172782_instances_templates_10.csv` — 10 R-band template images
* `20172782_instances_science.csv` — 54 science images
* `20172782_instances_science_two_images.csv` — 2 science images

(Template images where chosen based on their simulated date relative to when the simulated supernova was active.)

Copy the `.csv` files you'll need to your parent directory.  For this example, we're going to use `...templates_1...` and `...science_two...`.  The other files are there in case you want to try a run with more data.

If you look at these `.csv` files, there are four pieces of information on each line:
* The filename of the OpenUniverse image, relative to `$SIMS_DIR`.  Inside the container (below), `$SIMS_DIR` will be at `/sims_dir`.  On Perlmutter outside the container, `$SIMS_DIR` is `/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data`.
* The pointing of the image
* The SCA on which the supernova is present for this pointing
* The MJD of the pointing
* The band (filter) of the exposure

## Running interactively

The easiest way to just run something is to do it on an interactive node on Perlmutter.  (See [below](#runningonaslurmqueue) for running it with slurm.)

First, get yourself a session on an interactive GPU node with:
```
salloc -t 04:00:00 -A m4385 --constraint=gpu -q interactive
```
after a minute or so, that should log you into one of the nodes with a session that will last 4 hours.  (This is overkill; if you know it won't be that long, shorten the time after the `-t` flag.)  You can verify that you're on a compute node by running `nvidia-smi`; you should see four different GPUs listed each with either 40MB or 80GB of memory, but no GPU processes running.

cd into your "parent" directory (if you're not there already).

Copy the files `interactive_podman.sh` and `phrosty_config.yaml` from `phrosty/examples/perlmutter` to your parent directory:
```
cp -p phrosty/examples/perlmutter/interactive_podman.sh phrosty_config.yaml ./
```
Look at the `.sh` file.  (There's no need to edit it; this is so you can see what's going on.)  You'll see number of `--mount` parameters.  Each of these takes a directory on the host machine (the `source`) and maps it to a directory inside the podman container (the `target`).  For example, you will see your phrosty checkout goes to `/phrosty` inside the container.  In addition, a bunch of environment variables are set so that *phrosty* will be able to find all of these directories inside the container.

Now do
```
bash interactive_podman.sh
```

This will put you inside the container.  Your prompt will change to something like `root@56356f1a4b9b:/usr/src#` (where the hex barf will be different every time).  At any time, run `ls -F /`; if you see directories `phrosty`, `phrosty_temp`, `roman_imsim`, and the others that were mounted by `interactive_podman.sh`, then you know you're working inside the container, rather than on the host machine.  Verify that the GPUs are visible inside the container with `nvidia-smi`.

Go to the `/home` directory, which is where your work directory should be mounted:
```
cd /home
```

The main Python executable for running the pipeline is `phrosty/phrosty/pipeline.py`.  Run
```
python phrosty/phrosty/pipeline.py --help
```
to see how it works, and to see what the various parameters you can specify are.

Run this on your example lightcurve with:
```
python phrosty/phrosty/pipeline.py \
  -c phrosty_config.yaml \
  --oid 20172782 \
  -r 7.551093401915147 \
  -d -44.80718106491529 \
  -b R062 \
  -t 20172782_instances_templates_1.csv \
  -s 20172782_instances_science_two_images.csv \
  -p 3 \
  -w 3
```

(If you run with `.csv` files that have larger number of images, you probably want to pass a larger number to `-p`; this is a number of parallel CPU processes that will run at once, and is limited by how many CPUs and how much memory you have available.  The code will only run one GPU process at once.  You can also try increasing `-w`, but this is more limited by filesystem performance than the number of CPUs and the amount of memory you have available.  We've set these both to 3 right now because there are only 3 files being processed (one template and two science images).  Empirically, on Perlmutter nodes, you can go up to something like `-p 15`; while there are (many) more CPUs than that, memory is the limiting factor.  Also, empirically, on Perlmutter, you can go up to something like `-w 5` before you reach the point of diminishing returns.  This is more variable, because whereas you have the node's CPUs to yourself, you're sharing the filesystem with the rest of the users of the system.)

If all is well, you should see a final line that looks something like:
```
[2025-01-07 18:30:05 - phrosty - INFO] Results saved to /lc_out_dir/data/20172782/20172782_R062_all.csv
```

Outside the container (i.e. on Perlmutter), you should be able to find the file `data/20172782/20172782_R062_all.csv` underneath the `lc_out_dir` subdirectory of your parent directory.  Congratulations, this has the lightcurve!  (TODO: document the columns of this `.csv` file, but you can approximately guess what they are based on the column headers.)

You will also find new files in the `dia_out_dir` subdirectory, including several large `.fits` files.

### Running with the NSight Profiler

When developing/debugging the pipeline, it's useful to run with a profiler, so you can see where the code is spending most of its time.  The huge `roman-snpit-dev` Docker image includes the NVIDIA NSight Systems profiler, and (at least as of this writing) the *phrosty* code includes hooks to flag parts of the code to the nsight profiler.  You can generate a profile for your code by running it with:

```
nsys profile \
  --trace-fork-before-exec=true \
  --python-backtrace=cuda \
  --python-sampling=true \
  --trace=cuda,nvtx,cublas,cusparse,cudnn,cudla,cusolver,opengl,openacc,openmp,osrt,mpi,nvvideo,vulkan,python-gil \
  python phrosty/phrosty/pipeline.py \
  -c phrosty_config.yaml \
  --oid 20172782 \
  -r 7.551093401915147 \
  -d -44.80718106491529 \
  -b R062 \
  -t 20172782_instances_templates_1.csv \
  -s 20172782_instances_science_two_images.csv \
  -p 3 \
  -w 3
```

_Ideally_, this would create a file `report1.nsys-rep`, but something about that is broken; I'm not sure what.  It does create a file `report1.qdstrm`, which you can then somehow somewhere else convert to a `.nsys-rep` file.  On a Linux system, if you've installed the `nsight-compute` and `nsight-systems` packages (see [Nvidia's Nsight Systems installation guide](https://docs.nvidia.com/nsight-systems/InstallationGuide/index.html)), you can download the `.qdstrm` file to your system and run
```
/opt/nvidia/nsight-systems/2024.4.2/host-linux-x64/QdstrmImporter -i <name>.qdstrm
```
where `<name>.qstrm` is the file you downloaded.  (Note that the directory may have something other than `2024.4.2` in it, depending on what version you've installed.  For best comptibility with the version of Nsight in the `v0.0.1` docker image, I recommend trying to install something close to `nsight-compute-2024.3.1` and  `nsight-systems-2024.4.2`; exactly what is avialable seems to vary with time.)  This should produce a file `<name>.nsys-rep`.  Then, on your local desktop, run
```
nsys-ui <name>.nsys-rep
```

to look at the profile.


<a name="runningonaslurmqueue"></a>

## Running on a slurm queue

For reference, see [the NERSC documentation on running jobs on Perlmutter](https://docs.nersc.gov/systems/perlmutter/running-jobs/).

You need to [set up your enviroment](#settingupenvironment) the same way you would for an interactive job.  You want to be working in your [parent directory](#parentdirectory) (the directory where you cloned the `phrosty` and `sfft` archives).

### Creating a job script

To submit a job to a batch queue, you need to write a slrum script, which is just a shell script with some directives in the comments at the top.  An example script may be found in the file `20172782_2sci1templ.sh` in the same directory where this `README.md` file exists.  If you look at this script, you will see that it contains mostly a combination of the `podman-hpc` and `python phrosty/phrosty/pipeline.py` commands above under "running interactively".  Instead of starting a shell with `/bin/bash`, the `podman-hpc` command just runs the job directly.  It also adds a `-w /home` flag so it will be working in the right location.

Copy this script to the directory where you are working, and edit it if necessary.

You can adjust this job for whatever OpenUniverse candidate and images you want to run by editing the values of the `--oid`, `-r`, `-d`, `-b`, `-t`, and `-s` flags passed to `python phrosty/phrosty/pipeline.py`.  (The meanings of these flags are documented by `python phrosty/phrosty/pipeline.py --help`.  See [secure lists of images for your supernova](#securelistsofimagesforyoursupernova) above.)

At the top are the directives that control how the job is submitted.  Many of these you can leave as is.  (If you're morbidly curious, see [full documentation on the `sbatch` command](https://slurm.schedmd.com/sbatch.html).)  The ones you are most likely to want to change are

* `#SBATCH --output <filename>` : this is the filename that will hold all of the output written to the console for your job.  It will be written in the directory where you run `slurm`.
* `#SBATCH --qos shared` : this tells slurm which queue to submit to.  See [NERSC's information on Perlmutter queues](https://docs.nersc.gov/jobs/policy/).  By default, you want to submit to the `shared` queue.  Phrosty only currently uses a single GPU.  Each Perlmutter node has 4 GPUs, so if you submit to a queue that gives you an entire node, you're wasting it.  The shared queue has the advantage that _usually_ jobs will start faster than they will on node-exclusive queues.  (You can sometimes wait days for a job on the regular queue to start!)  Additionally, our NERSC allocation will only be charged for the fraction of the node that we used.  However, when you're first testing, and you're only running a very small number of images, you might want to submit to the `debug` queue.  That allocates an entire node for the job, but _might_ start faster than jobs on the shared queue start.  (Try the shared queue first, though, because the job may well start within a few minutes.)
* `#SBATCH --time 00:10:00` : This is how long the job will run before the queue manager kills it.  The default, 10 minutes, is more than enough for the sample script that has two science and one template image.  If you're running a bigger job, then you need to specify more time.  Right now, assume that you need a couple of minutes times the number of science images times the number of template images, plus a few minutes of overhead.  Because phrosty is under heavy development and things are changing, this number will be highly variable; it may well be that 2 minutes per science/template pair is overkill, but it may be that because of where phrosty is in development it takes a lot longer than that.

You can probably leave the rest of the flags as is.  The `--cpus-per-task` and `--gpus-per-task` flags are set so that it will only ask for a quarter of a node.  (The queue manager is very particular about numbers passed to GPU nodes on the shared queue.  It needs you to ask for exactly 32 CPU cores for each GPU, and it needs you to ask for _exactly_ the right amount of memory.  The extra comment marks on the `####SBATCH --mem` line tell slurm to ignore it, as it seems to get the default right, and it's not worth fiddling with it to figure out what you should ask for.  A simple calculation would suggest that 64GB per GPU is what you should ask for, but when you do that, slurm thinks you're asking for 36 CPUs worth of memory, not 32 CPUs worth of memory.  The actual number is something like 56.12GB, but again, since the default seems to do the right thing, it's not worth fiddling with this.)

If look look at the bottom of the script, you will see that the number of parallel worker jobs that phrosty uses is set to 15 (`-p 15` as a flag to `python phrosty/phrosty/pipeline.py`).  The total number of processes that the python program runs at once is this, plus the number of FITS writer threads (given by `-w`), plus one for the master process that launches all of the others.   You will notice that this total is less than the 32 CPUs that we nominally have.  To be safe, assume that each of the `-p` processes will use ~3GB of memory.  By limiting ourselves to 15 processes, we should safely fit within the amount of memory allocated to the job (allowing for some overhead for the driver process and the FITS writer processes).   Based on performance, you might want to play with the number of FITS writing threads (the number after `-w`); assume that each FITS writer process will use ~1GB of memory.  (TODO: investigate how much they really use.)

### Submitting your job

Once you've are satisfied with your job script, submit it with
```
sbatch 20172782_2csi1temp.sh
```
(replacing the argument with the actual name of your script).  If all is well, you should see an output something like:
```
Submitted batch job 35680404
```
That number is the job id of your job.  If you see other things, they are probably error messages, and you need to fix what went wrong.

### Monitoring your job

If you run
```
squeue --me
```
you will see all jobs you have submitted that are either pending or still running.  In the `ST` (short for "state") column, if you see `PD`, it means your job is still pending, and hasn't started yet.  If you see `R`, it means your job is running; in this case, the `TIME` column will tell you how long your job has been running.  If you see `CG`, it means your job has recently finished (either succesfully, or with an error), and the system is currently cleaning it up.  If you see nothing, it means either that your job failed to submit (in which case you should have gotten an error message after your `sbatch` command above), or that it's finished.  ("Finished" may mean "exited right away with an error".)  Look at your output file to see what happened.

While the job is running, you can look at the output file to see how far it's gone and how it's doing.  (This is the file you specified on the `#SBATCH --output` line of your slurm script)

If you want to see the status of jobs that have completed, there are a few jobs you can run; try each of:
```
scontrol show job <jobid>
sacct -j <jobid>
sacct -j <jobid> -o jobid,jobname,maxvmsize,reqmem,cputime --units=G
seff <jobid>
```
(For more things you can pass to `sacct`, see [its documentation](https://slurm.schedmd.com/sacct.html).)  For all of those, `<jobid>` is the ID of your job on the slurm system.  While the job is still running you can see that job id in the left column of the output of `squeue --me`.  After your job is over, you can look at the output file.  Assuming you used the example slurm script from this directory, you should see the jobid near the top of the output file.


### Checking job results

Look at your output file.  The last line should be something like
```
[2025-02-10 15:43:32 - phrosty - INFO] Results saved to /lc_out_dir/data/20172782/20172782_R062_all.csv
```
Note that `/lc_out_dir/...` is the absolute path _inside_ the container; it maps to `lc_out_dir/...` underneath your working directory where you ran `sbatch`.  You will find the lightcurve in that `.csv` file.  There will also be a number of files written to the `dia_out_dir` subdirectory.
