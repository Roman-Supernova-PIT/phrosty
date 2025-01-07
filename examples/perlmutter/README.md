# Running phrosty on Perlmutter

## Preparatory warning

As of this writing, the `podman-hpc` installed on Perlmutter is broken.  Run
```
podman-hpc infohpc | grep -i version
```
If you see version 1.1.3, then the broken version is there.  Until this is on version 1.1.4 or later, in place of all `podman-hpc` commands below, run `/global/common/shared/tig/podman-hpc/bin/podman-hpc`.  (Strictly speaking, this is only necessary for pulling and migrating images, so you will not need to edit the `.sh` scripts where you actually create and run containers.)


## Set up environment.

### podman image

This example will use `podman-hpc` to run phrosty on Perlmutter.  You need to pull the Docker image with:
```
podman-hpc pull registry.nersc.gov/m4385/rknop/roman-snpit-dev:v0.0.1
```

Verify that you got it with the command `podman-hpc images`, which should show (at least):
```
REPOSITORY                                      TAG         IMAGE ID      CREATED      SIZE        R/O
registry.nersc.gov/m4385/rknop/roman-snpit-dev  v0.0.1      22f0c5ffe9ba  2 weeks ago  13.3 GB     false
registry.nersc.gov/m4385/rknop/roman-snpit-dev  v0.0.1      22f0c5ffe9ba  2 weeks ago  13.3 GB     true
```

(The Image ID may not match exactly, and of course the Created time will probably be different.)  The `R/O false` version is the image you pulled directly, and will only be available on the login node where you ran the pull.  More important is the `R/O true` (the readonly, or "squashed", image).  This one should be available on all login nodes and all compute nodes.  (You can verify this by running `podman-hpc images` on another login node; you should see only the `R/O true` image there.)

### Pick a place to work

Work in one of two places.  You make yourself a subdirectory underneath `/pscratch/sd/<u>/<username>`, where `<username>` is your NERSC username and `<u>` is the first letter of your username.  (You can get to this directory with `cd $SCRATCH`; this is your top-level scratch directory, and NERSC sets the `SCRATCH` environment variable to point to it.)  Alternatively, you can create yourself a subdirectory somewhere underneath `/global/cfs/cdirs/m4385/users`.  This is the shared SNPIT space on the NERSC community file system, so if you're going to work there, be aware that you're using up our shared file allocation.  At the moment, that's not a worry.

I'm going to call the place you've picked to work your "parent" directory.

### Get sfft

In your parent directory, do
```
git clone https://github.com/Roman-Supernova-PIT/sfft.git
cd sfft
git checkout fixes_20241022
cd ..
```

(Note: you may prefer to clone `git@github.com:Roman-Supernova-PIT/sfft.git` instead of `https:...`.  If you use github with ssh keys, then you already probably know that you want to do this.  If not, don't worry about it.)

In the future, we will merge the changes we need back to the `master` branch of sfft, but for now, you need the `fixes_20241022` branch.

### Get phrosty

```
git clone https://github.com/Roman-Supernova-PIT/phrosty.git
cd phrosty
git checkout main
cd ..
```

(Again, you may prefer to clone `git@github.com:Roman-Supernova-PIT/phrosty.git`, but you will already know if you want to do that.)

(ASIDE: while under code review, you actually need to checkout `rob_examples`, but that will change to `main`, and then we should remove this aside.)


### Locate existing directories

*phrosty* currently reads data from the OpenUniverse sims.  On NERSC, you can find the necessary information at the following directories:

* `/dvs_ro/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data` ("sims dir")
* `/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/PQ+HDF5_ROMAN+LSST_LARGE` ("snana_pq_dir")
* `/dvs_ro/cfs/cdirs/lsst/www/DESC_TD_PUBLIC/Roman+DESC/ROMAN+LSST_LARGE_SNIa-normal` ("snid_lc_dir")

### Create needed directories

You need to make the following directories.  (They don't have to have exactly these names.  However, for purposes of the example, create these directories with these names as subdirectories under your parent directory.)

* `sn_info_dir`
* `dia_out_dir`
* `lc_out_dir`

In addition, create a directory `phrosty_temp` somewhere underneath `$SCRATCH`, e.g.:
```mkdir $SCRATCH/phrosty_temp```.  (The further examples below will assume that this is where you made it.)

### Populate your sn_info_dir

This is where you put information that *phrosty* needs to find information about the OpenUniverse sims, and about any supernova from that sim you want to run it on.  The first file you need is `tds.yaml`; copy that file from this directory (i.e. `phrosty/examples/perlmutter/tds.yaml`) into your `sn_info_dir`.  (This file is a modified version of the standard OpenUniverse sims file `/global/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data/RomanTDS`, fixing some paths for our own purposes.)

### Secure lists of images for your supernova.

Pick a supernova to run on.  TODO: more information.

For this example, we're going to run on the object with id 20172782.  In the example directory, you can find three `.csv` files that have information about the template and/or science images we're going to use:
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


## Running interactively

The easiest way to just run something is to do it on an interactive node on Perlmutter.  (See below for running it with slurm.)

First, get yourself a session on an interactive GPU node with:
```
salloc -t 04:00:00 -A m4385 --constraint=gpu -q interactive
```
after a minute or so, that should log you into one of the nodes with a session that will last 4 hours.  (This is overkill; if you know it won't be that long, shorten the time after the `-t` flag.)  You can verify that you're on a compute node by running `nvidia-smi`; you should see four different GPUs listed each with either 40MB or 80GB of memory, but no GPU processes running.

Go to your "parent" directory

Copy the file `interactive_podman.sh` from this directory to your parent directory.  Look at this file.  You'll see number of `--mount` parameters.  Each of these takes a directory on the host machine (the `source`) and maps it to a directory inside the podman container (the `target`).  For example, you will see your phrosty checkout goes to `/phrosty` inside the container.  In addition, a bunch of environment variables are set so that *phrosty* will be able to find all of these directories inside the container.

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
  --oid 20172782 \
  -r 7.551093401915147 \
  -d -44.80718106491529 \
  -b R062 \
  -t 20172782_instances_templates_1.csv \
  -s 20172782_instances_science_two_images.csv \
  -p 4 \
  -w 4
```

(If you run with `.csv` files that have larger number of images, you probably want to pass a larger number to `-p`; this is a number of parallel CPU processes that will run at once, and is limited by how many CPUs and how much memory you have available.  The code will only run one GPU process at once.  You can also try increasing `-w`, but this is more limited by filesystem performance than the number of CPUs and the amount of memory you have available.)

If all is well, you should see a final line that looks something like:
```
[2025-01-07 18:30:05 - phrosty - INFO] Results saved to /lc_out_dir/data/20172782/20172782_R062_all.csv
```

Outside the container (i.e. on Perlmutter), you should be able to find the file `data/20172782/20172782_R062_all.csv` underneath the `lc_out_dir` subdirectory of your parent directory.  Congratulations, this has the lightcurve!  (TODO: document the columns of this `.csv` file, but you can approximately guess what they are based on the column headers.)


### Running with the NSight Profiler

When developing/debugging the pipeline, it's useful to run with a profiler, so you can see where the code is spending most of its time.  The huge `roman-snpit-dev` Docker image includes the NVIDIA NSight Systems profiler, and (at least as of this writing) the *phrosty* code includes hooks to flag parts of the code to the nsight profiler.  You can generate a profile for your code by running it with:

```
nsys profile \
  --trace-fork-before-exec=true \
  --python-backtrace=cuda \
  --python-sampling=true \
  --trace=cuda,nvtx,cublas,cusparse,cudnn,cudla,cusolver,opengl,openacc,openmp,osrt,mpi,nvvideo,vulkan,python-gil \
  python phrosty/phrosty/pipeline.py \
  --oid 20172782 \
  -r 7.551093401915147 \
  -d -44.80718106491529 \
  -b R062 \
  -t 20172782_instances_templates_1.csv \
  -s 20172782_instances_science_two_images.csv \
  -p 4 \
  -w 4
```

_Ideally_, this would create a file `report1.nsys-rep`, but something about that is broken; I'm not sure what.  It does create a file `report1.qdstrm`, which you can then somehow somewhere else convert to a `.nsys-rep` file.  On a Linux system, if you've installed the `night-compute` and `nsight-systems` packages (see https://docs.nvidia.com/nsight-systems/InstallationGuide/index.html), you can download the `.qdstrm` file to your system and run
```
/opt/nvidia/nsight-systems/2024.4.2/host-linux-x64/QdstrmImporter -i <name>.qdstrm
```
where `<name>.qstrm` is the file you downloaded.  (Note that the directory may have something other than `2024.4.2` in it, depending on what version you've installed.  For best comptibility with the version of Nsight in the `v0.0.1` docker image, I recommend installing `nsight-compute-2024.3.1` and  `nsight-systems-2024.4.2`.)  This should produce a file `<name>.nsys-rep`.  Then, on your local desktop, run
```
nsys-ui <name>.nsys-rep
```

to look at the profile.


## Running on a slurm queue

TODO
