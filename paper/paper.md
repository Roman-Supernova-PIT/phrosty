---
title: 'phrosty: A difference imaging pipeline for Roman'

tags:
  - Python
  - astronomy
  - Roman
  - difference imaging
  - space-based
  - supernova
  - cosmology
  - SN Ia
  - SNe Ia
  - GPU

authors:
  - name: Lauren N. Aldoroty
    orcid: 0000-0003-0183-451X
    corresponding: true
    affiliation: 1
  - name: Lei Hu
    orcid: 0000-0001-7201-1938
    affiliation: 2
  - name: Rob A. Knop
    orcid: 0000-0002-3803-1641
    affiliation: 3
  - name: Cole Meldorf
    orcid: 0000-0002-6920-1498
    affiliation: 6
  - name: Daniel Scolnic
    orcid: 0000-0002-4934-5849
    affiliation: 1
  - name: Shu Liu 
    affiliation: 4
  - name: W. Michael Wood-Vasey
    affiliation: 4
  - name: Marcus Manos
    orcid: 0009-0004-0015-2052
    affiliation: 5
  - name: Lucas Erlandson
    orcid: 0000-0003-4544-6148
    affiliation: 5
  - name: Rebekah Hounsell
    orcid: 0000-0002-0476-4206
    affiliation: "7, 8"
  - name: Ben Rose
    orcid: 0000-0002-1873-8973
    affiliation: 9
  - name: Masao Sako
    orcid: 0000-0003-2764-7093
    affiliation: 6
  - name: Michael Troxel
    orcid: 0000-0002-5622-5212
    affiliation: 1
    
affiliations:
 - name: Department of Physics, Duke University,  Durham, NC 27708, USA
   index: 1
 - name: Department of Physics, Carnegie Mellon University, Pittsburgh, PA 15213, USA
   index: 2
 - name: Lawrence Berkeley National Laboratory, 1 Cyclotron Road, MS 50B-4206, Berkeley, CA 94720, USA
   index: 3
 - name: University of Pittsburgh
   index: 4
 - name: NVIDIA, 2788 San Tomas Expressway, Santa Clara, CA 95051
   index: 5
 - name: Department of Physics and Astronomy, University of Pennsylvania, 209 South 33rd Street, Philadelphia, PA 19104, USA
   index: 6
 - name: University of Maryland, Baltimore County, Baltimore, MD 21250, USA
   index: 7
 - name: NASA Goddard Space Flight Center, Greenbelt, MD 20771, USA
   index: 8
 - name: Department of Physics and Astronomy, Baylor University, One Bear Place 97316, Waco, TX 76798-7316, USA
   index: 9
   
date: 13 August 2025
bibliography: paper.bib

---

# Summary
NASA’s Nancy Grace Roman Space Telescope (_Roman_) will provide an opportunity to study dark energy with unprecedented precision using several techniques, including measurements of Type Ia Supernovae (SNe Ia). Here, we present `phrosty` (PHotometry for ROman with SFFT for tYpe Ia supernovae): a difference imaging pipeline for measuring the brightness of transient point sources in the sky, primarily SNe Ia, using _Roman_ data. `phrosty` is written in Python. We implement a GPU-accelerated version of the Saccadic Fast Fourier Transform (SFFT) method for difference imaging. Our GPU-optimized pipeline takes approximately 2 seconds per science image to process raw images and output measured aperture and PSF photometry. 

# Statement of need

Difference imaging analysis (DIA) is a technique in transient analysis that allows both transient detection and photometric measurements. Broadly, this method involves subtracting two images, one with a transient and a ``reference'' image without a transient, to isolate the object of interest. 

The Nancy Grace Roman Space Telescope (_Roman_) is NASA’s first survey-oriented flagship mission. Its design is optimized for cosmological studies, which includes extraction and precision analysis of Type Ia Supernova (SN Ia) light curves. Analysis of these data pose several challenges: (1) the enormous quantity of data that will be returned on a regular basis, and (2) the spatially-dependent nature of the image properties as a function of location on the detector, which directly affects ease of (3) obtaining the required <1\% flux precision on SN Ia measurements.

_Roman_'s High Latitude Time Domain Survey (HLTDS) is specifically designed for cosmological studies with SNe Ia. The core component of this survey will last two years, and have a cumulative 158 days of exposure time. The survey will yield an average of 241 Level 2 (rate files) from imaging mode observations per day, which totals to approximately 0.4 TB of data [@rotacreport]. Other DIA and forced photometry survey pipelines are insufficient to process this quantity of data. For example, the Dark Energy Survey (DES) pipeline, which uses the HOTPANTS algorithm [@hotpants], takes 10 minutes to create a difference image from an image from a single detector, with 8.3 million pixels per detector. [@kessler2015]. Using this pipeline, it would take 80 hours to make difference images from all 241 _Roman_ HLTDS images obtained within one 24 hour period, which have nearly 17 million pixels per detector. In reality, the amount of time spent processing these data will be much more than 80 hours because (1) there are many steps associated with SN Ia light curve analysis in addition to creating the difference images, (2) data from other concurrent surveys will also be processed for transient studies, and (3) this amount of time also excludes re-processing data when pipelines or calibration information are updated.

`phrosty` is a DIA and forced photometry pipeline that addresses these challenges. We integrate the saccadic fast fourier transform (SFFT) method for difference imaging analysis (DIA, [@sfft, @sfftjwst]), which uses a flexible delta basis function that can accommodate Roman’s spatially-varying point spread function (PSF). `phrosty` uses both GPU and CPU computation; this architecture has been carefully chosen to address the challenge of processing the large quantity of survey data within a reasonable amount of time. 

Because Roman software is currently under-development, there is no similar publicly available all-in-one pipeline that takes _Roman_ images as input and outputs light curves. The architecture choices made in development have substantially accelerated processing time such that `phrosty` is among the fastest astronomical DIA pipelines available. 

# Pipeline overview

![Visualization of pipeline steps using a sample image from the [@openuniverse] simulations. Green outlines represent the science image, purple outlines represent the template image, and black outlines represent the difference image. The small inset green and purple squares in panel 1 represent 100 x 100 pixel cutouts.\label{fig:pipeline}(pipelinefigure.png)]

![An example light curve, generated using the output of `phrosty`. Circles indicate measured output from `phrosty`, and solid lines are the known simulated light curve from [@openuniverse]. \label{fig:lc}(20172782_lc.png)]

## Pipeline steps

Our pipeline steps are as follows, with "`[GPU]`" indicating a GPU-accelerated step (Figure \autoref{fig:pipeline}):

1.  Identify any number of appropriate science (contains SN Ia at a specified coordinate) and template images (does not contain SN Ia at the same specified coordinate) for a given SN. All input template images will be paired with all science images such that if there are $S$ science images and $T$ template images, the total number of difference images is $S \times T$. 
2. The pipeline is then called in the command line by the user. For example, 

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

The inputs are: a configuration file (`-c`), the ID number assigned to the object (`--oid`), the object's RA in degrees (`-r`), the object's declination in degrees (`-d`), the _Roman_ filter matching the data being processed (`-b`), a list of template images (`-t`), a list of science images (`-s`), the number of parallel CPU processes for everything except file writing (`-p`), and the number of parallel CPU process for file writing (`-w`). 

3. Subtract the sky background and generate detection masks [@sourceextractor] for one pair of science and template images. Retrieve appropriate PSFs. 
4. `[GPU]` Align the reference image, reference PSF, and reference detection mask to the science image. 
5. `[GPU]` Cross-convolve the template PSF with the science image, and the science PSF with the template. 
6. `[GPU]` Apply SFFT subtraction [@sfft] to the cross-convolved images to produce a difference image.
7. `[GPU]` Fit and apply the decorrelation kernel to the difference image. Also apply the decorrelation kernel to the science PSF and cross-convolved science image to obtain the correct PSF for fitting, and the correct image for deriving a zero point from field stars for photometry.
8. Fit PSF from previous step to field stars in cross-convolved science image, cross-match stars to truth catalog, and record median of magnitudes for stars with $19 < m_{truth} < 21.5$ (i.e., not saturated and not noisy, [@aldoroty2025])
9. Fit the same PSF to the decorrelated difference image at the SN coordinates. 
10. Steps 3--9 are repeated in serial for each science and template image pair. Light curve data are collated into a human- and machine-readable tabular *.csv file for subsequent analysis. A file is generated for each command line call. Thus, each SN has a separate file for each filter's data that is processed. These tables are sufficient to plot the measured light curve of a SN (Figure \autoref{fig:lc}). 

Steps 4--7 are part of the SFFT algorithm. The steps that are not GPU-accelerated, including file writing, are run in parallel processes on CPUs. The number of parallel processes is controlled by the user. 

# GPU acceleration and performance
![NVIDIA Nsight Systems profile for `phrosty`, run on OpenUniverse SN 20172782. This object has 55 simulated _Roman_ observations. Red blocks indicate external libraries that cannot be accelerated, including sky subtraction using Source Extractor ([@sourceextractor], ~15.6 s), retrieving PSFs from `galsim` ([@galsim], ~25.5 s), making stamps of each observation (~535 ms), and photometry (~4.3 s). Purple blocks are image alignment and convolution (~670 ms). Blue blocks are subtraction (~600 ms). Grey blocks are fitting the decorrelation kernel, applying the decorrelation kernel, and writing the FITS file, which together take ~125.5 ms. The total run time, from raw images to photometry, is ~2 minutes. This code profile was generated using the NVIDIA Curiosity Cluster, which features 10 NVIDIA DGX SuperPODs in the NVIDIA DGX Cloud. Although the run time on the NVIDIA Curiosity Cluster is 2 minutes, the run time on NERSC's Perlmutter cluster is [PERLMUTTER TIME] due to file I/O limitations. \label{fig:codeprofile}(codeprofile.png)]

`phrosty` was developed in part during a hackathon hosted by NASA and Open Hackathons. During this time, `phrosty` developers, the SFFT developer [@sfft, @sfftjwst], and NVIDIA engineers began porting more of the original SFFT package to GPU using CUDA [@cuda]. The final version of SFFT used by `phrosty` carries out all matrix operations on GPU, and holds all data in-memory to minimize time spent on file I/O processes. The NVIDIA Tools Extension SDK (NVTX) is incorporated throughout the pipeline in order to enable easy code profiling, which is a unique feature compared to similar software suites. 

`phrosty` currently runs primarily on NVIDIA A100 GPUs with 40 GB memory at the National Energy Research Scientific Computing Center's (NERSC) Perlmutter supercomputer, and has also successfully been run at the NVIDIA Curiosity Cluster. Each SN fully occupies its GPU's capacity. Thus, the number of SNe that can be run in parallel is equal to the number of GPUs available. However, due to `phrosty`'s unique GPU-optimized architecture, each image takes approximately 2 seconds to process, including initial compilation and overhead time (Figure \ref{fig:codeprofile}). We process the full $4088 \times 4088$ px image for each subtraction. It is possible to decrease the processed stamp size, in turn decreasing the required memory and computation time, down to a minimum size of $1000\times1000$ px and still obtain high-quality results. Thus, it is among the fastest astronomical difference imaging pipelines ever developed. 

# Future development

The version of `phrosty` presented in this work is a prototype. It will continue to be developed, and updated versions will be announced in subsequent publications. This includes making the pipeline fully compatible with _Roman_ Science Operations Center (SOC) products. Currently, compatibility with the [@openuniverse] simulations is integrated into the pipeline such that it is reliant on their metadata, and uses the FITS format versions of all files (where applicable). The final data products associated with _Roman_ will not be FITS format; the project phases out the FITS file format, and replaces it with the Advanced Scientific Data Format (ASDF, [@asdf]). Thus, later versions of `phrosty` will be fully compatible with ASDF files. 

The largest barrier to `phrosty`'s community accessibility is the computational resources it requires. Although it requires minimal computing time per SN, it uses nearly 40 GB of GPU memory for a single object because fast fourier transforms (FFTs) are inherently memory-heavy operations. Thus, it requires access to expensive computing resources, which are not always available to the entire astronomical community. Although this hardware becomes more readily available over time as technological advancements continue, and we do not expect `phrosty` to ever require more resources than it demands in its current status, one of our goals for future development is to reduce `phrosty`'s memory consumption as much as possible. 

# Acknowledgements

L. A. thanks Megan Sosey for the discussion about _Roman_ HLTDS data volume.

This work is supported by NASA under award number 80GSFC24M0006. Additionally, funding for the Roman Supernova Project Infrastructure Team has been provided by NASA under contract to 80NSSC24M0023. M. T. was funded by NASA under JPL Contract Task 70-711320, ``Maximizing Science Exploitation of Simulated Cosmological Survey Data Across Surveys''.

This research used resources of the National Energy Research Scientific Computing Center, which is supported by the Office of Science of the U.S. Department of Energy using award number HEP-ERCAP32751.

This work was completed in part at the NASA Open Hackathon, part of the Open Hackathons program. The authors would like to acknowledge OpenACC-Standard.org for their support. 

# References
