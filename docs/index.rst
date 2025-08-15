phrosty Documentation
---------------------


This is the documentation for phrosty.

**Version**: |version|


This package contains the Python software suite developed for use with the Roman Telescope project, as part of the Roman Supernova Project Implementation Team (PIT) project. 

Statement of Need
=================

Difference imaging analysis (DIA) is a technique in transient analysis that allows both transient detection and photometric measurements. Broadly, this method involves subtracting two images, one with a transient and a "reference" image without a transient, to isolate the object of interest. 

The Nancy Grace Roman Space Telescope (_Roman_) is NASA’s first survey-oriented flagship mission. Its design is optimized for cosmological studies, which includes extraction and precision analysis of Type Ia Supernova (SN Ia) light curves. Analysis of these data pose several challenges: (1) the enormous quantity of data that will be returned on a regular basis, and (2) the spatially-dependent nature of the image properties as a function of location on the detector, which directly affects ease of (3) obtaining the required &lt;1\% flux precision on SN Ia measurements.

_Roman_'s High Latitude Time Domain Survey (HLTDS) is specifically designed for cosmological studies with SNe Ia. The core component of this survey will last two years, and have a cumulative 158 days of exposure time. The survey will yield an average of 241 Level 2 (rate files) from imaging mode observations per day, which totals to approximately 0.4 TB of data. Other DIA and forced photometry survey pipelines are insufficient to process this quantity of data. For example, the Dark Energy Survey (DES) pipeline, which uses the HOTPANTS algorithm, takes 10 minutes to create a difference image from an image from a single detector, with 8.3 million pixels per detector. Using this pipeline, it would take 80 hours to make difference images from all 241 _Roman_ HLTDS images obtained within one 24 hour period, which have nearly 17 million pixels per detector. In reality, the amount of time spent processing these data will be much more than 80 hours because (1) there are many steps associated with SN Ia light curve analysis in addition to creating the difference images, (2) data from other concurrent surveys will also be processed for transient studies, and (3) this amount of time also excludes re-processing data when pipelines or calibration information are updated.

`phrosty` is a DIA and forced photometry pipeline that addresses these challenges. We integrate the saccadic fast fourier transform (SFFT) method for difference imaging analysis (DIA), which uses a flexible delta basis function that can accommodate Roman’s spatially-varying point spread function (PSF). `phrosty` uses both GPU and CPU computation; this architecture has been carefully chosen to address the challenge of processing the large quantity of survey data within a reasonable amount of time. 

Because Roman software is currently under-development, there is no similar publicly available all-in-one pipeline that takes _Roman_ images as input and outputs light curves. The architecture choices made in development have substantially accelerated processing time such that `phrosty` is among the fastest astronomical DIA pipelines available. 

Contact and Support
===================

Individuals who wish to contribute to phrosty, report issues, or seek support are encouraged to submit an issue via [GitHub](https://github.com/Roman-Supernova-PIT/phrosty/issues) and use the pre-loaded templates for feature requests and issue reports. Please include as much detail as possible, including a description of the problem, any associated error messages, inputs, and details about the environment the user is running. Please adhere to the phrosty [code of conduct](https://github.com/Roman-Supernova-PIT/phrosty/blob/main/CODE_OF_CONDUCT.md).

--------------------------------

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation.rst
   usage.rst
   development.rst
   changes.rst
   api.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
