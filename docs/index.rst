phrosty Documentation
---------------------


This is the documentation for phrosty.

**Version**: |version|


This package contains the Python software suite developed for use with the Roman Telescope project, as part of the Roman Supernova Project Implementation Team (PIT) project. 

Statement of Need
=================

The Nancy Grace Roman Space Telescope (_Roman_) is NASA’s next flagship mission, designed for near-infrared (NIR) observations (Spergel et al. 2013, Spergel et al. 2015, Akeson et al. 2019). There will be several core community surveys that focus on different scientific goals. The High Latitude Time Domain Survey (HLTDS) is one such core survey and is optimized for transient astronomy, which is the study of astronomical objects that change over time. One of the survey's key goals is to collect data for a cosmological analysis with Type Ia Supernovae (SNe Ia). The core component of this survey will last two years, and have a cumulative 3792 hours of exposure time. The survey will yield an average of 241 files from imaging mode observations per day, which totals to approximately 0.4 TB of data (Roman Observations Time Allocation Committee: Final Report and Recommendations). 

Cosmological studies with SNe Ia require extraction of SN Ia light curves from images. Difference imaging analysis (DIA) is a commonly-used technique in transient analysis (e.g., SNe Ia) that allows both transient detection and photometric measurements. Broadly, this method involves subtracting two images, one with a transient and a ``reference'' image without a transient, to isolate the object of interest. 

Current DIA and forced photometry survey pipelines are insufficient for the volume of data expected from the Roman HLTDS survey over a 24 hour period. It would take 40 parallelized hours to make difference images (Kessler et al. 2015), resulting in a large backlog of processing, and detection inefficiencies that limit scientific discovery. There are many intermediate steps that are associated with SN Ia light curve extraction in addition to the creation of difference images; this means that the end-to-end time required to process these data from raw images to light curves will be much more that 40 hours.

_Roman_ images will pose several analytical challenges: (1) the large quantity of data that will be returned on a regular basis, and (2) the spatially-dependent nature of the image properties as a function of location on the detector, which directly affects ease of (3) obtaining the required <1% flux precision on SN Ia measurements. `phrosty` is a DIA and forced photometry pipeline that addresses these challenges. We integrate the saccadic fast fourier transform (SFFT) method for difference imaging analysis (DIA; Hu et al. 2022, Hu & Wang 2024), which uses a flexible delta basis function that can accommodate Roman’s spatially-varying point spread function (PSF). `phrosty` uses both GPU and CPU computation; this architecture has been carefully chosen to address the challenge of processing the large quantity of survey data within a reasonable amount of time.

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
