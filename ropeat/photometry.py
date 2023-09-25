import numpy as np
from scipy import stats
from glob import glob
from copy import deepcopy
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import ZScaleInterval, simple_norm
from astropy.wcs import WCS
import astropy.units as u
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats
from photutils.background import LocalBackground, MMMBackground
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFBuilder, extract_stars, PSFPhotometry, IterativePSFPhotometry, IntegratedGaussianPRF

class truth():
    # https://roman.ipac.caltech.edu/data/sims/sn_image_sims/galaxia_akari.fits.gz
    def __init__(self,filepath):
        self.filepath = filepath
        self.truthhdu = fits.open(self.filepath)

    def table(self):
        """Returns an astropy table with truth data. 
        :param filepath: Path to truth file. File must be in fits format. 
        :type filepath: str
        :return: astropy table with columns ['index', 'ra', 'dec].
        :rtype: astropy table
        """

        truthhead = self.truthhdu[1].header
        truth = self.truthhdu[1].data
        truth_ra = (truth['ra']*u.radian).to(u.deg)
        truth_dec = (truth['dec']*u.radian).to(u.deg)
        truth_tab = Table([np.arange(0,len(truth),1),truth_ra,truth_dec], names=('index','ra','dec'))
        return truth_tab

    def truthcoords(self):
        """
        :return: astropy SkyCoords pairs from truth file with 
                    (ra, dec) in degrees. 
        :rtype: astropy SkyCoord object
        """        
        return SkyCoord(self.table()['ra'],
                        self.table()['dec'],
                        unit=(u.deg,u.deg))
    

class scienceimg():
    """
    Define class for science image. Initialize with:
    :param filepath: Path to science image. File must be in fits format.
    :type filepath: str
    :param truthcoords: Truth coordinates from truth file. Must be an astropy
                        SkyCoord object. Can input from ropeat.photometry.truth.truthcoords().
    :type truthcoords: astropy SkyCoord object
    :param bkgstat: Statistical method for calculating background.
    :type bkgstat: photutils.background object, optional
    :param band: WFI bandpass. Valid inputs: ['Y106', 'J129', 'H158', 'F184']
    :type band: str
    """
    def __init__(self,filepath,
                 truthcoords,
                 band='unspecified',
                 pointing='unspecified',
                 chip='unspecified',
                 bkgstat=MMMBackground(),
                 bkg_annulus=(50.0,80.0)
                 ):
        
        acceptable_bands = ['Y106', 'J129', 'H158', 'F184']
        if band not in acceptable_bands:
            raise ValueError(f'Argument "band" must be in {acceptable_bands}.')

        if pointing == 'unspecified':
            raise ValueError('You need to specify a pointing ID.')

        if chip == 'unspecified':
            raise ValueError('You need to specify a chip number.')

        self.filepath = filepath
        self.sciencehdu = fits.open(self.filepath)
        self.data = self.sciencehdu[1].data
        self.wcs = WCS(self.sciencehdu[1].header)

        self.truthcoords = truthcoords
        self.footprint_mask = self.wcs.footprint_contains(self.truthcoords)
        self.pixel_coords = self.wcs.world_to_pixel(self.truthcoords)
        self.coords_in = (self.pixel_coords[0][self.footprint_mask],
                          self.pixel_coords[1][self.footprint_mask])

        self.bkgstat = bkgstat
        self.bkg_annulus = bkg_annulus
        self._localbkg = LocalBackground(min(self.bkg_annulus),max(self.bkg_annulus),self.bkgstat)
        self.localbkg = self._localbkg(data=self.data,x=self.coords_in[0],
                                                      y=self.coords_in[1])

        self.band = band
        self.pointing = pointing
        self.chip = chip


    def set_bkgstat(self,new_stat):
        """If you don't want to use MMMBackground for the background
        statistics, you can set a different background estimator from
        photutils. Accepted background classes are:
        photutils.background.MeanBackground
        photutils.background.MedianBackground
        photutils.background.ModeEstimatorBackground
        photutils.background.MMMBackground
        photutils.background.SExtractorBackground
        photutils.background.BiweightLocationBackground

        But I only tested MMMBackground. 
        
        :param new_stat: background estimator from photutils.
        :type new_stat: photutils background object
        """
        self.bkgstat = new_stat