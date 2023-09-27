import numpy as np
from scipy import stats
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
from photutils.psf import EPSFBuilder, extract_stars, PSFPhotometry, IntegratedGaussianPRF

roman_bands = ['Y106', 'J129', 'H158', 'F184']

class truth():
    # https://roman.ipac.caltech.edu/data/sims/sn_image_sims/galaxia_akari.fits.gz
    def __init__(self,filepath):
        """filepath should be to a truth file that only contains stars.
        No galaxies! 

        :param filepath: Path to fits file. 
        :type filepath: str
        """
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
        truth_mags = {}
        for band in roman_bands:
            truth_mags[band] = truth[band]
        truth_tab = Table([np.arange(0,len(truth),1),truth_ra,truth_dec], names=('index','ra','dec'))
        truth_tab |= truth_mags
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
        
        if band not in roman_bands:
            raise ValueError(f'Argument "band" must be in {roman_bands}.')

        if pointing == 'unspecified':
            raise ValueError('You need to specify a pointing ID.')

        if chip == 'unspecified':
            raise ValueError('You need to specify a chip number.')

        self.filepath = filepath
        self.sciencehdu = fits.open(self.filepath)
        self.data = self.sciencehdu[1].data
        self.wcs = WCS(self.sciencehdu[1].header)
        self.mean, self.median, self.std = sigma_clipped_stats(self.data, sigma=3.0)

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

        self.ap_r = 3.0

        self.ap_phot_results = None
        self.psf_phot_results = None

    def plot_truth(self):
        zscale=ZScaleInterval()
        z1,z2 = zscale.get_limits(self.data)
        plt.figure(figsize=(8,8))
        plt.imshow(self.data,vmin=z1,vmax=z2,cmap='Greys',origin='lower')
        plt.plot(self.coords_in[0],self.coords_in[1],linestyle='',marker='o',
                fillstyle='none',markeredgecolor='red',alpha=0.5)
        plt.xlabel('x [px]')
        plt.ylabel('y [px]')
        plt.colorbar()
        plt.show()

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
        self._localbkg = LocalBackground(min(self.bkg_annulus),max(self.bkg_annulus),self.bkgstat)
        self.localbkg = self._localbkg(data=self.data,x=self.coords_in[0],
                                                      y=self.coords_in[1])
        
    def ap_phot(self,ap_r,method='subpixel',subpixels=5):
        """Do basic aperture photometry. 

        :param ap_r: Circular aperture radius.
        :type ap_r: float
        :param method: See documentation for photutils.aperture.aperture_photometry, defaults to 'subpixel'
        :type method: str, optional
        :param subpixels: See documentation for photutils.aperture.aperture_photometry, defaults to 5
        :type subpixels: int, optional
        :return: Table with results from aperture photometry.
        :rtype: astropy Table object
        """        
        self.ap_r = ap_r
        bkgval = self.bkgstat.calc_background(self.data)
        apertures = CircularAperture(np.transpose(self.coords_in), r=ap_r)
        ap_results = aperture_photometry(self.data-bkgval,apertures,method=method,subpixels=subpixels)
        apstats = ApertureStats(self.data, apertures)
        ap_results['max'] = apstats.max
        ap_results['x_init'], ap_results['y_init'] = ap_results['xcenter'].value, ap_results['ycenter'].value
        ap_results[self.band] = -2.5*np.log10(ap_results['aperture_sum'].value)
        
        self.ap_phot_results = ap_results

        return ap_results
    
    def psf_phot(self, psf, ap_results=None, saturation=8e4, noise=10**4.2,
                fwhm=3.0, fit_shape=(5,5), method='subpixel', subpixels=5,
                oversampling=3, maxiters=10,plot_epsf=False):
        """Do PSF photometry. This is essentially a wrapper for photutils.psf.PSFPhotometry.
        Right now, because we're using the truth coordinates, we use PSFPhotometry, not
        IterativePSFPhotometry. 

        :param psf: Type of PSF to use. Can either assume a Gaussian PSF ('gaussian'),
                        or construct PSF from field stars ('ePSF'). 
        :type psf: string
        :param ap_results: Output directly from scienceimg.ap_phot(). If this is not
                            specified, ap_phot() will be rerun. defaults to None
        :type ap_results: astropy Table object, optional
        :param saturation: Maximum pixel value to use for choosing stars to generating ePSF. 
                            Only necessary if psf='epsf', defaults to 8e4
        :type saturation: float, optional
        :param noise: Noise floor to use for choosing stars to generate ePSF.
                        Only necessary if psf='epsf', defaults to 10**4.2
        :type noise: float, optional
        :param fwhm: _description_, defaults to 3.0
        :type fwhm: float, optional
        :param fit_shape: _description_, defaults to (5,5)
        :type fit_shape: tuple, optional
        :param method: _description_, defaults to 'subpixel'
        :type method: str, optional
        :param subpixels: _description_, defaults to 5
        :type subpixels: int, optional
        :raises ValueError: _description_
        """        
        daofind = DAOStarFinder(fwhm=fwhm,threshold = 5.*(self.std))

        psf_types = ['gaussian','epsf']
        if ap_results is None:
            ap_results = self.ap_phot(self.ap_r,method=method,subpixels=subpixels)

        if psf not in psf_types:
            raise ValueError(f'Argument psf must be in {psf_types}.')
        elif psf == 'gaussian':
            psf_func = IntegratedGaussianPRF(sigma=fwhm)
        elif psf == 'epsf':
            psfstars = Table({'x': self.coords_in[0], 'y': self.coords_in[1],
                                'flux': ap_results['aperture_sum'], 'max': ap_results['max']})
            psfstars = psfstars[psfstars['max'] < saturation]
            psfstars = psfstars[psfstars['flux'] > noise]

            stampsize=25
            median_subtracted_data = self.data-self.median
            nddata = NDData(data=median_subtracted_data,wcs=self.wcs)
            extracted_stars = extract_stars(nddata, psfstars, size=stampsize)

            # Get rid of stamps with more than one source.
            exclude_coords = []
            for i in range(len(extracted_stars)):
                try:
                    stampsources = daofind(extracted_stars[i] - self.median)
                    if len(stampsources) > 1 or len(stampsources) < 1:
                        exclude_coords.append(extracted_stars.center_flat[i])
                except:
                    pass

            exclude_rows = []
            for c in exclude_coords:
                exclude_rows.append(psfstars[psfstars['x'] == c[0]])

            new_psfstars_rows = [x for x in psfstars if x not in exclude_rows]
            new_psfstars = Table(rows=new_psfstars_rows, names=psfstars.colnames)
            extracted_stars = extract_stars(nddata, new_psfstars, size=stampsize)

            # Build ePSF.
            epsf_builder = EPSFBuilder(oversampling=oversampling, maxiters=maxiters)
            psf_func, fitted_stars = epsf_builder(extracted_stars)

            if plot_epsf:
                norm = simple_norm(psf_func.data, 'log', percent=99.0)
                plt.imshow(psf_func.data, norm=norm, origin='lower', cmap='Greys')
                plt.colorbar()
                plt.title('ePSF')
                plt.show()

        psfphot = PSFPhotometry(psf_func, fit_shape, localbkg_estimator=self._localbkg,
                                finder=daofind, aperture_radius=self.ap_r)
        psf_results = psfphot(self.data, init_params=ap_results)
    
        # This will only return one table of values. So, we will probably want to
        # run this a bunch of times, save the outputs to csv files, and then run the
        # coordinate-matching stuff on the saved files. We really don't want a list
        # of tables stored in temporary memory. That's not great. 

        self.psf_phot_results = psf_results

        return psf_results
    
    def save_results(self,savepath,overwrite=False,truth_zpt=True,
                     truth_table=None):
        """_summary_

        :param savepath: path to save file. Include filename.
        :type savepath: str
        :param overwrite: Set True if you want to overwrite previously saved files, defaults to False
        :type overwrite: bool, optional
        :param truth_zpt: Zero point photometry to truth table, defaults to True
        :type truth_zpt: bool, optional
        :param truth_table: _description_, defaults to None
        :type truth_table: _type_, optional
        :raises Exception: _description_
        :raises ValueError: _description_
        """                     
        if self.ap_phot_results is None or self.psf_phot_results is None:
            raise Exception('You need to run either self.ap_phot() or self.psf_phot() to have results to save!')
            # If you ran self.psf_phot(), then self.ap_phot() was run. :)

        if truth_zpt == True and truth_table is None:
            raise ValueError('If truth_zpt==True, then you must provide a table from truth.table() in argument truth_table.')

        results_table = Table()
        index = [i for i in range(len(self.ap_phot_results))]
        results_table['index'] = index
        results_table['source_ID'] = [f'{self.band}_{self.pointing}_{self.chip}_{i}' for i in range(len(self.ap_phot_results))]
        results_table['x_init'], results_table['y_init'] = self.ap_phot_results['x_init'], self.ap_phot_results['y_init']
        
        ra_dec = self.wcs.pixel_to_world(self.ap_phot_results['x_init'], self.ap_phot_results['y_init'])
        results_table['ra'], results_table['dec'] = ra_dec.ra.value, ra_dec.dec.value
        
        results_table[f'{self.band}_ap_flux'] = self.ap_phot_results['aperture_sum']
        results_table[f'{self.band}_ap_mag'] = -2.5*np.log10(self.ap_phot_results['aperture_sum'])

        if self.psf_phot_results is not None:
            results_table[f'{self.band}_psf_flux'] = self.psf_phot_results['flux_fit']
            results_table[f'{self.band}_psf_flux_err'] = self.psf_phot_results['flux_err']
            results_table[f'{self.band}_psf_mag'] = -2.5*np.log10(results_table[f'{self.band}_psf_flux'])
            results_table[f'{self.band}_psf_mag_err'] = np.sqrt((1.09/results_table[f'{self.band}_psf_flux'])**2*results_table[f'{self.band}_psf_flux_err']**2)

        if truth_zpt == True:
            ap_zpt_mask = np.logical_and(results_table[f'{self.band}_ap_mag']>-11, results_table[f'{self.band}_ap_mag']<-9)
            psf_zpt_mask = np.logical_and(results_table[f'{self.band}_psf_mag']>-11, results_table[f'{self.band}_psf_mag']<-9)

            truthmag = truth_table[self.band][self.footprint_mask]
            ap_zpt = np.median(results_table[f'{self.band}_ap_mag'][ap_zpt_mask] - truthmag[ap_zpt_mask])
            psf_zpt = np.median(results_table[f'{self.band}_psf_mag'][psf_zpt_mask] - truthmag[psf_zpt_mask])

            results_table[f'{self.band}_ap_mag'] -= ap_zpt
            results_table[f'{self.band}_psf_mag'] -= psf_zpt

        results_table['band'] = self.band
        results_table['pointing'] = self.pointing
        results_table['chip'] = self.chip

        results_table.write(savepath, format='csv', overwrite=overwrite)

def crossmatch_truth(truth_filepath,results_filepaths,seplimit=0.1,psf=True):
    """
    This will handle a list of science files with mixed bands! Hopefully!

    :param truth_filepath: Filepath to truth file.
    :type truth_filepath: str
    :param results_filepaths: List of filepaths to csv files from scienceimg.save_results().
    :type results_filepaths: list
    :param seplimit: Angular separation limit to find matches, defaults to 0.1
    :type seplimit: float, optional
    :param psf: Include PSF photometry results in the matching table, defaults to True
    :type psf: bool, optional
    :return: _description_
    :rtype: _type_
    """    
    prefixes = ['ap']
    if psf:
        prefixes.append('psf')
        
    suffixes = ['_flux', '_flux_err', '_mag', '_mag_err']
    match_vals = []
    for b in roman_bands:
        for p in prefixes:
            for s in suffixes:
                if p != 'ap' and s != '_flux_err':
                    match_vals.append(b+'_'+p+s)
                else:
                    pass

    tr = truth(truth_filepath)
    tr_tab = deepcopy(tr.table())
    for i, file in enumerate(results_filepaths):
        check = Table.read(file, format='csv')
        band = check['band'][0]
        check_coords = SkyCoord(ra=check['ra']*u.degree, dec=check['dec']*u.degree)
        tr_coords = SkyCoord(ra=tr_tab['ra'], dec=tr_tab['dec']) # Already in degrees
        check_idx, tr_idx, angsep, dist3d = search_around_sky(check_coords,tr_coords,seplimit=seplimit*u.arcsec)

        # Create a mask so values go in correct rows.
        tr_mask = np.empty(len(tr_tab))*np.nan
        for val in match_vals:
            if band in val:
                tr_mask[tr_idx] = check[val]
                tr_tab[f'{val}_{i}'] = np.empty(len(tr_tab))*np.nan
                tr_tab[f'{val}_{i}'] = tr_mask
    
    # Put all the matching data in the same column. 
    n_tables = len(results_filepaths)
    matchescol = {x: [] for x in match_vals}

    # There absolutely must be a way to merge those columns with
    # array operations instead of loops, but I can't make it behave.
    # Also, for some reason it puts an empty list in the cells where
    # I input band data, but leaves cells blank if I don't put in the 
    # band data. Which I understand sounds correct, but it's not; empty
    # cells should be empty in the same way, uniformly. I guess TBD
    # on if this is a problem. 
    for key in matchescol:
        for row in tr_tab:
            vals = []
            for n in range(n_tables):
                try:
                    if not np.isnan(row[f'{key}_{n}']):
                        vals.append(row[f'{key}_{n}'])
                except:
                    pass
            matchescol[key].append(vals)

    tr_tab |= matchescol

    return tr_tab