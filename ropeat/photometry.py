import numpy as np
import os
from operator import methodcaller
from collections import defaultdict
from scipy import stats
from copy import deepcopy
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.convolution import convolve
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, hstack, MaskedColumn
from astropy.visualization import ZScaleInterval, simple_norm
from astropy.wcs import WCS
import astropy.units as u
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats
from photutils.background import LocalBackground, MMMBackground, Background2D
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFBuilder, extract_stars, PSFPhotometry, IntegratedGaussianPRF
from photutils.segmentation import make_2dgaussian_kernel, detect_sources, deblend_sources, SourceCatalog
from galsim import roman

roman_bands = ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'W146', 'K213']

class truth():
    # https://roman.ipac.caltech.edu/data/sims/sn_image_sims/galaxia_akari.fits.gz
    def __init__(self,filepath):
        """filepath should be to a truth file that only contains stars.
        No galaxies! The fits file should have columns ['ra', 'dec', 'Y106', 'J129', 'H158', 'F184'].

        :param filepath: Path to fits file. 
        :type filepath: str
        """
        self.filepath = filepath
        self.truthhdu = fits.open(self.filepath)

    def table(self,units='degrees'):
        """Returns an astropy table with truth data. 
        :param units: Input units of ra and dec. Currently accepts ['degrees','radians'].
        :type units: str        
        :return: astropy table
        :rtype: astropy table
        """
        
        acceptable_units = ['degrees','radians']
        if units not in acceptable_units:
            raise ValueError(f'Argument units must be in {acceptable_units}.')
            
        truth_tab = Table(self.truthhdu[1].data)
        
        if units=='degrees':
            truth_ra = truth_tab['ra']*u.deg
            truth_dec = truth_tab['dec']*u.deg
        elif units=='radians':
            truth_ra = (truth_tab['ra']*u.radian).to(u.deg)
            truth_dec = (truth_tab['dec']*u.radian).to(u.deg)

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
    :param pointing: Pointing number. Used for assigning a source ID in the output table.
    :type pointing: str, float, int
    :param chip: Chip/SCA number. Used for assigning a source ID in the output table.
    :type chip: str, float, int
    :param bkgstat: Background estimator from photutils. If you don't want to use MMMBackground for the background
        statistics, you can set a different background estimator from
        photutils. Accepted background classes are:
        photutils.background.MeanBackground
        photutils.background.MedianBackground
        photutils.background.ModeEstimatorBackground
        photutils.background.MMMBackground
        photutils.background.SExtractorBackground
        photutils.background.BiweightLocationBackground
    :type bkgstat: photutils background object
    :param bkg_annulus: Define annulus for computing local background with (inner radius, outer radius) in pixels. 
    :type bkg_annulus: tuple
    """
    def __init__(self,filepath,
                 truthcoords,
                 band='unspecified',
                 pointing='unspecified',
                 chip='unspecified',
                 bkgstat=MMMBackground(),
                 bkg_annulus=(50.0,80.0),
                 box_size=(50,50),
                 filter_size=(3,3)
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
        self.initcoords = None
        self.footprint_mask = self.wcs.footprint_contains(self.truthcoords)
        self.pixel_coords = self.wcs.world_to_pixel(self.truthcoords)
        self.coords_in = (self.pixel_coords[0][self.footprint_mask],
                          self.pixel_coords[1][self.footprint_mask])

        self.bkgstat = bkgstat
        self.bkg = Background2D(self.data, box_size, filter_size=filter_size, bkg_estimator=self.bkgstat)
        self.bkgimg = self.bkg.background
        self.bkg_annulus = bkg_annulus
        self._localbkg = LocalBackground(min(self.bkg_annulus),max(self.bkg_annulus),self.bkgstat)
        self.localbkg = self._localbkg(data=self.data,x=self.coords_in[0],
                                                      y=self.coords_in[1])
        self.bkg_sub_data = self.data - self.bkgimg
        
        self.convolution_kernel = make_2dgaussian_kernel(fwhm=3.0, size=5)
        self.convolved_data = None   

        self.band = band
        self.pointing = pointing
        self.chip = chip

        self.ap_r = 3.0

        self.ap_phot_results = None
        self.psf_phot_results = None
        
    def plot_img(self):
        """
        Plot the science image. 
        
        """

        zscale=ZScaleInterval()
        z1,z2 = zscale.get_limits(self.data)
        plt.figure(figsize=(8,8))
        plt.imshow(self.data,vmin=z1,vmax=z2,cmap='Greys',origin='lower')
        plt.xlabel('x [px]')
        plt.ylabel('y [px]')
        plt.colorbar()
        plt.show()
        
    def plot_truth(self):
        """"
        Plot the locations of the truth coordinates over the science image. 
        """
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

    def ap_phot(self,ap_r,method='subpixel',subpixels=5,truthcoordsonly=False):
        """Do basic aperture photometry and get morphological properties of the data. 

        :param ap_r: Circular aperture radius.
        :type ap_r: float
        :param method: See documentation for photutils.aperture.aperture_photometry, defaults to 'subpixel'
        :type method: str, optional
        :param subpixels: See documentation for photutils.aperture.aperture_photometry, defaults to 5
        :type subpixels: int, optional
        :return: Table with results from aperture photometry.
        :rtype: astropy Table object
        """        
        if not truthcoordsonly:
            # First, detect sources and get morphology info. 
            convolved_data = convolve(self.bkg_sub_data, self.convolution_kernel) 
            threshold = 1.5*self.bkg.background_rms
            segment_map = detect_sources(convolved_data, threshold, npixels=10)
            cat = SourceCatalog(self.bkg_sub_data, segment_map, convolved_data=convolved_data, wcs=self.wcs)
            cattab = cat.to_table()
            cattab['cxx'] = cat.cxx
            cattab['cxy'] = cat.cxy
            cattab['cyy'] = cat.cyy
            cattab['ellipticity'] = cat.ellipticity
            cattab['fwhm'] = cat.fwhm

            splitcoordsfunc = lambda x: x.to_string().split(' ')
            splitcoords = np.array(list(map(splitcoordsfunc, cattab['sky_centroid'])), dtype=float).T
            cattab['ra'] = splitcoords[0]
            cattab['dec'] = splitcoords[1]

            dropcols = ['kron_flux', 'kron_fluxerr', 'min_value', 'max_value', 'segment_flux', 'segment_fluxerr','local_background']
            cattab.remove_columns(dropcols)
        
        # Now, aperture photometry. 
        self.ap_r = ap_r
        if truthcoordsonly:
            self.initcoords = np.transpose(self.coords_in)
        elif not truthcoordsonly:
            x = np.array(cattab['xcentroid'])
            y = np.array(cattab['ycentroid'])
            self.initcoords = np.transpose(np.vstack([x,y]))
            
        apertures = CircularAperture(self.initcoords, r=ap_r)   
        ap_results = aperture_photometry(self.bkg_sub_data,apertures,method=method,subpixels=subpixels)
        apstats = ApertureStats(self.data, apertures)
        ap_results['max'] = apstats.max

        ap_results[f'{self.band}_ap_mag'] = -2.5*np.log10(ap_results['aperture_sum'].value)
        if not truthcoordsonly:
            self.ap_phot_results = hstack([ap_results, cattab])
        elif truthcoordsonly:
            ap_results.rename_column('xcenter','xcentroid')
            ap_results.rename_column('ycenter','ycentroid')
            
            for col in ap_results.colnames:
                ap_results[col] = ap_results[col].value
                
            self.ap_phot_results = ap_results
        
        return self.ap_phot_results
    
    def psf_phot(self, psf, ap_results=None, saturation=8e4, noise=10**4.2,
                fwhm=3.0, fit_shape=(5,5), method='subpixel', subpixels=5,
                oversampling=3, maxiters=10, exclude_duplicates=False, plot_epsf=False):
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
            psfstars = Table({'x': self.initcoords.T[0], 'y': self.initcoords.T[1],
                                'flux': ap_results['aperture_sum'], 'max': ap_results['max']})
            psfstars = psfstars[psfstars['max'] < saturation]
            psfstars = psfstars[psfstars['flux'] > noise]

            stampsize=25
            median_subtracted_data = self.data-self.median
            nddata = NDData(data=median_subtracted_data,wcs=self.wcs)
            extracted_stars = extract_stars(nddata, psfstars, size=stampsize)

            if exclude_duplicates:
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
    
    def save_results(self,savepath,zpt,overwrite=False,
                     truth_table=None):
        """Save photometry results to a csv. 

        :param savepath: path to save file. Include filename.
        :type savepath: str
        :param overwrite: Set True if you want to overwrite previously saved files, defaults to False
        :type overwrite: bool, optional
        :param zpt: Zero point ['galsim','truth',None]
        :type zpt: str or NoneType, optional
        :param truth_table: _description_, defaults to None
        :type truth_table: _type_, optional
        :raises Exception: _description_
        :raises ValueError: _description_
        """                     
        if self.ap_phot_results is None or self.psf_phot_results is None:
            raise Exception('You need to run either self.ap_phot() or self.psf_phot() to have results to save!')
            # If you ran self.psf_phot(), then self.ap_phot() was run. :)

        if zpt == 'truth' or zpt == 'galsim' and truth_table is None:
            raise ValueError('If zpt == "truth", then you must provide a table from truth.table() in argument truth_table.')

        results_table = Table()
        index = [i for i in range(len(self.ap_phot_results))]
        results_table['index'] = index
        results_table['source_ID'] = [f'{self.band}_{self.pointing}_{self.chip}_{i}' for i in range(len(self.ap_phot_results))]
        
        
        for col in self.ap_phot_results.colnames:
            results_table[col] = self.ap_phot_results[col]
        
        x_truth, y_truth = self.coords_in[0], self.coords_in[1]
        results_table['x_truth'], results_table['y_truth'] = x_truth, y_truth
        results_table['x_fit'], results_table['y_fit'] = self.psf_phot_results['x_fit'], self.psf_phot_results['y_fit']
        
        ra_dec = self.wcs.pixel_to_world(self.psf_phot_results['x_fit'], self.psf_phot_results['y_fit'])
        results_table['ra'], results_table['dec'] = ra_dec.ra.value, ra_dec.dec.value
        
        results_table[f'{self.band}_ap_max'] = self.ap_phot_results['max']
        
        results_table[f'{self.band}_ap_flux'] = self.ap_phot_results['aperture_sum']

        if self.psf_phot_results is not None:
            results_table[f'{self.band}_psf_flux'] = self.psf_phot_results['flux_fit']
            results_table[f'{self.band}_psf_flux_err'] = self.psf_phot_results['flux_err']
            results_table[f'{self.band}_psf_mag'] = -2.5*np.log10(results_table[f'{self.band}_psf_flux'])
            results_table[f'{self.band}_psf_mag_err'] = np.sqrt((1.09/results_table[f'{self.band}_psf_flux'])**2*results_table[f'{self.band}_psf_flux_err']**2)

        if zpt == 'truth':
            ap_zpt_mask = np.logical_and(results_table[f'{self.band}_ap_mag']>-11, results_table[f'{self.band}_ap_mag']<-9)
            psf_zpt_mask = np.logical_and(results_table[f'{self.band}_psf_mag']>-11, results_table[f'{self.band}_psf_mag']<-9)

            truthmag = truth_table[self.band][self.footprint_mask]
            results_table[f'{self.band}_truth'] = truthmag
            ap_zpt = np.median(results_table[f'{self.band}_ap_mag'][ap_zpt_mask] - truthmag[ap_zpt_mask])
            psf_zpt = np.median(results_table[f'{self.band}_psf_mag'][psf_zpt_mask] - truthmag[psf_zpt_mask])
            
            results_table[f'{self.band}_ap_mag'] -= ap_zpt
            results_table[f'{self.band}_psf_mag'] -= psf_zpt
            
        elif zpt == 'galsim':            
            maglims = [-5,-7.5] # This is the median of the un-zeropointed data +/- 1 standard deviation, roughly. 
            truthmag = truth_table[self.band][self.footprint_mask]
            results_table[f'{self.band}_truth'] = truthmag
            
            ap_zpt_mask = np.logical_and(results_table[f'{self.band}_ap_mag'] > maglims[1],
                                          results_table[f'{self.band}_ap_mag'] < maglims[0])
            psf_zpt_mask = np.logical_and(results_table[f'{self.band}_psf_mag'] > maglims[1],
                                          results_table[f'{self.band}_psf_mag'] < maglims[0])
            
            ap_zpt = np.median(results_table[f'{self.band}_ap_mag'][ap_zpt_mask] - truthmag[ap_zpt_mask])
            psf_zpt = np.median(results_table[f'{self.band}_psf_mag'][psf_zpt_mask] - truthmag[psf_zpt_mask])
            
            results_table[f'{self.band}_ap_mag'] -= ap_zpt
            results_table[f'{self.band}_psf_mag'] -= psf_zpt
            
        elif zpt is None:
            truthmag = truth_table[self.band][self.footprint_mask]
            results_table[f'{self.band}_truth'] = truthmag
            
        results_table['band'] = self.band
        results_table['pointing'] = self.pointing
        results_table['chip'] = self.chip
        
        dropcols = ['id','label','xcenter','ycenter','xcentroid','ycentroid']
        for col in dropcols:
            try:
                results_table.remove_column(col)
            except:
                pass

        results_table['psfphot_flags'] = self.psf_phot_results['flags']

        results_table.write(savepath, format='csv', overwrite=overwrite)

def crossmatch_truth(truth_filepath,results_filepaths,savename,overwrite=True,seplimit=0.1,psf=True,verbose=True,temp_file_path=None):
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
    all_suffixes = []
    
    for p in prefixes:
        for s in suffixes:
            if p == 'ap' and 'err' in s:
                pass
            else:
                all_suffixes.append('_'+p+s)
            for b in roman_bands:
                if p == 'ap' and 'err' in s:
                    pass
                else:
                    match_vals.append(b+'_'+p+s)
                    
                if f'{b}_max' not in match_vals:
                    match_vals.append(f'{b}_max')

    match_vals.append('ra_all')
    match_vals.append('dec_all')
    match_vals.append('x_truth')
    match_vals.append('y_truth')
    match_vals.append('x_fit')
    match_vals.append('y_fit')
    match_vals.append('pointing_all')
    match_vals.append('sca_all')
    match_vals.append('psfphot_flags_all')
    
    tr = truth(truth_filepath)
    tr_tab = tr.table()
    if temp_file_path is None:
        temp_file_name = 'tempfile_DELETEME.fits'
    elif isinstance(temp_file_path, str) and '.fits' in temp_file_path:
        temp_file_name = temp_file_path
    else:
        print('Something is wrong with your temporary file path.')
        print('Maybe it doesnt end in .fits? Or is not a string?')

    for col in match_vals:
        tr_tab.add_column('empty', name=col)
        
    for col in tr_tab.colnames:
        if col in roman_bands:
            tr_tab.rename_column(col, f'{col}_truth')
            
    tr_tab.rename_column('ra','ra_truth')
    tr_tab.rename_column('dec','dec_truth')
    
    tr_tab.write(temp_file_name, format='fits', overwrite=True)

    if verbose:
        print(f'We need to get through matching {len(results_filepaths)} files.')
        nfiles = 0
        nfail = 0
    for i, file in enumerate(results_filepaths):
        if os.path.exists(file):
            if verbose:
                print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
                print(file)
            with fits.open(temp_file_name) as f:
                tr_tab = Table(f[1].data)
                check = Table.read(file, format='csv')
                band = check['band'][0]
                if verbose:
                    print('Succesfully opened file. Next, crossmatch and merge tables.')
                check_coords = SkyCoord(ra=check['ra']*u.degree, dec=check['dec']*u.degree)
                tr_coords = SkyCoord(ra=tr_tab['ra_truth']*u.degree, dec=tr_tab['dec_truth']*u.degree) 
                check_idx, tr_idx, angsep, dist3d = search_around_sky(check_coords,tr_coords,seplimit=seplimit*u.arcsec)

                tr_tab = tr_tab.to_pandas()
                appendvals = lambda x,y : str(x) + ',' + str(y)
                
                # Collect all magnitudes into a string into one column. 
                for s in all_suffixes:
                    tlist = list(tr_tab[f'{band}{s}'])
                    tlist_reduced = list(np.array(tlist)[tr_idx])
                    clist = list(check[f'{band}{s}'])
                    clist_reduced = list(np.array(clist)[check_idx])
                    strcol = list(map(appendvals,tlist_reduced,clist_reduced))
                    strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
                    tr_tab.loc[tr_idx, f'{band}{s}'] = strcol
                    
                # Collect RA/dec into a string into one column. 
                for c in ['ra','dec']:
                    tlist = list(tr_tab[f'{c}_all'])
                    tlist_reduced = list(np.array(tlist)[tr_idx])
                    clist = list(check[c])
                    clist_reduced = list(np.array(clist)[check_idx])
                    strcol = list(map(appendvals,tlist_reduced,clist_reduced))
                    strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
                    tr_tab.loc[tr_idx, f'{c}_all'] = strcol
                    
                # Collect all x/y positions into a string into one column.
                for c in ['x','y']:
                    for a in ['fit','truth']:
                        tlist = list(tr_tab[f'{c}_{a}'])
                        tlist_reduced = list(np.array(tlist)[tr_idx])
                        clist = list(check[f'{c}_{a}'])
                        clist_reduced = list(np.array(clist)[check_idx])
                        strcol = list(map(appendvals,tlist_reduced,clist_reduced))
                        strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
                        tr_tab.loc[tr_idx, f'{c}_{a}'] = strcol
                    
                # Collect the max. flux pixel value in an aperture into one column. 
                tlist = list(tr_tab[f'{band}_max'])
                tlist_reduced = list(np.array(tlist)[tr_idx])
                clist = list(check['max'])
                clist_reduced = list(np.array(clist)[check_idx])
                strcol = list(map(appendvals,tlist_reduced,clist_reduced))
                strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
                tr_tab.loc[tr_idx, f'{band}_max'] = strcol
                    
                # Collect all SCA IDs into a string in one column.
                tlist = list(tr_tab['sca_all'])
                tlist_reduced = list(np.array(tlist)[tr_idx])
                clist = list(check['chip'])
                clist_reduced = list(np.array(clist)[check_idx])
                strcol = list(map(appendvals,tlist_reduced,clist_reduced))
                strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
                tr_tab.loc[tr_idx, 'sca_all'] = strcol
                
                # Collect all pointings into a string in one column.
                tlist = list(tr_tab['pointing_all'])
                tlist_reduced = list(np.array(tlist)[tr_idx])
                clist = list(check['pointing'])
                clist_reduced = list(np.array(clist)[check_idx])
                strcol = list(map(appendvals,tlist_reduced,clist_reduced))
                strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
                tr_tab.loc[tr_idx, 'pointing_all'] = strcol
                
                # Collect all psfphot flags into a string in one column.
                tlist = list(tr_tab['psfphot_flags_all'])
                tlist_reduced = list(np.array(tlist)[tr_idx])
                clist = list(check['psfphot_flags'])
                clist_reduced = list(np.array(clist)[check_idx])
                strcol = list(map(appendvals,tlist_reduced,clist_reduced))
                strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
                tr_tab.loc[tr_idx, 'psfphot_flags_all'] = strcol

            tr_tab = Table.from_pandas(tr_tab)
            tr_tab.write(temp_file_name, format='fits', overwrite=True)
            if verbose:
                print(f'Wrote the crossmatched data to the main catalog in a temporary file, {temp_file_name}.')
                nfiles += 1
                print(f'Made it through crossmatching {nfiles}/{len(results_filepaths)} files.')

        else:
            if verbose:
                print(f'Oops! {file} does not seem to exist.')
                nfail += 1
                print(f'In total, {nfail}/{len(results_filepaths)} filepaths have failed.')
        
    # Drop empty columns. 
    for col in [c for c in tr_tab.colnames if c not in ['object_id','config']]:
        try:
            tr_tab[col] = MaskedColumn(data=tr_tab[col].value, dtype=np.float64)
            tr_tab[col].mask = np.isinf(tr_tab[col].value)
            if all(tr_tab[col].mask):
                tr_tab.remove_column(col)

        except:
            if all(tr_tab[col] == 'empty'):
                tr_tab.remove_column(col)
            
    # Write table to file. 
    tr_tab.write(savename, format='csv', overwrite=overwrite)
    if verbose:
        print('Final crossmatched file is written.')
        print('Finally, deleting the temporary file.')
    os.remove(temp_file_name)