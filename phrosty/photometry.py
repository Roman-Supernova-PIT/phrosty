# IMPORTS Standard:
import numpy as np
import matplotlib.pyplot as plt

# IMPORTS Astro:
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, hstack, join
from astropy.visualization import simple_norm
import astropy.units as u
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats
from photutils.background import LocalBackground, MMMBackground, Background2D
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFBuilder, extract_stars, PSFPhotometry
from galsim import roman

# IMPORTS Internal:
from .utils import get_object_instances 

roman_bands = ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'K213']

def ap_phot(scienceimage,coords,wcs,
            ap_r=3, bkg_estimator=MMMBackground(), box_size=(50,50),
            filter_size=(3,3), method='subpixel',subpixels=5,
            merge_tables=True):
    """Does aperture photometry on the input science image and 
        specified coordinates. 

    :param scienceimage: Array containing science image.
    :type scienceimage: array-like
    :param coords: Table with columns 'x' and 'y' representing pixel
                                    coordinates for the detected sources. Table can contain
                                    more than just 'x' and 'y', and this may be useful if you 
                                    use merge_tables=True. 
    :type coords: astropy.table.Table
    :param ap_r: Aperture radius to use for aperture photometry. Defaults to 3.
    :type ap_r: int, optional
    :param bkg_estimator: _description_. Defaults to MMMBackground().
    :type bkg_estimator: _type_, optional
    :param box_size: _description_. Defaults to (50,50).
    :type box_size: tuple, optional
    :param filter_size: _description_. Defaults to (3,3).
    :type filter_size: tuple, optional
    :param method: _description_. Defaults to 'subpixel'.
    :type method: str, optional
    :param subpixels: _description_. Defaults to 5.
    :type subpixels: int, optional
    :param merge_tables: If true, output is merged coords and results from aperture
                        photometry. Defaults to True.
    :type merge_tables: boolean, optional

    Returns:
        astropy.table.Table: Table containing results from aperture photometry. 
    """
    x = np.array(coords['x'])
    y = np.array(coords['y'])
    photcoords = np.transpose(np.vstack([x,y]))
    apertures = CircularAperture(photcoords, r=ap_r)

    bkg = Background2D(scienceimage, box_size=box_size, filter_size=filter_size, bkg_estimator=bkg_estimator)
    bkgimg = bkg.background
    img_sub = scienceimage - bkgimg

    ap_results = aperture_photometry(img_sub,apertures,method=method,subpixels=subpixels)
    apstats = ApertureStats(scienceimage, apertures)
    ap_results['max'] = apstats.max

    # Needs to be 'xcentroid' and 'ycentroid' for PSF photometry. 
    ap_results.rename_column('xcenter','xcentroid')
    ap_results.rename_column('ycenter','ycentroid')
        
    for col in ap_results.colnames:
        ap_results[col] = ap_results[col].value

    if merge_tables:
        ap_results = hstack([ap_results,coords])
        # These are duplicates: 
        ap_results.remove_columns(['x','y'])
            
    return ap_results

def psf_phot(scienceimage,coords,wcs, 
            bkg_annulus=(50.0,80.0), saturation=0.9e5, noise=0.5e5,
            fwhm=3.0, fit_shape=(5,5), oversampling=3, maxiters=10, 
            exclude_duplicates=False, plot_epsf=False, 
            ap_r=3, bkg_estimator=MMMBackground(), box_size=(50,50),
            filter_size=(3,3), method='subpixel',subpixels=5):

    mean, median, stddev = sigma_clipped_stats(scienceimage)
    daofind = DAOStarFinder(fwhm=fwhm,threshold = 5.*(stddev))
    
    ap_results = ap_phot(scienceimage,coords,wcs,
                        ap_r=ap_r, bkg_estimator=bkg_estimator,
                        box_size=box_size, filter_size=filter_size,
                        method=method, subpixels=subpixels, merge_tables=True)

    psfstars = Table({'x': ap_results['xcentroid'], 'y': ap_results['ycentroid'],
                        'flux': ap_results['aperture_sum'], 'max': ap_results['max']})
    # NOTE: Need to make star and galaxy separation work in order to make this work. 
    print('len psfstars before saturation and flux', len(psfstars))
    psfstars = psfstars[psfstars['max'] < saturation]
    psfstars = psfstars[psfstars['flux'] > noise]
    print('len psfstars after saturation and flux', len(psfstars))

    stampsize=25
    median_subtracted_data = scienceimage - median
    nddata = NDData(data=median_subtracted_data)
    extracted_stars = extract_stars(nddata, psfstars, size=stampsize)

    if exclude_duplicates:
        # Get rid of stamps with more than one source.
        exclude_coords = []
        for i in range(len(extracted_stars)):
            try:
                stampsources = daofind(extracted_stars[i] - median)
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
    print('number of extracted stars:', len(extracted_stars))
    epsf_builder = EPSFBuilder(oversampling=oversampling, maxiters=maxiters)
    psf_func, fitted_stars = epsf_builder(extracted_stars)

    if plot_epsf:
        norm = simple_norm(psf_func.data, 'log', percent=99.0)
        plt.imshow(psf_func.data, norm=norm, origin='lower', cmap='Greys')
        plt.colorbar()
        plt.title('ePSF')
        plt.show()

    _localbkg = LocalBackground(min(bkg_annulus),max(bkg_annulus),bkg_estimator)
    localbkg = _localbkg(data=scienceimage,x=ap_results['xcentroid'],y=ap_results['ycentroid'])
    psfphot = PSFPhotometry(psf_func, fit_shape, localbkg_estimator=_localbkg,
                            finder=daofind, aperture_radius=ap_r)

    psf_results = psfphot(scienceimage, init_params=ap_results)
    psf_results['localbkg'] = localbkg

    radec = wcs.pixel_to_world(psf_results['x_fit'], psf_results['y_fit'])

    psf_results['ra'] = radec.ra.value
    psf_results['dec'] = radec.dec.value

    return psf_results

def crossmatch(pi,ti,seplimit=0.1):
    """ Cross-match the truth files from each image (TI) to the corresponding photometry
    file from that image (PI).

    :param ti: Astropy table from a truth file ('Roman_TDS_index_{band}_{pointing}_{sca}.txt')
    :type ti: Astropy table
    :param pi: Astropy table from a photometry file generated from the image with the same band,
                pointing, and sca as ti. 
    :type pi: Astropy table

    :return: Joined truth catalog and measured photometry catalog, 
            so that the measured objects are correctly crossmatched
            to their corresponding rows in the truth catalog. 
    :rtype: astropy.table 
    """

    if 'ra_truth' not in ti.colnames:
        ti['ra'].name = 'ra_truth'
    
    if 'dec_truth' not in ti.colnames:
        ti['dec'].name = 'dec_truth'

    if 'flux_truth' not in ti.colnames:
        ti['flux'].name = 'flux_truth'

    if 'mag_truth' not in ti.colnames:
        ti['mag'].name = 'mag_truth'

    tc = SkyCoord(ra=ti['ra_truth']*u.degree, dec=ti['dec_truth']*u.degree)
    pc = SkyCoord(ra=pi['ra']*u.degree, dec=pi['dec']*u.degree)
    ti_idx, pi_idx, angsep, dist3d = search_around_sky(tc,pc,seplimit=seplimit*(u.arcsec))

    ti_reduced = ti[ti_idx]
    pi_reduced = pi[pi_idx]

    ti_pi_reduced = hstack([ti_reduced,pi_reduced], join_type='exact')
    ti_x_pi = join(ti,ti_pi_reduced,join_type='outer')
    
    return ti_x_pi

def zpt_object_to_truth(object_data,ra,dec,
                        colnames=[['object_id','ra','dec','mag_truth','flux_truth','flux_fit','flux_err']]):

# Truth mag -2.5*np.log10(object_tab['flux_truth']) + 2.5*np.log10(exptime[band]*area_eff) + galsim_zp
    object_data['mag_truth'] = -2.5*np.log10(object_data['flux_truth']) + 2.5*np.log10()


def convert_flux_to_mag(ti_x_pi, band, config, zpt=False):
    """
    Input the astropy table from crossmatch.
    decide if you want to use this function or the one above!!!!!!!!!!!!
    """

    exptime = {'F184': 901.175,
               'J129': 302.275,
               'H158': 302.275,
               'K213': 901.175,
               'R062': 161.025,
               'Y106': 302.275,
               'Z087': 101.7}

    ti_x_pi['mag_fit'] = -2.5*np.log10(ti_x_pi['flux_fit'])
    ti_x_pi['mag_err'] = np.sqrt((1.09/ti_x_pi['flux_fit'])**2*ti_x_pi['flux_err']**2)

    if zpt:
        ti_x_pi['zpt'] = np.zeros(len(ti_x_pi))
        area_eff = roman.collecting_area
        galsim_zp = roman.getBandpasses()[band].zeropoint

        ti_x_pi['truth_mag'] = -2.5*np.log10(ti_x_pi['flux_truth']) + 2.5*np.log10(exptime[band]*area_eff) + galsim_zp
        
        for i, row in enumerate(ti_x_pi):
            objtab = get_object(row['object_id'], config)
            mean_mag = np.nanmean(objtab['mag_fit'])
            print('mean mag')
            print(mean_mag)
            zpt = np.unique(objtab['truth_mag'] - mean_mag)
            ti_x_pi['zpt'][i] = zpt
            ti_x_pi['mag_fit'] += zpt 

    return ti_x_pi

        #     if zpt == 'truth':
#         ap_zpt_mask = np.logical_and(results_table[f'{self.band}_ap_mag']>-11, results_table[f'{self.band}_ap_mag']<-9)
#         psf_zpt_mask = np.logical_and(results_table[f'{self.band}_psf_mag']>-11, results_table[f'{self.band}_psf_mag']<-9)

#         truthmag = truth_table[self.band][self.footprint_mask]
#         results_table[f'{self.band}_truth'] = truthmag
#         ap_zpt = np.median(results_table[f'{self.band}_ap_mag'][ap_zpt_mask] - truthmag[ap_zpt_mask])
#         psf_zpt = np.median(results_table[f'{self.band}_psf_mag'][psf_zpt_mask] - truthmag[psf_zpt_mask])
        
#         results_table[f'{self.band}_ap_mag'] -= ap_zpt
#         results_table[f'{self.band}_psf_mag'] -= psf_zpt