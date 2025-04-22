# IMPORTS Standard:
import numpy as np
import matplotlib.pyplot as plt

# IMPORTS Astro:
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, hstack, join
from astropy.visualization import simple_norm
from astropy.modeling.fitting import NonFiniteValueError
import astropy.units as u
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats
# from photutils.background import LocalBackground, MMMBackground, Background2D
from photutils.detection import DAOStarFinder
from photutils.psf import PSFPhotometry, FittableImageModel # EPSFBuilder, extract_stars, 
from photutils.background import LocalBackground, MMMBackground
from galsim import roman


"""
NOTE: This module assumes images have already had background subtracted. 
"""

roman_bands = ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'K213']

def ap_phot(scienceimage,coords,
            ap_r=9, method='subpixel',subpixels=5,
            merge_tables=True,**kwargs):
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

    ap_results = aperture_photometry(scienceimage,apertures,method=method,subpixels=subpixels,**kwargs)
    apstats = ApertureStats(scienceimage, apertures)
    ap_results['max'] = apstats.max

    # Needs to be 'xcentroid' and 'ycentroid' for PSF photometry. 
    # Same with 'flux'. 
    ap_results.rename_column('xcenter','xcentroid')
    ap_results.rename_column('ycenter','ycentroid')
    # ap_results.rename_column('aperture_sum','flux')

    for col in ap_results.colnames:
        ap_results[col] = ap_results[col].value

    if merge_tables:
        ap_results = hstack([ap_results,coords])
        # These are duplicates: 
        ap_results.remove_columns(['x','y'])

    return ap_results

def psfmodel(psfimg):
    psf = FittableImageModel(psfimg)
    return psf
    
def psf_phot(scienceimage,psf,init_params,wcs=None,
             forced_phot=True, fwhm=3.0, fit_shape=(5,5), 
             oversampling=3, maxiters=10):

    # mean, median, stddev = sigma_clipped_stats(scienceimage)
    # daofind = DAOStarFinder(fwhm=fwhm,threshold = 5.*(stddev))

    if 'flux_init' not in init_params.colnames:
        raise Exception('Astropy table passed to kwarg init_params must contain column \"flux_init\".')

    if forced_phot:
        print('x, y are fixed!')
        psf.x_0.fixed = True
        psf.y_0.fixed = True
    else:
        print('x, y are fitting parameters!')
        psf.x_0.fixed = False
        psf.x_0.fixed = False

    try: 
        bkgfunc = LocalBackground(15,25,MMMBackground())
        psfphot = PSFPhotometry(psf,fit_shape,localbkg_estimator=bkgfunc)
        psf_results = psfphot(scienceimage, init_params=init_params)

        if wcs is not None:
            init_radec = wcs.pixel_to_world(psf_results['x_init'], psf_results['y_init'])
            radec = wcs.pixel_to_world(psf_results['x_fit'], psf_results['y_fit'])

            psf_results['ra_init'] = init_radec.ra.value
            psf_results['dec_init'] = init_radec.dec.value

            psf_results['ra'] = radec.ra.value
            psf_results['dec'] = radec.dec.value

        return psf_results

    except NonFiniteValueError:
        print('fit_shape overlaps with edge of image, and therefore encloses NaNs! Photometry cancelled.')

