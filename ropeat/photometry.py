import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, hstack, MaskedColumn, join
from astropy.visualization import ZScaleInterval, simple_norm
import astropy.units as u
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats
from photutils.background import LocalBackground, MMMBackground, Background2D
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFBuilder, extract_stars, PSFPhotometry

roman_bands = ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'W146', 'K213']

def ap_phot(scienceimage,coords,wcs,
            ap_r=3, bkg_estimator=MMMBackground(), box_size=(50,50),
            filter_size=(3,3), method='subpixel',subpixels=5,
            merge_tables=True):
    """_summary_

    Args:
        scienceimage (array-like): Array containing science image.
        coords (astropy.table.Table): Table with columns 'x' and 'y' representing pixel
                                    coordinates for the detected sources. Table can contain
                                    more than just 'x' and 'y', and this may be useful if you 
                                    use merge_tables=True. 
        ap_r (int, optional): Aperture radius to use for aperture photometry. Defaults to 3.
        bkg_estimator (_type_, optional): _description_. Defaults to MMMBackground().
        box_size (tuple, optional): _description_. Defaults to (50,50).
        filter_size (tuple, optional): _description_. Defaults to (3,3).
        method (str, optional): _description_. Defaults to 'subpixel'.
        subpixels (int, optional): _description_. Defaults to 5.
        merge_tables (bool, optional): If true, output is merged coords and results from aperture
                                    photometry. Defaults to True.

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

def convert_flux_to_mag(fluxtab, zpt=False, truth=None):
    """
    Input the astropy table from psf_results.

    """

    fluxtab['mag_fit'] = -2.5*np.log10(fluxtab['flux_fit'])
    fluxtab['mag_err'] = np.sqrt((1.09/fluxtab['flux_fit'])**2*fluxtab['flux_err']**2)

    if zpt:
        if truth is None:
            raise ValueError('You need to provide a truth catalog if you want to zero point the magnitudes.')

        else:
            print('YOU HAVE TO CROSSMATCH BEFORE YOU ZERO POINT')
    return fluxtab

        #     if zpt == 'truth':
#         ap_zpt_mask = np.logical_and(results_table[f'{self.band}_ap_mag']>-11, results_table[f'{self.band}_ap_mag']<-9)
#         psf_zpt_mask = np.logical_and(results_table[f'{self.band}_psf_mag']>-11, results_table[f'{self.band}_psf_mag']<-9)

#         truthmag = truth_table[self.band][self.footprint_mask]
#         results_table[f'{self.band}_truth'] = truthmag
#         ap_zpt = np.median(results_table[f'{self.band}_ap_mag'][ap_zpt_mask] - truthmag[ap_zpt_mask])
#         psf_zpt = np.median(results_table[f'{self.band}_psf_mag'][psf_zpt_mask] - truthmag[psf_zpt_mask])
        
#         results_table[f'{self.band}_ap_mag'] -= ap_zpt
#         results_table[f'{self.band}_psf_mag'] -= psf_zpt

    return fluxtab

def crossmatch(pi,ti,seplimit=0.1):
    """ Cross-match the truth files from each image (TI) to the corresponding photometry
    file from that image (PI).

    :param ti: Astropy table from a truth file ('Roman_TDS_index_{band}_{pointing}_{sca}.txt')
    :type ti: Astropy table
    :param pi: Astropy table from a photometry file generated from the image with the same band,
                pointing, and sca as ti. 
    :type pi: Astropy table

    Returns:
        _type_: _description_ 
    """

    tc = SkyCoord(ra=ti['ra_truth']*u.degree, dec=ti['dec']*u.degree)
    pc = SkyCoord(ra=pi['ra']*u.degree, dec=pi['dec']*u.degree)
    ti_idx, pi_idx, angsep, dist3d = search_around_sky(tc,pc,seplimit=seplimit(u.arcsec))

    # ti_reduced = ti[ti_idx]
    pi_reduced = pi[pi_idx]

    ti[[pi_reduced.colnames]][ti_idx] = pi_reduced

    tixpi = join(ti,pi,join_type='outer')
    return ti


# def crossmatch_truth(truth_filepath,results_filepaths,savename,overwrite=True,seplimit=0.1,psf=True,verbose=True,temp_file_path=None):
#     """
#     This will handle a list of science files with mixed bands! Hopefully!

#     :param truth_filepath: Filepath to truth file.
#     :type truth_filepath: str
#     :param results_filepaths: List of filepaths to csv files from scienceimg.save_results().
#     :type results_filepaths: list
#     :param seplimit: Angular separation limit to find matches, defaults to 0.1
#     :type seplimit: float, optional
#     :param psf: Include PSF photometry results in the matching table, defaults to True
#     :type psf: bool, optional
#     :return: _description_
#     :rtype: _type_
#     """    
#     prefixes = ['ap']
#     if psf:
#         prefixes.append('psf')
        
#     suffixes = ['_flux', '_flux_err', '_mag', '_mag_err']
#     match_vals = []
#     all_suffixes = []
    
#     for p in prefixes:
#         for s in suffixes:
#             if p == 'ap' and 'err' in s:
#                 pass
#             else:
#                 all_suffixes.append('_'+p+s)
#             for b in roman_bands:
#                 if p == 'ap' and 'err' in s:
#                     pass
#                 else:
#                     match_vals.append(b+'_'+p+s)
                    
#                 if f'{b}_max' not in match_vals:
#                     match_vals.append(f'{b}_max')

#     match_vals.append('ra_all')
#     match_vals.append('dec_all')
#     match_vals.append('x_truth')
#     match_vals.append('y_truth')
#     match_vals.append('x_fit')
#     match_vals.append('y_fit')
#     match_vals.append('localbkg_all')
#     match_vals.append('pointing_all')
#     match_vals.append('sca_all')
#     match_vals.append('psfphot_flags_all')
    
#     tr = truth(truth_filepath)
#     tr_tab = tr.table()
#     if temp_file_path is None:
#         temp_file_name = 'tempfile_DELETEME.fits'
#     elif isinstance(temp_file_path, str) and '.fits' in temp_file_path:
#         temp_file_name = temp_file_path
#     else:
#         print('Something is wrong with your temporary file path.')
#         print('Maybe it doesnt end in .fits? Or is not a string?')

#     for col in match_vals:
#         tr_tab.add_column('empty', name=col)
        
#     for col in tr_tab.colnames:
#         if col in roman_bands:
#             tr_tab.rename_column(col, f'{col}_truth')
            
#     tr_tab.rename_column('ra','ra_truth')
#     tr_tab.rename_column('dec','dec_truth')
    
#     tr_tab.write(temp_file_name, format='fits', overwrite=True)

#     if verbose:
#         print(f'We need to get through matching {len(results_filepaths)} files.')
#         nfiles = 0
#         nfail = 0
#     for i, file in enumerate(results_filepaths):
#         if os.path.exists(file):
#             if verbose:
#                 print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
#                 print(file)
#             with fits.open(temp_file_name) as f:
#                 tr_tab = Table(f[1].data)
#                 check = Table.read(file, format='csv')
#                 band = check['band'][0]
#                 if verbose:
#                     print('Succesfully opened file. Next, crossmatch and merge tables.')
#                 check_coords = SkyCoord(ra=check['ra']*u.degree, dec=check['dec']*u.degree)
#                 tr_coords = SkyCoord(ra=tr_tab['ra_truth']*u.degree, dec=tr_tab['dec_truth']*u.degree) 
#                 check_idx, tr_idx, angsep, dist3d = search_around_sky(check_coords,tr_coords,seplimit=seplimit*u.arcsec)

#                 tr_tab = tr_tab.to_pandas()
#                 appendvals = lambda x,y : str(x) + ',' + str(y)
                
#                 # Collect all magnitudes into a string into one column. 
#                 for s in all_suffixes:
#                     tlist = list(tr_tab[f'{band}{s}'])
#                     tlist_reduced = list(np.array(tlist)[tr_idx])
#                     clist = list(check[f'{band}{s}'])
#                     clist_reduced = list(np.array(clist)[check_idx])
#                     strcol = list(map(appendvals,tlist_reduced,clist_reduced))
#                     strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
#                     tr_tab.loc[tr_idx, f'{band}{s}'] = strcol
                    
#                 # Collect RA/dec into a string into one column. 
#                 for c in ['ra','dec']:
#                     tlist = list(tr_tab[f'{c}_all'])rm 
#                     tlist_reduced = list(np.array(tlist)[tr_idx])
#                     clist = list(check[c])
#                     clist_reduced = list(np.array(clist)[check_idx])
#                     strcol = list(map(appendvals,tlist_reduced,clist_reduced))
#                     strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
#                     tr_tab.loc[tr_idx, f'{c}_all'] = strcol
                    
#                 # Collect all x/y positions into a string into one column.
#                 for c in ['x','y']:
#                     for a in ['fit','truth']:
#                         tlist = list(tr_tab[f'{c}_{a}'])
#                         tlist_reduced = list(np.array(tlist)[tr_idx])
#                         clist = list(check[f'{c}_{a}'])
#                         clist_reduced = list(np.array(clist)[check_idx])
#                         strcol = list(map(appendvals,tlist_reduced,clist_reduced))
#                         strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
#                         tr_tab.loc[tr_idx, f'{c}_{a}'] = strcol
                    
#                 # Collect the max. flux pixel value in an aperture into one column. 
#                 tlist = list(tr_tab[f'{band}_max'])
#                 tlist_reduced = list(np.array(tlist)[tr_idx])
#                 clist = list(check['max'])
#                 clist_reduced = list(np.array(clist)[check_idx])
#                 strcol = list(map(appendvals,tlist_reduced,clist_reduced))
#                 strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
#                 tr_tab.loc[tr_idx, f'{band}_max'] = strcol
                    
#                 # Collect all SCA IDs into a string in one column.
#                 tlist = list(tr_tab['sca_all'])
#                 tlist_reduced = list(np.array(tlist)[tr_idx])
#                 clist = list(check['chip'])
#                 clist_reduced = list(np.array(clist)[check_idx])
#                 strcol = list(map(appendvals,tlist_reduced,clist_reduced))
#                 strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
#                 tr_tab.loc[tr_idx, 'sca_all'] = strcol
                
#                 # Collect all pointings into a string in one column.
#                 tlist = list(tr_tab['pointing_all'])
#                 tlist_reduced = list(np.array(tlist)[tr_idx])
#                 clist = list(check['pointing'])
#                 clist_reduced = list(np.array(clist)[check_idx])
#                 strcol = list(map(appendvals,tlist_reduced,clist_reduced))
#                 strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
#                 tr_tab.loc[tr_idx, 'pointing_all'] = strcol
                
#                 # Collect all psfphot flags into a string in one column.
#                 tlist = list(tr_tab['psfphot_flags_all'])
#                 tlist_reduced = list(np.array(tlist)[tr_idx])
#                 clist = list(check['psfphot_flags'])
#                 clist_reduced = list(np.array(clist)[check_idx])
#                 strcol = list(map(appendvals,tlist_reduced,clist_reduced))
#                 strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
#                 tr_tab.loc[tr_idx, 'psfphot_flags_all'] = strcol
                
#                 # Collect all localbkg values into a string in one column.
#                 tlist = list(tr_tab['localbkg_all'])
#                 tlist_reduced = list(np.array(tlist)[tr_idx])
#                 clist = list(check['localbkg'])
#                 clist_reduced = list(np.array(clist)[check_idx])
#                 strcol = list(map(appendvals,tlist_reduced,clist_reduced))
#                 strcol = [strcol[i][6:] if strcol[i][0] == 'e' else strcol[i] for i in range(len(strcol))]
#                 tr_tab.loc[tr_idx, 'localbkg_all'] = strcol

#             tr_tab = Table.from_pandas(tr_tab)
#             tr_tab.write(temp_file_name, format='fits', overwrite=True)
#             if verbose:
#                 print(f'Wrote the crossmatched data to the main catalog in a temporary file, {temp_file_name}.')
#                 nfiles += 1
#                 print(f'Made it through crossmatching {nfiles}/{len(results_filepaths)} files.')

#         else:
#             if verbose:
#                 print(f'Oops! {file} does not seem to exist.')
#                 nfail += 1
#                 print(f'In total, {nfail}/{len(results_filepaths)} filepaths have failed.')
        
#     # Drop empty columns. 
#     for col in [c for c in tr_tab.colnames if c not in ['object_id','config']]:
#         try:
#             tr_tab[col] = MaskedColumn(data=tr_tab[col].value, dtype=np.float64)
#             tr_tab[col].mask = np.isinf(tr_tab[col].value)
#             if all(tr_tab[col].mask):
#                 tr_tab.remove_column(col)

#         except:
#             if all(tr_tab[col] == 'empty'):
#                 tr_tab.remove_column(col)
            
#     # Write table to file. 
#     tr_tab.write(savename, format='csv', overwrite=overwrite)
#     if verbose:
#         print('Final crossmatched file is written.')
#         print('Finally, deleting the temporary file.')
#     os.remove(temp_file_name)