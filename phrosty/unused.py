# ruff: noqa

# These are things that may be out of date and may no longer work as is,
# but we may need to use something like them at some point in the
# future.

# If you want to use this, move them somewhere better and write tests.
# If we decide we really don't need them, delete them from this file.

# def build_psf(scienceimage,coords,wcs,ap_r=9,plot_epsf=False,
#             saturation=0.9e5, noise=1e4, method='subpixel',subpixels=5, 
#             fwhm=3.0, oversampling=3, maxiters=3, forced_photometry=True,
#             exclude_duplicates=False):

#     """
#     Build PSF from field stars. 
#     """

#     mean, median, stddev = sigma_clipped_stats(scienceimage)
#     daofind = DAOStarFinder(fwhm=fwhm,threshold = 5.*(stddev))
#     ap_results = ap_phot(scienceimage,coords,
#                         ap_r=ap_r, method=method, 
#                         subpixels=subpixels, merge_tables=True)

#     psfstars = Table({'x': ap_results['xcentroid'], 'y': ap_results['ycentroid'],
#                             'flux': ap_results['flux'], 'max': ap_results['max']})
#     # NOTE: Need to make star and galaxy separation work in order to make this work. 
#     print('Number of stars before saturation and flux cuts for ePSF:', len(psfstars))
#     psfstars = psfstars[psfstars['max'] < saturation]
#     psfstars = psfstars[psfstars['flux'] > noise]
#     print('Number of stars after saturation and flux cuts for ePSF:', len(psfstars))

#     stampsize=25
#     nddata = NDData(data=scienceimage)
#     extracted_stars = extract_stars(nddata, psfstars, size=stampsize)

#     if exclude_duplicates:
#         # Get rid of stamps with more than one source.
#         exclude_coords = []
#         for i in range(len(extracted_stars)):
#             try:
#                 stampsources = daofind(extracted_stars[i] - median)
#                 if len(stampsources) > 1 or len(stampsources) < 1:
#                     exclude_coords.append(extracted_stars.center_flat[i])
#             except:
#                 pass

#         exclude_rows = []
#         for c in exclude_coords:
#             exclude_rows.append(psfstars[psfstars['x'] == c[0]])

#         new_psfstars_rows = [x for x in psfstars if x not in exclude_rows]
#         new_psfstars = Table(rows=new_psfstars_rows, names=psfstars.colnames)
#         extracted_stars = extract_stars(nddata, new_psfstars, size=stampsize)

#     # Build ePSF.
#     print('Number of extracted stars for ePSF:', len(extracted_stars))
#     epsf_builder = EPSFBuilder(oversampling=oversampling, maxiters=maxiters)
#     psf_func, fitted_stars = epsf_builder(extracted_stars)
    
#     if plot_epsf:
#         norm = simple_norm(psf_func.data, 'log', percent=99.0)
#         plt.imshow(psf_func.data, norm=norm, origin='lower', cmap='Greys')
#         plt.colorbar()
#         plt.title('ePSF')
#         plt.show()

#     if forced_photometry:
#         psf_func.x_0.fixed = True
#         psf_func.y_0.fixed = True

#     return psf_func

# def crossmatch(pi,ti,seplimit=0.1):
#     """ Cross-match the truth files from each image (TI) to the corresponding photometry
#     file from that image (PI).

#     :param ti: Astropy table from a truth file ('Roman_TDS_index_{band}_{pointing}_{sca}.txt')
#     :type ti: Astropy table
#     :param pi: Astropy table from a photometry file generated from the image with the same band,
#                 pointing, and sca as ti. 
#     :type pi: Astropy table

#     :return: Joined truth catalog and measured photometry catalog, 
#             so that the measured objects are correctly crossmatched
#             to their corresponding rows in the truth catalog. 
#     :rtype: astropy.table 
#     """

#     if 'ra_truth' not in ti.colnames:
#         ti['ra'].name = 'ra_truth'
    
#     if 'dec_truth' not in ti.colnames:
#         ti['dec'].name = 'dec_truth'

#     if 'flux_truth' not in ti.colnames:
#         ti['flux'].name = 'flux_truth'

#     if 'mag_truth' not in ti.colnames:
#         ti['mag'].name = 'mag_truth'

#     tc = SkyCoord(ra=ti['ra_truth']*u.degree, dec=ti['dec_truth']*u.degree)
#     pc = SkyCoord(ra=pi['ra']*u.degree, dec=pi['dec']*u.degree)
#     ti_idx, pi_idx, angsep, dist3d = search_around_sky(tc,pc,seplimit=seplimit*(u.arcsec))

#     ti_reduced = ti[ti_idx]
#     pi_reduced = pi[pi_idx]

#     ti_pi_reduced = hstack([ti_reduced,pi_reduced], join_type='exact')
#     ti_x_pi = join(ti,ti_pi_reduced,join_type='outer')
    
#     return ti_x_pi

# def convert_flux_to_mag(ti_x_pi, band, zpt=False):
#     """Convert all fluxes to magnitudes from the crossmatched table. 

#     :param ti_x_pi: Astropy table, directly output from photometry.crossmatch.
#     :type ti_x_pi: astropy.table.Table
#     :param band: Roman filter.
#     :type band: str
#     :param zpt: Set to True if you want to zeropoint the fit flux values from PSF photometry
#                 to the truth catalog, as well as apply the galsim zeropoint to the truth magnitudes. Defaults to False.
#     :type zpt: bool, optional
#     :return: Input ti_x_pi with additional columns. 
#     :rtype: astropy.table.Table

#     """

#     exptime = {'F184': 901.175,
#                'J129': 302.275,
#                'H158': 302.275,
#                'K213': 901.175,
#                'R062': 161.025,
#                'Y106': 302.275,
#                'Z087': 101.7}

#     ti_x_pi['mag_fit'] = -2.5*np.log10(ti_x_pi['flux_fit'])
#     ti_x_pi['mag_err'] = np.sqrt((1.09/ti_x_pi['flux_fit'])**2*ti_x_pi['flux_err']**2)

#     if zpt:
#         ti_x_pi['zpt'] = np.zeros(len(ti_x_pi))
#         area_eff = roman.collecting_area
#         galsim_zp = roman.getBandpasses()[band].zeropoint

#         ti_x_pi['truth_mag'] = -2.5*np.log10(ti_x_pi['flux_truth']) + 2.5*np.log10(exptime[band]*area_eff) + galsim_zp
        
#         # This is going to be slow. Should parallelize. 
#         # First of all, looping through the entire crossmatched object list is inefficient and can be
#         # parallelized. Second of all, get_object_instances has a slow part in it that should also
#         # be parallelized. 
#         for i, row in enumerate(ti_x_pi):
#             print(row['object_id'], row['ra_truth'], row['dec_truth'])
#             objtab = get_object_instances(row['object_id'], row['ra_truth'], row['dec_truth'], bands=band)
#             objdata = get_object_data(row['object_id'], objtab)
#             objdata['mag_fit'] = -2.5*np.log10(objdata['flux_fit'])
#             mean_mag = np.nanmedian(objdata['mag_fit'])
#             zpt = row['mag_fit'] - mean_mag
#             # zpt = np.unique(objdata['mag_truth'] - mean_mag)
#             ti_x_pi['zpt'][i] = zpt
#             ti_x_pi['mag_fit'] += zpt 

#     return ti_x_pi

#         #     if zpt == 'truth':
# #         ap_zpt_mask = np.logical_and(results_table[f'{self.band}_ap_mag']>-11, results_table[f'{self.band}_ap_mag']<-9)
# #         psf_zpt_mask = np.logical_and(results_table[f'{self.band}_psf_mag']>-11, results_table[f'{self.band}_psf_mag']<-9)

# #         truthmag = truth_table[self.band][self.footprint_mask]
# #         results_table[f'{self.band}_truth'] = truthmag
# #         ap_zpt = np.median(results_table[f'{self.band}_ap_mag'][ap_zpt_mask] - truthmag[ap_zpt_mask])
# #         psf_zpt = np.median(results_table[f'{self.band}_psf_mag'][psf_zpt_mask] - truthmag[psf_zpt_mask])
        
# #         results_table[f'{self.band}_ap_mag'] -= ap_zpt
# #         results_table[f'{self.band}_psf_mag'] -= psf_zpt


# FUNCTION IS OBSOLETE.
#
# Delete it once we feel safe to do so.
#
# def _coord_transf(ra,dec):
#     """
#     Helper function for _sca_check and get_object_instances.get p
#     Inputs must be in radians. 

#     :return: Transformed x, y, z coordinates of given RA, dec. 
#     :rtype: tuple  

#     """
#     x = np.cos(dec) * np.cos(ra)
#     y = np.cos(dec) * np.sin(ra)
#     z = np.sin(dec)

#     return x,y,z

# FUNCTION IS OBSOLETE.
#
# Delete it once we feel safe to do so.
#
# def _distance(x0,x1,y0,y1,z0,z1):
#     """Distance formula. Helper function for _sca_check and get_object_instances. 

#     :param x0: _description_
#     :type x0: _type_
#     :param x1: _description_
#     :type x1: _type_
#     :param y0: _description_
#     :type y0: _type_
#     :param y1: _description_
#     :type y1: _type_
#     :param z0: _description_
#     :type z0: _type_
#     :param z1: _description_
#     :type z1: _type_
#     :return: _description_
#     :rtype: _type_
#     """    
#     return np.sqrt((x0 - x1)**2 + (y0 - y1)**2 + (z0 - z1)**2)

# FUNCTION IS OBSOLETE.
#
# Delete it once we feel safe to do so.
#
# def _sca_check(sca_ra, sca_dec, oid_x, oid_y, oid_z):

#     """
#     This is a helper function for get_object_instances and is 
#     not meant to be called directly. 

#     :return: Indices of which specific SCAs contain a given object.
#             Returns indices in columns, where the first column is 
#             SCA-1 and the second is the pointing that these SCAs
#             belong to. 
#     :rtype: np.array
    
#     """
#     # Each SCA is 0.11 arcsec/px and 4088 x 4088 px. 
#     # So, at most, the object will be the distance between the
#     # center coordinate and the corner of the SCA away from the 
#     # center. 

#     sca_ra = (sca_ra * u.deg).to(u.rad).value
#     sca_dec = (sca_dec * u.deg).to(u.rad).value

#     sca_halfsize = ((0.11 * 4088 / 2.) * u.arcsec).to(u.rad).value
#     d_req = np.sqrt(2)*sca_halfsize

#     sca_x, sca_y, sca_z = _coord_transf(sca_ra, sca_dec)

#     d_actual = _distance(sca_x, oid_x, sca_y, oid_y, sca_z, oid_z)

#     idx = np.stack(np.where(d_actual < d_req)).T

#     return idx

# FUNCTION IS OBSOLETE.
#
# Delete it once we feel safe to do so.
#
# def _obj_in(oid,df):
#     """Helper function for get_object_instances. Tells you if
#         a row in a table has the specified object ID.     

#     :param oid: Unique object ID. 
#     :type oid: int
#     :param df: _description_
#     :type df: _type_
#     :return: _description_
#     :rtype: boolean
#     """    
#     if int(oid) in df['object_id'].astype(int):
#         return True
#     else:
#         return False    

# FUNCTION IS OBSOLETE.
#
# Delete it once we feel safe to do so.
#
# def get_object_instances(ra,dec,oid=None,bands=get_roman_bands(),
#                         pointings=np.arange(0,57365,1),
#                         mjd_start=-np.inf,mjd_end=np.inf,
#                         obseq_path=obseq_path,
#                         obseq_radec_path=obseq_radec_path):

#     """
#     Retrieves all images that a unique object or set of coordinates is in. There are three steps to this, because
#     I think it will make the code run faster:
#     1. First cut (coarse): Check object's proximity to boresight coordinates for all pointings.
#     2. Second cut (fine): Of the pointings that passed the first cut, check the proximity of the center
#        of each individual SCA to the object's coordinates. If an object ID is not provided, this is the 
#        final step. 
#     3. Of the SCAs that passed the second cut, open each truth file and check if the object ID is in it. 

#     RA/dec arguments should be in degrees, as they are in the obseq file. 

#     :param ra: 
#     :type ra: float
#     :param dec:
#     :type dec: float
#     :param oid: If None, return table containing SCAs that contain the provided RA, dec. If an object ID is
#                 provided in this field, this function checks each truth file associated with each image 
#                 listed in the table it assembles to ensure that the particular object is present in those 
#                 images. WARNING: This is slow. 
#     :type oid: int or None, optional
#     :param band: Filters to include in search. Default ['F184', 'H158', 'J129', 'K213', 'R062', 'Y106', 'Z087'].
#     :type band: list or str, optional
#     :param pointings: Limit search to particular pointings.  
#     :type pointings: list or np.ndarray, optional
#     :param mjd_start: Start MJD to include in search. 
#     :type mjd_start: float, optional
#     :param mjd_end: End MJD to include in search. 
#     :type mjd_end: float, optional
#     :return: Astropy table with columns filter, pointing, SCA identifying images that contain the input RA 
#             and dec. If and object ID is provided in argument oid, this function checks each truth file associated 
#             with each image listed in the table it assembles to ensure that the particular object is present in those 
#             images. WARNING: This is slow. Note that if an object ID is not provided, the final list returned will
#             contain some images where the provided RA and dec are slightly outside its bounds. This is because for
#             speed, it takes the center of the SCA from the obseq files, and circumscribes a circle around the SCA for
#             the search zone. 
#     :rtype: astropy.table.Table
#     """
    
#     with fits.open(obseq_path) as osp:
#         obseq_orig = Table(osp[1].data)

#     with fits.open(obseq_radec_path) as osradecp:
#         obseq_radec_orig = Table(osradecp[1].data)

#     pointing_idx = np.array(pointings)
#     obseq = obseq_orig[pointing_idx]
#     obseq_radec = obseq_radec_orig[pointing_idx]

#     band_idx = np.where(np.in1d(obseq['filter'],bands))[0]
#     obseq = obseq[band_idx]
#     obseq_radec = obseq_radec[band_idx]
#     pointing_idx = pointing_idx[band_idx]

#     mjd_idx = np.where((obseq['date'] > mjd_start) & (obseq['date'] < mjd_end))[0]
#     obseq = obseq[mjd_idx]
#     obseq_radec = obseq_radec[mjd_idx]
#     pointing_idx = pointing_idx[mjd_idx]

#     ra_oid = (ra * u.deg).to(u.rad).value
#     dec_oid = (dec * u.deg).to(u.rad).value

#     x_oid, y_oid, z_oid = _coord_transf(ra_oid, dec_oid)

#     ra_obseq = (obseq['ra'] * u.deg).to(u.rad).value 
#     dec_obseq = (obseq['dec'] * u.deg).to(u.rad).value

#     x, y, z = _coord_transf(np.array(ra_obseq), np.array(dec_obseq))

#     d_actual = _distance(x, x_oid, y, y_oid, z, z_oid)
#     d_req = np.sin(0.009/2.) # Why? Because in roman_imsim.telescope.near_pointing(), 
#                              # this is the required distance. See where it says
#                              # self.sbore2 = np.sin(max_rad_from_boresight/2.), and
#                              # earlier, defines max_rad_from_boresight = 0.009 by default.
#                              # This part of my code is also basically a copy of near_pointing. 

#     distance_idx = np.where(d_actual/2. <= d_req)[0]
#     idx_firstcut = pointing_idx[distance_idx]

#     sca_tab = obseq_radec_orig[idx_firstcut]
#     sca_tab_ra = np.array(sca_tab['ra'])
#     sca_tab_dec = np.array(sca_tab['dec'])
#     sca_tab_filter = np.array(sca_tab['filter'])

#     sca_rd = np.array([sca_tab_ra,sca_tab_dec]).reshape((2,len(sca_tab),18)).T
#     sca_check = _sca_check(sca_rd[:,:,0], sca_rd[:,:,1],x_oid,y_oid,z_oid)
#     pointing_idx = sca_check[:,1]
#     sca_idx = sca_check[:,0]

#     secondcut_tab = Table([sca_tab_filter[pointing_idx], idx_firstcut[pointing_idx], sca_idx+1], names=('filter','pointing','sca'))

#     # Now, make sure the provided RA, dec are actually in these images. 
#     # This is slow because it opens each file for its WCS. 
#     inimg = list(map(radec_isin, [ra]*len(secondcut_tab), [dec]*len(secondcut_tab), \
#                     [None]*len(secondcut_tab), secondcut_tab['filter'], secondcut_tab['pointing'], secondcut_tab['sca']))

#     thirdcut_tab = secondcut_tab[inimg]

#     if oid is not None:
#         # This part of the code is really slow because I'm opening files. 
#         # Want to parallelize in future to speed up. 
#         final_idx = []
#         for i, row in enumerate(thirdcut_tab):
#             df = read_truth_txt(band=row['filter'],pointing=row['pointing'],sca=row['sca'])
#             if _obj_in(oid, df):
#                 final_idx.append(i)
#             del df

#         tab = thirdcut_tab[final_idx]

#     else:
#         tab = thirdcut_tab

#     return tab

# FUNCTION IS OBSOLETE.
#
# Delete it once we feel safe to do so.
#
# def get_object_data(oid, metadata,
#                     colnames=['object_id','ra','dec','mag_truth','flux_truth','flux_fit','flux_err'],
#                     crossmatch_dir='/hpc/group/cosmology/lna18/roman_sim_imgs/Roman_Rubin_Sims_2024/preview/crossmatched_truth'):

#     # NOTE: THESE PATHS HAVE NOT BEEN GENERALIZED. 

#     """Retrieves all information from crossmatched photmetry files about a particular object, specified with its
#     unique object ID. 

#     :param oid: Object ID. 
#     :type oid: int
#     :param metadata: Output from get_object_instance. Table with filter, pointing, and SCA numbers for this
#                     function to search through. 
#     :type metadata: astropy.table
#     :param colnames: The columns that this function should retrieve from the crossmatched photometry files.
#     :type colnames: list, optional
#     :return: Astropy table with all instances of argument oid found in the crossmatched photometry files
#             that are present in the metadata table. 
#     :rtype: astropy.table.Table

#     """
    
#     object_tab = Table(names=colnames)
#     for row in metadata:
#         band = row['filter']
#         p = row['pointing']
#         sca = row['sca']

#         filepath = pa.join(crossmatch_dir,f'{band}/{p}/Roman_TDS_xmatch_{band}_{p}_{sca}.txt')
#         phot = Table.read(filepath, format='csv')
#         object_row = phot[phot['object_id'] == int(oid)]
#         object_row_reduced = object_row[colnames]
#         object_tab = vstack([object_tab, object_row_reduced], join_type='exact')

#     return object_tab
