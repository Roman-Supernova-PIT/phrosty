import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, vstack
from astropy import units as u

def read_truth_txt(band,pointing,sca):
    """
    Reads in the txt versions of the truth files as convenient astropy tables. 

    """
    truthpath = f'/cwork/mat90/RomanDESC_sims_2024/RomanTDS/truth/{band}/{pointing}/Roman_TDS_index_{band}_{pointing}_{sca}.txt'
    truth_colnames = ['object_id', 'ra', 'dec', 'x', 'y', 'realized_flux', 'flux', 'mag', 'obj_type']
    truth_pd = pd.read_csv(truthpath, comment='#', skipinitialspace=True, sep=' ', names=truth_colnames)
    truth = Table.from_pandas(truth_pd)

    return truth

def get_corners(band,pointing,sca):
    imgpath = f'/cwork/mat90/RomanDESC_sims_2024/RomanTDS/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz'
    with fits.open(imgpath) as hdu:
        wcs = WCS(hdu[1].header)
    corners = [[0,0],[0,4088],[4088,0],[4088,4088]]
    wcs_corners = wcs.pixel_to_world_values(corners)

    return wcs_corners

def get_mjd(pointing):
    obseq_path = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/Roman_TDS_obseq_11_6_23.fits'
    with fits.open(obseq_path) as obs:
        obseq = Table(obs[1].data)
    mjd = obseq['date'][pointing]

    return mjd

def _coord_transf(ra,dec):
    """
    Helper function for sca_check and get_object_instances.
    Inputs must be in radians. 

    """
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    return x,y,z

def _distance(x0,x1,y0,y1,z0,z1):
    return np.sqrt((x0 - x1)**2 + (y0 - y1)**2 + (z0 - z1)**2)

def _sca_check(sca_ra, sca_dec, oid_x, oid_y, oid_z):

    """
    This is a helper function for get_object_instances and is not meant to be called 
    directly. 
    
    """
    # Each SCA is 0.11 arcsec/px and 4088 x 4088 px. 
    # So, at most, the object will be the distance between the
    # center coordinate and the corner of the SCA away from the 
    # center. 

    sca_ra = (sca_ra * u.deg).to(u.rad).value
    sca_dec = (sca_dec * u.deg).to(u.rad).value

    sca_halfsize = ((0.11 * 4088 / 2.) * u.arcsec).to(u.rad).value
    d_req = np.sqrt(2)*sca_halfsize

    sca_x, sca_y, sca_z = _coord_transf(sca_ra, sca_dec)

    d_actual = _distance(sca_x, oid_x, sca_y, oid_y, sca_z, oid_z)

    idx = np.stack(np.where(d_actual < d_req)).T
    # Note: returns indices in columns, where the first column is the
    # sca-1 and the second is pointing
    return idx

def _obj_in(oid,df):
    if int(oid) in df['object_id'].astype(int):
        return True
    else:
        return False

def get_object_instances(oid,ra,dec,
                        mjd_start=None,mjd_end=None):

    """
    Retrieves all images that a unique object is in. There are three steps to this, because
    I think it will make the code run faster:
    1. First cut (coarse): Check object's proximity to boresight coordinates for all pointings.
    2. Second cut (fine): Of the pointings that passed the first cut, check the proximity of the center
       of each individual SCA to the object's coordinates.
    3. Of the SCAs that passed the second cut, open each truth file and check if the object ID is in it. 


    RA/dec arguments should be in degrees, as they are in the obseq file. 

    """
    
    obseq_path = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/Roman_TDS_obseq_11_6_23.fits'
    with fits.open(obseq_path) as osp:
        obseq = Table(osp[1].data)

    obseq_radec_path = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/Roman_TDS_obseq_11_6_23_radec.fits'
    with fits.open(obseq_radec_path) as osradecp:
        obseq_radec = Table(osradecp[1].data)

    if mjd_start is not None and mjd_end is not None:
        mjd_idx = np.where((obseq['date'] > mjd_start) & (obseq['date'] < mjd_end))[0]
    elif mjd_start is not None and mjd_end is None:
        mjd_idx = np.where(obseq['date'] > mjd_start)[0]
    elif mjd_end is not None and mjd_start is None:
        mjd_idx = np.where(obseq['date'] < mjd_end)[0]
    elif mjd_start is None and mjd_end is None:
        mjd_idx = np.arange(0,len(obseq),1)

    obseq = obseq[mjd_idx]
    obseq_radec = obseq_radec[mjd_idx]

    ra_oid = (ra * u.deg).to(u.rad).value
    dec_oid = (dec * u.deg).to(u.rad).value

    x_oid, y_oid, z_oid = _coord_transf(ra_oid, dec_oid)

    ra_obseq = (obseq['ra'] * u.deg).to(u.rad).value 
    dec_obseq = (obseq['dec'] * u.deg).to(u.rad).value

    x, y, z = _coord_transf(np.array(ra_obseq), np.array(dec_obseq))

    d_actual = _distance(x, x_oid, y, y_oid, z, z_oid)
    d_req = np.sin(0.009/2.) # Why? Because in roman_imsim.telescope.near_pointing(), 
                             # this is the required distance. See where it says
                             # self.sbore2 = np.sin(max_rad_from_boresight/2.), and
                             # earlier, defines max_rad_from_boresight = 0.009 by default.
                             # This part of my code is also basically a copy of near_pointing. 
    
    idx_pointing = np.where(d_actual/2. <= d_req)[0]
    idx_firstcut = mjd_idx[idx_pointing]

    sca_tab = obseq_radec[idx_firstcut]
    sca_tab_ra = np.array(sca_tab['ra'])
    sca_tab_dec = np.array(sca_tab['dec'])
    sca_tab_filter = np.array(sca_tab['filter'])

    sca_rd = np.array([sca_tab_ra,sca_tab_dec]).reshape((2,len(sca_tab),18)).T
    sca_check = _sca_check(sca_rd[:,:,0], sca_rd[:,:,1],x_oid,y_oid,z_oid)
    pointing_idx = sca_check[:,1]
    sca_idx = sca_check[:,0]

    secondcut_tab = Table([sca_tab_filter[pointing_idx], idx_firstcut[pointing_idx], sca_idx+1], names=('filter','pointing','sca'))

    # This part of the code is really slow because I'm opening files. 
    # Want to parallelize in future to speed up. 
    final_idx = []
    for i, row in enumerate(secondcut_tab):
        df = read_truth_txt(row['filter'],row['pointing'],row['sca'])
        if _obj_in(oid, df):
            final_idx.append(i)
        del df

    return secondcut_tab[final_idx]

def get_object_data(oid, metadata,
                    colnames=['object_id','ra','dec','mag_truth','flux_truth','flux_fit','flux_err']):

    """
    metadata is the output from get_object_instance. 
    """
    
    object_tab = Table(names=colnames)
    for row in metadata:
        band = row['filter']
        p = row['pointing']
        sca = row['sca']

        filepath = f'/hpc/group/cosmology/lna18/roman_sim_imgs/Roman_Rubin_Sims_2024/preview/crossmatched_truth/{band}/{p}/Roman_TDS_xmatch_{band}_{p}_{sca}.txt'
        phot = Table.read(filepath, format='csv')
        object_row = phot[phot['object_id'] == int(oid)]
        object_row_reduced = object_row[colnames]
        object_tab = vstack([object_tab, object_row_reduced], join_type='exact')

    return object_tab

def train_config(objtype):
    if objtype == 'galaxy':
        return 0
    elif objtype == 'star':
        return 1
    elif objtype == 'transient':
        return 2
    else: 
        return 3

def predict_config(objtype):
    if objtype == 0:
        return 'galaxy'
    elif objtype == 1:
        return 'star'
    elif objtype == 2:
        return 'transient'
    else:
        return 'other'

# def get_obj_type_from_ID(ID):
# THIS IS NO LONGER TRUE. CHANGE ID BOUNDS.
# But also may not be necessary because object type is now
# a column in the truth txt files. 
#     if ID > 9921000000000 and ID < 10779202101973:
#         return 'galaxy'
#     elif ID > 30328699913 and ID < 50963307358:
#         return 'star'
#     elif ID > 20000001 and ID < 120026449:
#         return 'transient'
#     else:
#         return 'unknown'
