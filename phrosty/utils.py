# IMPORTS Standard:
import os
import os.path as pa
import numpy as np
import pandas as pd
import warnings
from glob import glob

# IMPORTS Astro:
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import skycoord_to_pixel
from astropy.table import Table, vstack
from astropy import units as u

# CHANGE THIS TO FIT YOUR USE: 
rootdir = os.getenv('ROOTDIR', '/cwork/mat90/RomanDESC_sims_2024/RomanTDS')

# The FITSFixedWarning is consequenceless and it will get thrown every single time you deal with a WCS. 
warnings.simplefilter('ignore',category=FITSFixedWarning)

def _build_filepath(path,band,pointing,sca,filetype,rootdir=rootdir):
    """_summary_

    :param path: _description_
    :type path: str
    :param band: _description_
    :type band: str
    :param pointing: _description_
    :type pointing: str
    :param sca: _description_
    :type sca: str
    :param filetype: _description_
    :type filetype: str
    :raises ValueError: _description_
    :raises ValueError: _description_
    :raises ValueError: _description_
    :return: _description_
    :rtype: _type_
    """

    # First, what kind of file are we looking for? 
    if filetype not in ['image', 'truth', 'truthtxt']:
        raise ValueError(f'filetype must be in {filetype}.')
    elif filetype == 'image':
        subdir = 'images/simple_model'
        prefix = 'simple_model'
        extension = 'fits.gz'
    elif filetype == 'truth':
        subdir = 'images/truth'
        prefix = 'truth'
        extension = 'fits.gz'
    elif filetype == 'truthtxt':
        subdir = 'truth'
        prefix = 'index'
        extension = 'txt'

    # Did you already provide a path? 
    if path is not None:
        return path
    elif (band is None) or (pointing is None) or (sca is None):
        raise ValueError('You need to specify band, pointing, and sca if you do not provide a full filepath.')
    elif (band is not None) and (pointing is not None) and (sca is not None):
        path = pa.join(rootdir,subdir,band,str(pointing),f'Roman_TDS_{prefix}_{band}_{str(pointing)}_{str(sca)}.{extension}')
        return path 

    elif (path is None) and (band is None) and (pointing is None) and (sca is None):
        raise ValueError('You need to provide either the full image path, or the band, pointing, and SCA.')

def get_roman_bands():
    """
    :return: List of bands included in the Roman-DESC TDS simulations. 
    :rtype: list
    """    
    return ['F184', 'H158', 'J129', 'K213', 'R062', 'Y106', 'Z087']

def read_truth_txt(path=None,band=None,pointing=None,sca=None):
    """
    Reads in the txt versions of the truth files as convenient astropy tables. 

    :param truthpath: Path to txt catalog version of truth file. If you do not 
                        provide this, you need to specify the arguments
                        band, pointing, and sca.
    :type truthpath: str, optional
    :param band: Roman filter. If you do not provide this, you need to provide truthpath.
    :type band: str, optional
    :param pointing: Pointing ID. If you do not provide this, you need to provide truthpath.
    :type pointing: str, optional
    :param sca: SCA ID. If you do not provide this, you need to provide truthpath.
    :type sca: str, optional 
    :return: Astropy table with contents of the specified catalog txt file. 
    :rtype: astropy.table.Table

    """

    _truthpath = _build_filepath(path=path,band=band,pointing=pointing,sca=sca,filetype='truthtxt')
    truth_colnames = ['object_id', 'ra', 'dec', 'x', 'y', 'realized_flux', 'flux', 'mag', 'obj_type']
    truth_pd = pd.read_csv(_truthpath, comment='#', skipinitialspace=True, sep=' ', names=truth_colnames)
    truth = Table.from_pandas(truth_pd)

    return truth

def get_fitsobj(path=None,band=None,pointing=None,sca=None):
    imgpath = _build_filepath(path,band,pointing,sca,'image')
    return fits.open(imgpath)

def radec_isin(ra,dec,path=None,band=None,pointing=None,sca=None):
    _imgpath = _build_filepath(path,band,pointing,sca,'image')
    with fits.open(_imgpath) as hdu:
        wcs = WCS(hdu[1].header)
    worldcoords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    x, y = skycoord_to_pixel(worldcoords,wcs)
    pxradec = np.array([x,y])
    if np.logical_or(any(pxradec < 0), any(pxradec > 4088)): 
        res = False
    else:
        res = True

    return res

def get_corners(path=None,band=None,pointing=None,sca=None):
    """Retrieves the RA, dec of the corners of the specified SCA in degrees. 

    :param band: Roman filter.
    :type band: str
    :param pointing: Pointing ID.
    :type pointing: str
    :param sca: SCA ID.
    :type sca: str
    :return: Tuple containing four numpy arrays, each with the RA and dec of the corner
            of the specified image in degrees. 
    :rtype: tuple
    """    
    _imgpath = _build_filepath(path,band,pointing,sca,'image')
    with fits.open(_imgpath) as hdu:
        wcs = WCS(hdu[1].header)
    corners = [[0,0],[0,4088],[4088,0],[4088,4088]]
    wcs_corners = wcs.pixel_to_world_values(corners)

    return wcs_corners

def get_transient_radec(oid):
    """
    Retrieve RA, dec of a transient based on its object ID. 
    """
    oid = int(oid)
    filepath = '/cwork/mat90/RomanDESC_sims_2024/roman_rubin_cats_v1.1.2_faint/snana*.parquet'
    file_list = glob(filepath)
    for file in file_list:
        # Read the Parquet file
        df = pd.read_parquet(file)
        if len(df[df['id'] == oid]) != 0:
            ra = df[df['id'] == oid]['ra'].values[0]
            dec = df[df['id'] == oid]['dec'].values[0]
    return ra, dec

def get_transient_mjd(oid):
    """
    Retrieve start and end dates of a transient based on its object ID. 
    """
    filepath = '/cwork/mat90/RomanDESC_sims_2024/roman_rubin_cats_v1.1.2_faint/snana*.parquet'
    file_list = glob(filepath)
    for file in file_list:
        # Read the Parquet file
        df = pd.read_parquet(file)
        if len(df[df['id'] == oid]) != 0:
            start = df[df['id'] == oid]['start_mjd'].values[0]
            end = df[df['id'] == oid]['end_mjd'].values[0]
    return start, end 

def get_mjd_limits(obseq_path=f'{rootdir}/Roman_TDS_obseq_11_6_23.fits'): 
    """
    Retrive the earliest and latest MJD in the simulations.
    """
    
    with fits.open(obseq_path) as obs:
        obseq = Table(obs[1].data)

    start = min(obseq['date'])
    end = max(obseq['date'])

    return start, end

def get_radec_limits(obseq_path=f'{rootdir}/Roman_TDS_obseq_11_6_23.fits'):
    """
    Retrieve RA, dec limits (boresight coordinates).
    """
    with fits.open(obseq_path) as obs:
        obseq = Table(obs[1].data)

    ra_min = min(obseq['ra'])
    ra_max = max(obseq['ra'])

    dec_min = min(obseq['dec'])
    dec_max = max(obseq['dec'])

    return {'ra': [ra_min,ra_max], 'dec': [dec_min, dec_max]}


def get_mjd(pointing,obseq_path=f'{rootdir}/Roman_TDS_obseq_11_6_23.fits'):
    """
    Retrieve MJD of a given pointing. 

    :param pointing: Pointing ID. 
    :type pointing: int
    :param obseq_path: Path to obseq file Roman_TDS_obseq_11_6_23.fits.
    :type obseq_path: str, optional
    :return: MJD of specified pointing. 
    :rtype: float
    """    

    with fits.open(obseq_path) as obs:
        obseq = Table(obs[1].data)
    mjd = float(obseq['date'][int(pointing)])

    return mjd

def pointings_near_mjd(mjd,window=3,obseq_path=f'{rootdir}/Roman_TDS_obseq_11_6_23.fits'):
    """
    Retrieve pointings near given MJD.

    :param mjd: Central MJD to search around.
    :type mjd: float
    :param window: Number of days around central MJD to include in search. 
    :type window: float
    :param obseq_path: Path to obseq file Roman_TDS_obseq_11_6_23.fits.
    :type obseq_path: str, optional
    :return: Pointings within specified MJD range. 
    :rtype: list
    """ 

    with fits.open(obseq_path) as obs:
        obseq = Table(obs[1].data)

    pointings = np.where(np.logical_and(obseq['date'] < mjd + window, obseq['date'] > mjd - window))[0]
    return pointings 

def get_mjd_info(mjd_start=-np.inf,mjd_end=np.inf,return_inverse=False,obseq_path = f'{rootdir}/Roman_TDS_obseq_11_6_23.fits'):
    """
    Get all pointings and corresponding filters between two MJDs.
    Returns an astropy table with columns 'filter' and 'pointing'. 
    Does not return an 'sca' column because every sca belonging to a
    given pointing satisfies an MJD requirement. 

    :param mjd_start: Start MJD, defaults to -np.inf
    :type mjd_start: float, optional
    :param mjd_end: End MJD, defaults to np.inf
    :type mjd_end: float, optional
    :param return_inverse: If true, returns all pointings outside the MJD range specified instead of inside. 
    :type return_inverse: bool
    :param obseq_path: Path to obseq file Roman_TDS_obseq_11_6_23.fits.
    :type obseq_path: str, optional
    :return: Astropy table with pointing numbers and corresponding filters that satisfy the
            MJD requirements. 
    :rtype: astropy.table.Table
    """
    with fits.open(obseq_path) as obs:
        obseq = Table(obs[1].data)

    if not return_inverse:
        mjd_idx = np.where((obseq['date'] > float(mjd_start)) & (obseq['date'] < float(mjd_end)))[0]
    elif return_inverse:
        mjd_idx = np.where(~((obseq['date'] > float(mjd_start)) & (obseq['date'] < float(mjd_end))))[0]

    mjd_tab = Table([obseq['filter'][mjd_idx], mjd_idx], names=('filter','pointing'))

    return mjd_tab

def get_exptime(band=None):

    exptime = {'F184': 901.175,
           'J129': 302.275,
           'H158': 302.275,
           'K213': 901.175,
           'R062': 161.025,
           'Y106': 302.275,
           'Z087': 101.7}

    if band in exptime.keys():
        return exptime[band]
    else:
        return exptime

def _coord_transf(ra,dec):
    """
    Helper function for _sca_check and get_object_instances.get p
    Inputs must be in radians. 

    :return: Transformed x, y, z coordinates of given RA, dec. 
    :rtype: tuple  

    """
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    return x,y,z

def _distance(x0,x1,y0,y1,z0,z1):
    """Distance formula. Helper function for _sca_check and get_object_instances. 

    :param x0: _description_
    :type x0: _type_
    :param x1: _description_
    :type x1: _type_
    :param y0: _description_
    :type y0: _type_
    :param y1: _description_
    :type y1: _type_
    :param z0: _description_
    :type z0: _type_
    :param z1: _description_
    :type z1: _type_
    :return: _description_
    :rtype: _type_
    """    
    return np.sqrt((x0 - x1)**2 + (y0 - y1)**2 + (z0 - z1)**2)

def _sca_check(sca_ra, sca_dec, oid_x, oid_y, oid_z):

    """
    This is a helper function for get_object_instances and is 
    not meant to be called directly. 

    :return: Indices of which specific SCAs contain a given object.
            Returns indices in columns, where the first column is 
            SCA-1 and the second is the pointing that these SCAs
            belong to. 
    :rtype: np.array
    
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

    return idx

def _obj_in(oid,df):
    """Helper function for get_object_instances. Tells you if
        a row in a table has the specified object ID.     

    :param oid: Unique object ID. 
    :type oid: int
    :param df: _description_
    :type df: _type_
    :return: _description_
    :rtype: boolean
    """    
    if int(oid) in df['object_id'].astype(int):
        return True
    else:
        return False    

def get_object_instances(ra,dec,oid=None,bands=get_roman_bands(),
                        pointings=np.arange(0,57365,1),
                        mjd_start=-np.inf,mjd_end=np.inf):

    """
    Retrieves all images that a unique object or set of coordinates is in. There are three steps to this, because
    I think it will make the code run faster:
    1. First cut (coarse): Check object's proximity to boresight coordinates for all pointings.
    2. Second cut (fine): Of the pointings that passed the first cut, check the proximity of the center
       of each individual SCA to the object's coordinates. If an object ID is not provided, this is the 
       final step. 
    3. Of the SCAs that passed the second cut, open each truth file and check if the object ID is in it. 

    RA/dec arguments should be in degrees, as they are in the obseq file. 

    :param ra: 
    :type ra: float
    :param dec:
    :type dec: float
    :param oid: If None, return table containing SCAs that contain the provided RA, dec. If an object ID is
                provided in this field, this function checks each truth file associated with each image 
                listed in the table it assembles to ensure that the particular object is present in those 
                images. WARNING: This is slow. 
    :type oid: int or None, optional
    :param band: Filters to include in search. Default ['F184', 'H158', 'J129', 'K213', 'R062', 'Y106', 'Z087'].
    :type band: list or str, optional
    :param pointings: Limit search to particular pointings.  
    :type pointings: list or np.ndarray, optional
    :param mjd_start: Start MJD to include in search. 
    :type mjd_start: float, optional
    :param mjd_end: End MJD to include in search. 
    :type mjd_end: float, optional
    :return: Astropy table with columns filter, pointing, SCA identifying images that contain the input RA 
            and dec. If and object ID is provided in argument oid, this function checks each truth file associated 
            with each image listed in the table it assembles to ensure that the particular object is present in those 
            images. WARNING: This is slow. Note that if an object ID is not provided, the final list returned will
            contain some images where the provided RA and dec are slightly outside its bounds. This is because for
            speed, it takes the center of the SCA from the obseq files, and circumscribes a circle around the SCA for
            the search zone. 
    :rtype: astropy.table.Table
    """
    
    obseq_path = pa.join(rootdir,'Roman_TDS_obseq_11_6_23.fits')
    with fits.open(obseq_path) as osp:
        obseq_orig = Table(osp[1].data)

    obseq_radec_path = pa.join(rootdir,'Roman_TDS_obseq_11_6_23_radec.fits')
    with fits.open(obseq_radec_path) as osradecp:
        obseq_radec_orig = Table(osradecp[1].data)

    pointing_idx = np.array(pointings)
    obseq = obseq_orig[pointing_idx]
    obseq_radec = obseq_radec_orig[pointing_idx]

    band_idx = np.where(np.in1d(obseq['filter'],bands))[0]
    obseq = obseq[band_idx]
    obseq_radec = obseq_radec[band_idx]
    pointing_idx = pointing_idx[band_idx]

    mjd_idx = np.where((obseq['date'] > mjd_start) & (obseq['date'] < mjd_end))[0]
    obseq = obseq[mjd_idx]
    obseq_radec = obseq_radec[mjd_idx]
    pointing_idx = pointing_idx[mjd_idx]

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

    distance_idx = np.where(d_actual/2. <= d_req)[0]
    idx_firstcut = pointing_idx[distance_idx]

    sca_tab = obseq_radec_orig[idx_firstcut]
    sca_tab_ra = np.array(sca_tab['ra'])
    sca_tab_dec = np.array(sca_tab['dec'])
    sca_tab_filter = np.array(sca_tab['filter'])

    sca_rd = np.array([sca_tab_ra,sca_tab_dec]).reshape((2,len(sca_tab),18)).T
    sca_check = _sca_check(sca_rd[:,:,0], sca_rd[:,:,1],x_oid,y_oid,z_oid)
    pointing_idx = sca_check[:,1]
    sca_idx = sca_check[:,0]

    secondcut_tab = Table([sca_tab_filter[pointing_idx], idx_firstcut[pointing_idx], sca_idx+1], names=('filter','pointing','sca'))

    # Now, make sure the provided RA, dec are actually in these images. 
    # This is slow because it opens each file for its WCS. 
    inimg = list(map(radec_isin, [ra]*len(secondcut_tab), [dec]*len(secondcut_tab), \
                    [None]*len(secondcut_tab), secondcut_tab['filter'], secondcut_tab['pointing'], secondcut_tab['sca']))

    thirdcut_tab = secondcut_tab[inimg]

    if oid is not None:
        # This part of the code is really slow because I'm opening files. 
        # Want to parallelize in future to speed up. 
        final_idx = []
        for i, row in enumerate(thirdcut_tab):
            df = read_truth_txt(band=row['filter'],pointing=row['pointing'],sca=row['sca'])
            if _obj_in(oid, df):
                final_idx.append(i)
            del df

        tab = thirdcut_tab[final_idx]

    else:
        tab = thirdcut_tab

    return tab

def get_object_data(oid, metadata,
                    colnames=['object_id','ra','dec','mag_truth','flux_truth','flux_fit','flux_err'],
                    crossmatch_dir='/hpc/group/cosmology/lna18/roman_sim_imgs/Roman_Rubin_Sims_2024/preview/crossmatched_truth'):

    """Retrieves all information from crossmatched photmetry files about a particular object, specified with its
    unique object ID. 

    :param oid: Object ID. 
    :type oid: int
    :param metadata: Output from get_object_instance. Table with filter, pointing, and SCA numbers for this
                    function to search through. 
    :type metadata: astropy.table
    :param colnames: The columns that this function should retrieve from the crossmatched photometry files.
    :type colnames: list, optional
    :return: Astropy table with all instances of argument oid found in the crossmatched photometry files
            that are present in the metadata table. 
    :rtype: astropy.table.Table

    """
    
    object_tab = Table(names=colnames)
    for row in metadata:
        band = row['filter']
        p = row['pointing']
        sca = row['sca']

        filepath = pa.join(crossmatch_dir,f'{band}/{p}/Roman_TDS_xmatch_{band}_{p}_{sca}.txt')
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
