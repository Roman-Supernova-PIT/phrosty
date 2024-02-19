import numpy as np
import pandas as pd
from numba import vectorize, float64
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, vstack
from astropy import units as u

def read_truth_txt(truth_path):
    truth_colnames = ['object_id', 'ra', 'dec', 'x', 'y', 'realized_flux', 'flux', 'mag', 'obj_type']
    truth_pd = pd.read_csv(truth_path, comment='#', skipinitialspace=True, sep=' ', names=truth_colnames)
    truth = Table.from_pandas(truth_pd)

    return truth

def get_corners(band,pointing,sca):
    imgpath = f'/cwork/mat90/RomanDESC_sims_2024/RomanTDS/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz'
    hdu = fits.open(imgpath)
    wcs = WCS(hdu[1].header)
    corners = [[0,0],[0,4088],[4088,0],[4088,4088]]
    wcs_corners = wcs.pixel_to_world_values(corners)

    return wcs_corners

def get_mjd(pointing):
    obseq_path = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/Roman_TDS_obseq_11_6_23.fits'
    obseq_hdu = fits.open(obseq_path)
    obseq = Table(obseq_hdu[1].data)
    mjd = obseq['date'][pointing]

    return mjd

def get_object_instances(oid,ra,dec,colnames=[['object_id','ra','dec','mag_truth','flux_truth','flux_fit','flux_err']]):

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

    ra_oid = (ra * u.deg).to(u.rad)
    dec_oid = (dec * u.deg).to(u.rad)

    x_oid = np.cos(dec_oid)*np.cos(ra_oid)
    y_oid = np.cos(dec_oid)*np.sin(ra_oid)
    z_oid = np.sin(dec_oid)

    ra_obseq = (obseq['ra'] * u.deg).to(u.rad)
    dec_obseq = (obseq['dec'] * u.deg).to(u.rad)

    x = np.cos(np.array(dec_obseq))*np.cos(np.array(ra_obseq))
    y = np.cos(np.array(dec_obseq))*np.sin(np.array(ra_obseq))
    z = np.sin(np.array(dec_obseq))

    d_actual = np.sqrt((x - x_oid)**2 + (y - y_oid)**2 + (z - z_oid)**2)
    d_req = np.sin(0.009/2.) # Why? Because in roman_imsim.telescope.near_pointing(), 
                             # this is the required distance. See where it says
                             # self.sbore2 = np.sin(max_rad_from_boresight/2.), and
                             # earlier, defines max_rad_from_boresight = 0.009 by default.
                             # This part of my code is also basically a copy of near_pointing. 
    
    idx_firstcut = np.where(d_actual/2. <= d_req)[0]

    sca_coords = obseq_radec[idx_firstcut]


    return idx_firstcut 

def sca_check(sca_ra, sca_dec, oid_ra, oid_dec):
    # Each SCA is 0.11 arcsec/px and 4088 x 4088 px. 
    # So, at most, the object will be the distance between the
    # center coordinate and the corner of the SCA away from the 
    # center. 

    sca_ra_as = (sca_ra * u.deg).to(u.arcsec)
    sca_dec_as = (sca_dec * u.deg).to(u.arcsec)

    sca_halfsize = 0.11 * 4088 / 2

    d_req = (np.sqrt(2)*sca_halfsize * u.arcsec).to(u.rad)



# def get_object(oid,config,colnames=[['object_id','ra','dec','mag_truth','flux_truth','flux_fit','flux_err']]):
#     """
#     Retrieves all information for a particular object ID after you've crossmatched photometry with the corresponding
#     truth file. 

#     """
#     subconfig = config[config['object_id'] == oid]
#     print(subconfig)
#     object_tab = Table(names=colnames)

    # for row in subconfig:
    #     band = row['filter']
    #     p = row['pointing']
    #     sca = row['sca']
    #     filepath = f'/hpc/group/cosmology/lna18/roman_sim_imgs/Roman_Rubin_Sims_2024/preview/crossmatched_truth/{band}/{p}/Roman_TDS_xmatch_{band}_{p}_{sca}.txt'
    #     phot = Table.read(filepath, format='csv')
    #     object_row = phot[phot['object_id'] == int(oid)]
    #     object_row_reduced = object_row[colnames]
    #     object_tab = vstack([object_tab,object_row_reduced], join_type='exact')

    # return object_tab


# def get_obj_type_from_ID(ID):
# THIS IS NO LONGER TRUE. CHANGE ID BOUNDS.
#     if ID > 9921000000000 and ID < 10779202101973:
#         return 'galaxy'
#     elif ID > 30328699913 and ID < 50963307358:
#         return 'star'
#     elif ID > 20000001 and ID < 120026449:
#         return 'transient'
#     else:
#         return 'unknown'

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

# def make_truth_config_table(list_of_paths,n_jobs=20,verbose=False):
#     colnames = ['object_id', 'ra', 'dec', 'x', 'y', 'realized_flux', 'flux', 'mag', 'obj_type']
#     config_dict = {'object_id': [],
#                    'band':      [],
#                    'pointing':  [],
#                    'sca':       []}
#     if verbose:
#         print('We need to get through', len(list_of_paths), 'images.')
#         counter = 1
#     for path in list_of_paths:
#         band = path.split('_')[-3]
#         pointing = path.split('_')[-2]
#         sca = path.split('_')[-1].split('.')[0]
#         if verbose:
#             print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
#             print(f'This is image {counter}/{len(list_of_paths)}.')
#             print(band, pointing, sca)
#             f = open(path,'r')
#             lines = f.readlines()
#             for l in lines[1:]:
#                 oid = l.split()[0]
#                 config_dict['object_id'].append(oid)
#                 config_dict['band'].append(band)
#                 config_dict['pointing'].append(pointing)
#                 config_dict['sca'].append(sca)

#             if verbose:
#                 counter += 1
#                 print(f'There are {len(config_dict["object_id"])} rows in the table.')
#     if verbose:
#         print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
#         print('Done pulling OIDs from files.')
#         print('Now, add the config values.')

#     # Update this to do something that makes sense about the object IDs and the number of jobs. 
#     oid_config = pd.DataFrame(config_dict).sort_values('object_id',ignore_index=True)
#     bins = sorted(list(np.linspace(10**6,10**8,10)) + list(np.linspace(10**8,10**12,20)) + list(np.linspace(10**12,10**13,30)))
#     config_array = np.digitize(oid_config['object_id'].astype(int),bins)
#     oid_config['config'] = config_array

#     if verbose:
#         print('Added the config column!')
    
#     return oid_config