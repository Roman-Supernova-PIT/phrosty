import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, vstack

def read_truth_txt(truth_path):
    truth_colnames = ['object_id', 'ra', 'dec', 'x', 'y', 'realized_flux', 'flux', 'mag', 'obj_type']
    truth_pd = pd.read_csv(truth_path, comment='#', skipinitialspace=True, sep=' ', names=truth_colnames)
    truth = Table.from_pandas(truth_pd)

    return truth

def get_object(oid,config,colnames=[['object_id','ra','dec','mag_truth','flux_truth','flux_fit','flux_err']]):
    """
    Retrieves all information for a particular object ID after you've crossmatched photometry with the corresponding
    truth file. 

    """
    subconfig = config[config['object_id'] == oid]
    print(subconfig)
    object_tab = Table(names=colnames)

    for row in subconfig:
        band = row['filter']
        p = row['pointing']
        sca = row['sca']
        filepath = f'/hpc/group/cosmology/lna18/roman_sim_imgs/Roman_Rubin_Sims_2024/preview/crossmatched_truth/{band}/{p}/Roman_TDS_xmatch_{band}_{p}_{sca}.txt'
        phot = Table.read(filepath, format='csv')
        object_row = phot[phot['object_id'] == int(oid)]
        object_row_reduced = object_row[colnames]
        object_tab = vstack([object_tab,object_row_reduced], join_type='exact')

    return object_tab

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