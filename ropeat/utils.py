import numpy as np
import pandas as pd
from numba import njit

def get_obj_type_from_ID(ID):
    if ID > 9921000000000 and ID < 10779202101973:
        return 'galaxy'
    elif ID > 30328699913 and ID < 50963307358:
        return 'star'
    elif ID > 20000001 and ID < 120026449:
        return 'transient'
    else:
        return 'unknown'

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

def make_truth_config_table(list_of_paths):
    for path in list_of_paths:
        f = open(path, 'r')
        

@njit(cache=True)
def collate_truth(truthfile,fulltab):
    """
    This should be used in a loop, like this:
    colnames = ['object_id', 'ra', 'dec', 'x', 'y', 'realized_flux', 'flux', 'mag', 'obj_type']
    fulltab = pd.DataFrame(columns=colnames).to_numpy()
    pointings = ['11776']
    for p in pointings:
        for sca in np.arange(1,19,1):
            filepath = f'/hpc/group/cosmology/RomanDESC_sims_2024/truth/H158/{p}/Roman_TDS_index_H158_{p}_{sca}.txt'
            truthfile = pd.read_csv(filepath, comment='#', skipinitialspace=True, sep=' ', names=colnames).to_numpy()
            fulltab = collate_truth(truthfile,fulltab)
    """

    fulltab = np.vstack((fulltab,truthfile))
    fulltab = fulltab[fulltab[:,0].argsort()]
    unique_ids = np.unique(fulltab[:,0])

    for uid in unique_ids:
        idx = np.where(fulltab[:,0] == uid)[0]
        subtab = fulltab[idx]

        for colidx in np.arange(3,8,1):
            vallist = subtab[:,colidx]
            fulltab[idx,colidx] = vallist

    return fulltab