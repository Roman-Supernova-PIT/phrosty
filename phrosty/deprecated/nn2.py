# IMPORTS Standard:
import os
import numpy as np

# IMPORTS Astro:
from astropy.table import Table, join

r"""
See https://iopscience.iop.org/article/10.1086/491583.
Good matrix help PDF: https://www.stat.cmu.edu/~cshalizi/mreg/15/lectures/13/lecture-13.pdf

n := number of data points in an LC.
N := n(n-1)/2

v := Free parameter. Final LC. Length n.
A := vector of differences between LC data points, V_{i} - V_{j} of length N.
X := (N,n) matrix containing 0, 1, and -1 that turns V_{i} and V_{j} "on" and "off".

V = xA

Want to minimize:

chi^{2} = (A - Xv)^{T} C^{-1} (A - Xv)

Then,
d\chi^{2} = (C + C^{-1})(A - Xv)
...
v = (X^{T} C^{-1} X)^{-1} X^{T} C^{-1} A

1. Turn photometry csv files into differences.
2. Use step 1 to construct A.

NOTE: This method requires that all input LCs have the same length.

Thanks to Bastien Carreres for help with the matrix operations.

"""


def parse_csvs(filepaths):
    """Read in csv files. Merge tables.

    'filepaths' is a list of paths.
    Collects all epochs of observations.

    """

    overlap_cols = ['filter', 'pointing', 'sca', 'mjd', 'ra_init', 'dec_init', 'x_init', 'y_init']
    unique_cols = ['ra_fit', 'dec_fit', 'x_fit', 'y_fit', 'mag', 'magerr', 'flux', 'fluxerr', 'zpt']
    combined_table = Table(names=overlap_cols)

    for f in filepaths:
        basename = os.path.basename(f)
        basename_split = basename.split('_')
        band, pointing, sca = basename_split[1],  basename_split[2], basename_split[3]
        colname_extension = f'{band}_{pointing}_{sca}'
        lc_table = Table.read(f, format='csv')

        for col in unique_cols:
            f.rename_column(col, f'{col}_{colname_extension}')

        # Only combine on shared fixed information.
        combined_table = join(combined_table, lc_table, keys=overlap_cols, join_type='inner')

    return combined_table


def construct_A_and_C(combined_table):
    """Input table from parse_csvs.

    Constructs vector A and diagonal elements of C for all epochs of
    observation.

    """

    colnames = combined_table.colnames
    flux_table =  combined_table[[col for col in colnames if 'flux_' in col]]
    fluxerr_table = combined_table[[col for col in colnames if 'fluxerr_' in col]]

    flux_array = np.array([flux_table[col].data for col in flux_table.colnames])
    fluxerr_array = np.array([fluxerr_table[col].data for col in fluxerr_table.colnames])

    i, j = np.triu_indices(len(flux_table), k=1)

    A_all = []
    C_diag_all = []

    for epoch_flux, epoch_fluxerr in zip(flux_array, fluxerr_array):
        A_all.append(epoch_flux[i]-epoch_flux[j])
        C_diag_all.append(epoch_fluxerr[i], epoch_fluxerr[j])

    # for idx in i:
    #     for jdx in j:
    #         dv = flux_array[idx] - flux_array[jdx]
    #         sv = np.sqrt(fluxerr_array[idx]**2 + fluxerr_array[jdx]**2)

    #         A_all.append(dv)
    #         C_diag_all.append(sv)

    return A_all, C_diag_all


def construct_X(n, N):
    """Construct matrix X. Applies to one epoch at a time."""
    i, j = np.triu_indices(n, k=1)
    X = np.zeros((N, n))
    np.put_along_axis(X, i[:, None], 1, axis=1)
    np.put_along_axis(X, j[:, None], -1, axis=1)

    return X


def solve(A, X, C):
    """Solve one epoch.

    v = (X^{T} C^{-1} X)^{-1} X^{T} C^{-1} A

    Final solution will be (original vector of flux differences - v).

    """

    XTCXXT = np.linalg.lstsq(X.T @ np.linalg.inv(C) @ X, X.T, rcond=None)
    v = XTCXXT @ np.linalg.inv(C) @ A

    return v
