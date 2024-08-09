# IMPORTS Standard:
import os
import numpy as np

# IMPORTS Astro:
from astropy.io import fits
from astropy.table import Table

# IMPORTS Internal:
from phrosty.utils import get_mjd, get_transient_zcmb, get_transient_peakmjd

def test_get_mjd():
    pointing = 0
    test_path = os.path.join(os.path.dirname(__file__), 'testdata', 'Roman_TDS_obseq_11_6_23.fits')
    with fits.open(test_path) as tp:
        tobseq = Table(tp[1].data)
    tmjd = float(tobseq['date'][int(pointing)])
    mjd = get_mjd(pointing)

    assert tmjd == mjd

def test_get_transient_zcmb():
    oid = 20172782
    zcmb = np.float32(0.3601)

    tzcmb = get_transient_zcmb(oid)

    assert zcmb == tzcmb

def test_get_transient_peakmjd():
    oid = 20172782
    mjd = np.float32(62476.507812)

    tmjd = get_transient_peakmjd(oid)

    assert mjd == tmjd