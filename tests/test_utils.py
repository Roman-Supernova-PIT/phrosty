# IMPORTS Astro:
from astropy.io import fits
from astropy.table import Table

# IMPORTS Internal:
from phrosty.utils import get_mjd

def test_get_mjd():
    pointing = 0
    test_path = 'testdata/Roman_TDS_obseq_11_6_23.fits'
    with fits.open(test_path) as tp:
        tobseq = Table(tp[1].data)
    tmjd = float(tobseq['date'][int(pointing)])
    mjd = get_mjd(pointing)

    assert tmjd == mjd