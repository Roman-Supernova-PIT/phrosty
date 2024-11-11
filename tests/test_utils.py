# IMPORTS Standard:
import os
import os.path
import re
import pathlib
import pytest

import numpy as np

# IMPORTS Astro:
from astropy.io import fits
from astropy.table import Table

# IMPORTS Internal:
import phrosty.utils
from phrosty.utils import get_mjd, get_transient_zcmb, get_transient_peakmjd, _build_filepath

def test_utils_globals():
    assert pathlib.Path( phrosty.utils.rootdir ).is_dir()
    assert pathlib.Path( phrosty.utils.snana_pq_dir ).is_dir()

    tab = Table.read( phrosty.utils.obseq_path )
    assert len(tab) == 57365
    assert set( tab.columns ) == {'ra', 'pa', 'date', 'filter', 'exptime', 'dec'}

    tab = Table.read( phrosty.utils.obseq_radec_path )
    assert len(tab) == 57365
    assert set( tab.columns ) == {'ra', 'dec', 'filter'}


def test_build_filepath():
    band = 'R062'
    pointing = 51235
    sca = 11

    with pytest.raises( ValueError, match=r"filetype must be in \['image', 'truth', 'truthtxt'\]" ):
        _build_filepath( None, band, pointing, sca, 'foo' )

    for args in [ [ None, None, pointing, sca, 'image' ],
                  [ None, band, None, sca, 'image' ],
                  [ None, band, pointing, None, 'image' ] ]:
        with pytest.raises( ValueError, match="You need to specify band" ):
            _build_filepath( *args )

    types = [ 'image', 'truth', 'truthtxt' ]
    for typ in types:
        fp = _build_filepath( None, band, pointing, sca, typ )
        assert re.search( f"^Roman_TDS_.*_{band}_{pointing}_{sca}", os.path.basename(fp) )

    fp = _build_filepath( "/foo/bar", None, None, None, 'image' )
    assert fp == '/foo/bar'


def test_get_roman_bands():
    assert phrosty.utils.get_roman_bands() == [ 'R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'K213' ]


def test_read_truth_txt():
    band = 'R062'
    pointing = 51235
    sca = 11
    truth = phrosty.utils.read_truth_txt( None, band, pointing, sca )
    assert len(truth) == 24628
    assert set(truth.columns) == {'dec', 'ra', 'obj_type', 'flux', 'x', 'realized_flux', 'object_id', 'mag', 'y'}


def test_radec_isin():
    band = 'R062'
    pointing = 35083
    sca = 8
    rain = 7.551093401915147
    decin = -44.80718106491529
    raout = 7.4106
    decout = -44.7131
    rawayout = 8
    decwayout = 32

    assert phrosty.utils.radec_isin( rain, decin, path=None, band=band, pointing=pointing, sca=sca )
    assert not phrosty.utils.radec_isin( raout, decout, path=None, band=band, pointing=pointing, sca=sca )
    assert not phrosty.utils.radec_isin( rawayout, decwayout, path=None, band=band, pointing=pointing, sca=sca )


def test_get_corners():
    band = 'R062'
    pointing = 35083
    sca = 8

    corners = phrosty.utils.get_corners( path=None, band=band, pointing=pointing, sca=sca )
    expected = ( np.array([  7.58004938, -44.71290687]),
                 np.array([  7.5896352 , -44.83369681]),
                 np.array([  7.4075831 , -44.71566849]),
                 np.array([  7.4168088 , -44.83647236]) )
    assert( np.all( np.abs(c-e) < 0.00028 ) for c, e in zip( corners, expected ) )


# TEST _read_parquet?

# This one also implicitly tests make_object_table
def test_get_transient_radec():
    trns = 20172782
    coords = (7.551093401915147, -44.80718106491529)
    assert phrosty.utils.get_transient_radec(trns) == pytest.approx( coords, abs=0.000028 )


def test_get_transient_mjd():
    trns = 20172782
    endpoints = (62450.0, 62881.0)
    assert phrosty.utils.get_transient_mjd(trns) == pytest.approx( endpoints, abs=0.1 )

def test_get_transient_zcmb():
    oid = 20172782
    zcmb = np.float32(0.3601)

    tzcmb = get_transient_zcmb(oid)

    assert zcmb == pytest.approx( tzcmb, rel=1e-5 )

def test_get_transient_peakmjd():
    oid = 20172782
    mjd = np.float32(62476.507812)

    tmjd = get_transient_peakmjd(oid)

    assert mjd == pytest.approx( tmjd, abs=0.001 )


def test_get_transient_info():
    oid = 20172782
    coords = (7.551093401915147, -44.80718106491529, 62450.0, 62881.0)
    assert phrosty.utils.get_transient_info( oid ) == pytest.approx( coords, rel=1e-6 )


def test_transient_in_or_out():
    oid = 20172782
    mjd0 = 62450.0
    mjd1 = 62881.0
    band = 'R062'

    in_tab, out_tab = phrosty.utils.transient_in_or_out( oid, mjd0, mjd1, band )
    assert len(in_tab) == 53
    assert len(out_tab) == 82
    assert set(in_tab.columns) == { 'filter', 'pointing', 'sca' }
    assert set(out_tab.columns) == { 'filter', 'pointing', 'sca' }

@pytest.mark.skip( reason="Not used in pipeline.py, test is a TODO" )
def test_get_templates():
    assert False

@pytest.mark.skip( reason="Not used in pipeline.py, test is a TODO" )
def test_get_science():
    assert False

def test_get_mjd_limits():
    start, end = phrosty.utils.get_mjd_limits()
    assert ( start, end ) == pytest.approx( (62000.02139, 63563.0579), abs=0.1 )

def test_get_radec_limits():
    val = phrosty.utils.get_radec_limits()
    assert val['ra'] == pytest.approx( [6.97879, 12.0204], abs=0.001 )
    assert val['dec'] == pytest.approx( [-46.5199, -41.4786], abs=0.001 )


def test_get_mjd():
    pointing = 0
    test_path = os.path.join(os.path.dirname(__file__), 'testdata', 'Roman_TDS_obseq_11_6_23.fits')
    with fits.open(test_path) as tp:
        tobseq = Table(tp[1].data)
    tmjd = float(tobseq['date'][int(pointing)])
    mjd = get_mjd(pointing)

    assert tmjd == mjd


def test_pointings_near_mjd():
    obseq = Table.read( phrosty.utils.obseq_path )

    ptgs = phrosty.utils.pointings_near_mjd( 62000.04011 )
    assert len(ptgs) == 383
    assert np.all( obseq[ptgs]['date'] > 62000.04011 - 3 )
    assert np.all( obseq[ptgs]['date'] < 62000.04011 + 3 )
    ptgs = phrosty.utils.pointings_near_mjd( 62000.04011, window=1 )
    assert len(ptgs) == 220
    assert np.all( obseq[ptgs]['date'] > 62000.04011 - 1 )
    assert np.all( obseq[ptgs]['date'] < 62000.04011 + 1 )

    ptgs = phrosty.utils.pointings_near_mjd( 10000 )
    assert len(ptgs) == 0

def test_get_mjd_info():
    mjd0 = 62000.04011 - 3
    mjd1 = 62000.04011 + 3

    info = phrosty.utils.get_mjd_info( mjd0, mjd1 )
    assert len(info) == 383
    assert set( info['filter'] ) == { 'H158', 'F184', 'R062', 'K213', 'Y106', 'Z087', 'J129' }


def test_get_exptime():
    expected = {'F184': 901.175,
                'J129': 302.275,
                'H158': 302.275,
                'K213': 901.175,
                'R062': 161.025,
                'Y106': 302.275,
                'Z087': 101.7}

    assert phrosty.utils.get_exptime() == expected

    for band, t in expected.items():
        assert phrosty.utils.get_exptime( band ) == pytest.approx( t, abs=0.001 )
