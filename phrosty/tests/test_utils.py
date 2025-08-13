# IMPORTS Standard:
import numpy as np
import os
import os.path
import re
import pytest

# IMPORTS Astro:
from astropy.io import fits
from astropy.table import Table

# IMPORTS Internal:
import phrosty.utils
from snpit_utils.config import Config


def test_build_filepath( test_dia_image ):
    band = test_dia_image[ 'band' ]
    pointing = test_dia_image[ 'pointing' ]
    sca = test_dia_image[ 'sca' ]

    with pytest.raises( ValueError, match=r"filetype must be in \['image', 'truth', 'truthtxt'\]" ):
        phrosty.utils._build_filepath( None, band, pointing, sca, 'foo' )

    for args in [ [ None, None, pointing, sca, 'image' ],
                  [ None, band, None, sca, 'image' ],
                  [ None, band, pointing, None, 'image' ] ]:
        with pytest.raises( ValueError, match="You need to specify band" ):
            phrosty.utils._build_filepath( *args )

    types = [ 'image', 'truth', 'truthtxt' ]
    for typ in types:
        fp = phrosty.utils._build_filepath( None, band, pointing, sca, typ )
        assert re.search( f"^Roman_TDS_.*_{band}_{pointing}_{sca}", os.path.basename(fp) )

    fp = phrosty.utils._build_filepath( "/foo/bar", None, None, None, 'image' )
    assert fp == '/foo/bar'


def test_ou2024_obseq_path():
    test_path_none = os.path.join( Config.get().value('ou24.tds_base'), 'Roman_TDS_obseq_11_6_23.fits' )
    test_path_kwarg = '/foo/bar'

    path_none = phrosty.utils.ou2024_obseq_path()
    path_kwarg = phrosty.utils.ou2024_obseq_path( test_path_kwarg )

    assert path_none == test_path_none
    assert path_kwarg == test_path_kwarg


def test_get_roman_bands():
    assert phrosty.utils.get_roman_bands() == [ 'R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'K213' ]


def test_read_truth_txt( test_dia_image ):
    band = test_dia_image[ 'band' ]
    pointing = test_dia_image[ 'pointing' ]
    sca = test_dia_image[ 'sca' ]
    truth = phrosty.utils.read_truth_txt( None, band, pointing, sca )
    assert len(truth) == 25622
    assert set(truth.columns) == {'dec', 'ra', 'obj_type', 'flux', 'x', 'realized_flux', 'object_id', 'mag', 'y'}


def test_radec_isin( test_sn, test_dia_image ):
    band = test_dia_image[ 'band' ]
    pointing = test_dia_image[ 'pointing' ]
    sca = test_dia_image[ 'sca' ]
    rain = test_sn[ 'ra' ]
    decin = test_sn[ 'dec' ]
    # These next few are based on knowledge of what test_dia_image is
    raout = 7.4106
    decout = -44.7131
    rawayout = 8
    decwayout = 32

    assert phrosty.utils.radec_isin( rain, decin, path=None, band=band, pointing=pointing, sca=sca )
    assert not phrosty.utils.radec_isin( raout, decout, path=None, band=band, pointing=pointing, sca=sca )
    assert not phrosty.utils.radec_isin( rawayout, decwayout, path=None, band=band, pointing=pointing, sca=sca )


def test_get_corners( test_dia_image ):
    band = test_dia_image[ 'band' ]
    pointing = test_dia_image[ 'pointing' ]
    sca = test_dia_image[ 'sca' ]

    corners = phrosty.utils.get_corners( path=None, band=band, pointing=pointing, sca=sca )
    expected = ( np.array([  7.58004938, -44.71290687]),
                 np.array([  7.5896352 , -44.83369681]),
                 np.array([  7.4075831 , -44.71566849]),
                 np.array([  7.4168088 , -44.83647236]) )
    assert np.all( np.abs(c-e) < 0.00028 for c, e in zip( corners, expected ) )


# TEST _read_parquet?


# This one also implicitly tests make_object_table
def test_get_transient_radec( test_sn ):
    trns = test_sn[ 'oid' ]
    coords = ( test_sn['ra'], test_sn['dec'] )
    assert phrosty.utils.get_transient_radec(trns) == pytest.approx( coords, abs=0.000028 )


def test_get_transient_mjd( test_sn ):
    trns = test_sn[ 'oid' ]
    endpoints = (test_sn['mjd0'], test_sn['mjd1'])
    assert phrosty.utils.get_transient_mjd(trns) == pytest.approx( endpoints, abs=0.1 )


def test_get_transient_zcmb( test_sn ):
    oid = test_sn[ 'oid' ]
    zcmb = np.float32( test_sn['zcmb'] )

    tzcmb = phrosty.utils.get_transient_zcmb(oid)

    assert zcmb == pytest.approx( tzcmb, rel=1e-5 )


def test_get_transient_peakmjd( test_sn ):
    oid = test_sn[ 'oid' ]
    mjd = np.float32( test_sn[ 'peakmjd' ] )

    tmjd = phrosty.utils.get_transient_peakmjd(oid)

    assert mjd == pytest.approx( tmjd, abs=0.001 )


def test_get_transient_info( test_sn ):
    oid = test_sn[ 'oid' ]
    coords = (test_sn['ra'], test_sn['dec'], test_sn['mjd0'], test_sn['mjd1'])
    assert phrosty.utils.get_transient_info( oid ) == pytest.approx( coords, rel=1e-6 )


def test_transient_in_or_out( test_sn, test_dia_image ):
    oid = test_sn[ 'oid' ]
    mjd0 = test_sn[ 'mjd0' ]
    mjd1 = test_sn[ 'mjd1' ]
    band = test_dia_image[ 'band' ]

    in_tab, out_tab = phrosty.utils.transient_in_or_out( oid, mjd0, mjd1, band )
    # These numbers are based on the choice of object in test_sn
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
    mjd = phrosty.utils.get_mjd(pointing)

    assert tmjd == mjd


def test_pointings_near_mjd():
    obseq = Table.read( phrosty.utils.ou2024_obseq_path() )

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
    # These are the exposure times that were used in the OpenUniverse sims
    # (See Troxel et al. 2025 <put in the ref here when it exists>)
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
