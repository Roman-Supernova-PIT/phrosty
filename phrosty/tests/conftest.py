import pytest # noqa: F401
import os
import pathlib

import tox # noqa: F401
from tox.pytest import init_fixture # noqa: F401

import phrosty.imagesubtraction

direc = pathlib.Path( __file__ ).parent

@pytest.fixture( scope='session' )
def dia_out_dir():
    return pathlib.Path( os.getenv( "DIA_OUT_DIR", "/dia_out_dir" ) )

@pytest.fixture( scope='session' )
def sims_dir():
    return pathlib.Path( os.getenv( "SIMS_DIR", "/sims_dir" ) )

@pytest.fixture( scope='session' )
def sn_info_dir():
    return pathlib.Path( os.getenv( "SN_INFO_DIR", direc / "sn_info_dir" ) )

@pytest.fixture( scope='session' )
def template_csv():
    return direc / "20172782_instances_templates_1.csv"

@pytest.fixture( scope='session' )
def two_science_csv():
    return direc / "20172782_instances_science_two_images.csv"

@pytest.fixture( scope='session' )
def test_dia_image():
    # This is the first image from the csv file in two_science_csv
    return { 'relpath': 'RomanTDS/images/simple_model/R062/35083/Roman_TDS_simple_model_R062_35083_8.fits.gz',
             'pointing': 35083,
             'sca': 8,
             'mjd': 62455.174,
             'band': 'R062'
            }

@pytest.fixture( scope='session' )
def test_sn():
    # This object is on the science images in two_science_csv
    return { 'oid': 20172782,
             'mjd0': 62450.0,
             'mjd1': 62881.0,
             'ra': 7.551093401915147,
             'dec': -44.80718106491529,
             'zcmb': 0.3601,
             'peakmjd': 62476.507812
            }
             

@pytest.fixture( scope='session' )
def compressed_template_image_path( sims_dir, template_csv ):
    with open( template_csv ) as ifp:
        line = ifp.readline()
        img, point, sca, mjd = line.split()

    return sims_dir / img

@pytest.fixture( scope='session' )
def one_compressed_science_image_path( sims_dir, two_science_csv ):
    with open( two_science_csv ) as ifp:
        line = ifp.readline()
        img, point, sca, mjd = line.split()

    return sims_dir / img


@pytest.fixture( scope='session' )
def template_image_path( compressed_template_image_path, dia_out_dir ):
    out_path = dia_out_dir / "templ.fits"
    phrosty.imagesubtraction.gz_and_ext( compressed_template_image_path, out_path )

    try:
        yield out_path

    finally:
        out_path.unlink( missing_ok=True )
    
@pytest.fixture( scope='session' )
def one_science_image_path( one_compressed_science_image_path, dia_out_dir ):
    out_path = dia_out_dir / "sci.fits"
    phrosty.imagesubtraction.gz_and_ext( one_compressed_science_image_path, out_path )

    try:
        yield out_path

    finally:
        out_path.unlink( missing_ok=True )
