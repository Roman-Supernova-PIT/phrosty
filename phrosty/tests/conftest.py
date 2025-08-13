import pytest # noqa: F401
import pathlib

import tox # noqa: F401
from tox.pytest import init_fixture # noqa: F401

from snpit_utils.config import Config

import phrosty.imagesubtraction

direc = pathlib.Path( __file__ ).parent


@pytest.fixture( scope='session' )
def dia_out_dir():
    return pathlib.Path( Config.get().value( 'photometry.phrosty.paths.dia_out_dir' ) )


@pytest.fixture( scope='session' )
def sims_dir():
    return pathlib.Path( Config.get().value( 'ou24.tds_base' ) ) / 'images'


# @pytest.fixture( scope='session' )
# def sn_info_dir():
#     return pathlib.Path( os.getenv( "SN_INFO_DIR", direc / "sn_info_dir" ) )


@pytest.fixture( scope='session' )
def template_csv():
    return direc / "20172782_instances_templates_1.csv"


@pytest.fixture( scope='session' )
def two_science_csv():
    return direc / "20172782_instances_science_2.csv"


@pytest.fixture( scope='session' )
def test_dia_image():
    # This is the first image from the csv file in 20172782_instances_science_2.csv
    return { 'relpath': 'simple_model/Y106/35198/Roman_TDS_simple_model_Y106_35198_2.fits.gz',
             'pointing': 35198,
             'sca': 2,
             'mjd': 62455.669,
             'band': 'Y106'
            }

# @pytest.fixture( scope='session' )
# def test_truth_file():
#     # This is the truth file corresponding to the above image in test_dia_image().
#     return pathlib.Path( Config.get().value('ou24.tds_base') ) / 'truth' /
#                          'Y106' / '35198' / 'Roman_TDS_index_Y106_35198_2.txt'


@pytest.fixture( scope='session' )
def test_sn():
    # This object is on the science images in 20172782_instances_science_2.csv
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
        assert line.split() == [ 'path', 'pointing', 'sca', 'mjd', 'band' ]
        line = ifp.readline()
        img, _point, _sca, _mjd, _band = line.split()

    return sims_dir / img


@pytest.fixture( scope='session' )
def one_compressed_science_image_path( sims_dir, two_science_csv ):
    with open( two_science_csv ) as ifp:
        line = ifp.readline()
        assert line.split() == [ 'path', 'pointing', 'sca', 'mjd', 'band' ]
        line = ifp.readline()
        img, _point, _sca, _mjd, _band = line.split()

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
