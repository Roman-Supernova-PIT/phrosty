import pytest # noqa: F401
import pathlib

import tox # noqa: F401
from tox.pytest import init_fixture # noqa: F401

from snpit_utils.config import Config
from snappl.image import FITSImageOnDisk
from snappl.diaobject import DiaObject
from snappl.imagecollection import ImageCollection


direc = pathlib.Path( __file__ ).parent


@pytest.fixture( scope='session' )
def dia_out_dir():
    return pathlib.Path( Config.get().value( 'photometry.phrosty.paths.dia_out_dir' ) )


@pytest.fixture( scope='session' )
def template_csv():
    return direc / "20172782_instances_templates_1.csv"


@pytest.fixture( scope='session' )
def two_science_csv():
    return direc / "20172782_instances_science_2.csv"


# @pytest.fixture( scope='session' )
# def test_dia_image():
#     # This is the first image from the csv file in 20172782_instances_science_2.csv
#     return { 'relpath': 'simple_model/Y106/35198/Roman_TDS_simple_model_Y106_35198_2.fits.gz',
#              'pointing': 35198,
#              'sca': 2,
#              'mjd': 62455.669,
#              'band': 'Y106'
#             }

# @pytest.fixture( scope='session' )
# def test_truth_file():
#     # This is the truth file corresponding to the above image in test_dia_image().
#     return pathlib.Path( Config.get().value('ou24.tds_base') ) / 'truth' /
#                          'Y106' / '35198' / 'Roman_TDS_index_Y106_35198_2.txt'


# @pytest.fixture( scope='session' )
# def test_sn():
#     # This object is on the science images in 20172782_instances_science_2.csv
#     return { 'oid': 20172782,
#              'mjd0': 62450.0,
#              'mjd1': 62881.0,
#              'ra': 7.551093401915147,
#              'dec': -44.80718106491529,
#              'zcmb': 0.3601,
#              'peakmjd': 62476.507812
#             }


@pytest.fixture( scope="session" )
def object_for_tests():
    return DiaObject.find_objects( collection='ou2024', id=20172782 )[0]


@pytest.fixture( scope="session" )
def ou2024_image_collection():
    return ImageCollection.get_collection( "ou2024" )


# These next two are session scope fixtures, so be sure not to modify
#   the things that you get from them.
@pytest.fixture
def one_science_image( scope="session" ):
    try:
        img = FITSImageOnDisk( path=('/photometry_test_data/ou2024/images/simple_model/'
                                     'Y106/35198/Roman_TDS_simple_model_Y106_35198_2.fits.gz'),
                               imhdu=1,
                               pointing=35198,
                               sca=2 ).uncompressed_version()
        yield img
    finally:
        img.path.unlink( missing_ok=True )


@pytest.fixture
def one_template_image( scope="session" ):
    try:
        img = FITSImageOnDisk( path=('/photometry_test_data/ou2024/images/simple_model/'
                                     'Y106/5934/Roman_TDS_simple_model_Y106_5934_3.fits.gz' ),
                               imhdu=1,
                               pointing=5934,
                               sca=3 ).uncompressed_version()
        yield img
    finally:
        img.path.unlink( missing_ok=True )


@pytest.fixture
def one_ou2024_template_image( ou2024_image_collection ):
    return ou2024_image_collection.get_image( pointing=5934, sca=3, band='Y106' )


@pytest.fixture
def two_ou2024_science_images( ou2024_image_collection ):
    img1 = ou2024_image_collection.get_image( pointing=35198, sca=2, band='Y106' )
    img2 = ou2024_image_collection.get_image( pointing=39790, sca=15, band='Y106' )
    return img1, img2
