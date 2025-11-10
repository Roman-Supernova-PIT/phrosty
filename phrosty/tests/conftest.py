import pytest # noqa: F401
import pathlib

import tox # noqa: F401
from tox.pytest import init_fixture # noqa: F401

from snappl.config import Config
from snappl.image import FITSImageOnDisk
from snappl.diaobject import DiaObject
from snappl.imagecollection import ImageCollection


_direc = pathlib.Path( __file__ ).parent.resolve()


@pytest.fixture( scope='session', autouse=True )
def config():
    """Set the global config for phrosty tests.

    Tests ignore SNPIT_CONFIG env var, but will always use
    phrosty/tests/phrosty_test_config.

    Output directories are made underneath test_output.  This fixture
    empties out all those output directories when the tests finish.

    """

    # Directories we'll use for test outputs
    temp_dir = _direc / 'test_output/temp'
    scratch_dir = temp_dir
    dia_out_dir = _direc / 'test_output/dia_out_dir'
    ltcv_dir = _direc / 'test_output/lc_out_dir'
    # Make sure they exist
    for path in [ temp_dir, scratch_dir, dia_out_dir, ltcv_dir ]:
        path.mkdir( exist_ok=True, parents=True )

    try:
        # We're going to mangle the config to move the output
        #  directories underneath tests so that we can clean them up
        #  when done.  Ideally, each test should clean itself up.
        # This means cheating with the config so we can modify it.  The
        #  first call to Config.get() gets a pointer to the global
        #  singleton config object (as it defaults to static=True .
        # (We don't just do this in the phrosty_test_config.yaml file
        #   because some examples use that and need access to the
        #   directories that the tests don't destroy.)
        cfg = Config.get( _direc / 'phrosty_test_config.yaml', setdefault=True )
        # Now we're going to poke inside the Config object so we can
        #   modify this global singleton.  We're not supposed to do
        #   that.  If this was Java, it totally wouldn't let us.  But
        #   here we go.
        cfg._static = False
        # Edit the config.
        cfg.set_value( 'photometry.snappl.temp_dir', str(temp_dir) )
        cfg.set_value( 'photometry.phrosty.paths.scratch_dir', str(temp_dir) )
        cfg.set_value( 'photometry.phrosty.paths.temp_dir', str(temp_dir) )
        cfg.set_value( 'photometry.phrosty.paths.dia_out_dir', str(dia_out_dir) )
        cfg.set_value( 'photometry.phrosty.paths.ltcv_dir', str(ltcv_dir) )
        # Reset the config to static
        cfg._static = True

        yield cfg

    finally:
        # Clean up the output directories.  (Don't delete the
        #   directories themselves, but do delete all files and
        #   subdirectories.)
        # TODO: we could put a check here that they're actually empty,
        #   and yell if they're not.  That means some test didn't clean
        #   up after itself.
        def nukedir( path ):
            for f in path.iterdir():
                if f.is_dir():
                    nukedir( f )
                    f.rmdir()
                else:
                    f.unlink( missing_ok=True )

        nukedir( temp_dir )
        nukedir( scratch_dir )
        nukedir( dia_out_dir )
        nukedir( ltcv_dir )


@pytest.fixture( scope='session' )
def template_csv():
    return _direc / "20172782_instances_templates_1.csv"


@pytest.fixture( scope='session' )
def two_science_csv():
    return _direc / "20172782_instances_science_2.csv"


@pytest.fixture( scope="session" )
def dia_out_dir( config ):
    return pathlib.Path( config.value( 'photometry.phrosty.paths.dia_out_dir' ) )

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
    return DiaObject.find_objects( collection='ou2024', name=20172782 )[0]


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
                               imagehdu=1,
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
