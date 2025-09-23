import os
import numpy as np
import pytest
from phrosty.pipeline import Pipeline

# TODO : separate tests for PipelineImage, for all the functions in#
#   PipelineImage and Pipeline.  Right now we just have this regression
#   test.


@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0 ), reason="SKIP_GPU_TESTS is set" )
def test_pipeline_run( object_for_tests, ou2024_image_collection,
                       one_ou2024_template_image, two_ou2024_science_images ):
    pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                    science_images=two_ou2024_science_images,
                    template_images=[one_ou2024_template_image],
                    nprocs=2, nwrite=3 )
    ltcv = pip()

    with open( ltcv ) as ifp:
        hdrline = ifp.readline().strip()
        assert hdrline == ( 'ra,dec,mjd,filter,pointing,sca,template_pointing,template_sca,zpt,'
                            'aperture_sum,flux_fit,flux_fit_err,mag_fit,mag_fit_err,success' )
        kws = hdrline.split( "," )
        pairs = []
        for line in ifp:
            data = line.strip().split( "," )
            pairs.append( { kw: val for kw, val in zip( kws, data ) } )

    assert len(pairs) == 2
    for pair, img in zip( pairs, two_ou2024_science_images ):
        assert float(pair['ra']) == pytest.approx( object_for_tests.ra, abs=0.1/3600. )
        assert float(pair['dec']) == pytest.approx( object_for_tests.dec, abs=0.1/3600. )
        assert float(pair['mjd']) == pytest.approx( img.mjd, abs=1e-3 )
        assert pair['filter'] == 'Y106'
        assert int(pair['pointing']) == int(img.pointing)
        assert int(pair['sca']) == int(img.sca)
        assert int(pair['template_pointing']) == int(one_ou2024_template_image.pointing)
        assert int(pair['template_sca']) == int(one_ou2024_template_image.sca)
        # TODO : fix the zeropoint! This tolerance is huge...
        zpt = float( pair['zpt'] )
        assert zpt == pytest.approx( 32, abs=1 )
        assert pair['success'] == 'True'

    # Tests aren't exactly reproducible from one run to the next,
    #   because some classes (including the galsim PSF that we use right
    #   now) have random numbers in them, and at the moment we aren't
    #   controlling the seed.  So, we can only test for approximately
    #   consistent results.  Going to do 0.3 times the uncertainty,
    #   because a difference by that much is not all that meaningful
    #   change, and empirically they vary by that much.  (Which is
    #   alarming, but what can you do.)

    dflux = float( pairs[0]['flux_fit_err'] )
    assert dflux == pytest.approx( 25., rel=0.3 )
    dmag = float( pairs[0]['mag_fit_err'] )
    assert dmag == pytest.approx( 0.15, abs=0.1 )
    assert float( pairs[0]['aperture_sum'] ) == pytest.approx( 1006.8969554640036, abs=0.3*dflux )
    assert float( pairs[0]['flux_fit'] ) == pytest.approx( 181.9182196835094, abs=0.3*dflux )
    assert float( pairs[0]['mag_fit'] ) == pytest.approx( -5.64, abs=max( 0.3*dmag, 0.01 ) )

    dflux = float( pairs[1]['flux_fit_err'] )
    assert dflux == pytest.approx( 27., rel=0.3 )
    dmag = float( pairs[1]['mag_fit_err'] )
    assert dmag == pytest.approx( 0.05, abs=0.1 )
    assert float( pairs[1]['aperture_sum'] ) == pytest.approx( 4112, abs=0.3*dflux )
    assert float( pairs[1]['flux_fit'] ) == pytest.approx( 721, abs=0.3*dflux )
    assert float( pairs[1]['mag_fit'] ) == pytest.approx( -7.14, abs=max( 0.3*dmag, 0.01 ) )

    # TODO : cleanup output directories!  This is scary if you're using the same
    #   directories for tests and for running... so don't do that... but the
    #   way we're set up right now, you probably are.


@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0 ), reason="SKIP_GPU_TESTS is set" )
def test_pipeline_failures( object_for_tests, ou2024_image_collection,
                            one_ou2024_template_image, two_ou2024_science_images,
                            nan_image ):
    pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                    science_images=two_ou2024_science_images,
                    template_images=[one_ou2024_template_image],
                    nprocs=2, nwrite=3 )

    # First, check the images as-is. Make sure there are no failures.
    for key in pip.failures:
        assert len(pip.failures[key]) == 0
    
    new_test_imgs = [nan_image, two_ou2024_science_images[1]]

    pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                science_images=new_test_imgs,
                template_images=[one_ou2024_template_image],
                nprocs=2, nwrite=3 )

    for key in pip.failures:
        print(key)
        print(len(pip.failures[key]))