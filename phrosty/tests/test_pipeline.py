import os
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
                            'aperture_sum,flux_fit,flux_fit_err,mag_fit,mag_fit_err' )
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
        # TODO : fix the zeropoint!
        assert pair['zpt'] == ''

    # NOTE : tests fail right now.  Results don't seem to be stable
    #   from one run to the next.  Investigate.

    assert float( pairs[0]['aperture_sum'] ) == pytest.approx( 1006.8969554640036, rel=1e-5 )
    assert float( pairs[0]['flux_fit'] ) == pytest.approx( 180.9182196835094, rel=1e-5 )
    assert float( pairs[0]['flux_fit_err'] ) == pytest.approx( 25.35943500011416, rel=1e-5 )
    assert float( pairs[0]['mag_fit'] ) == pytest.approx( -5.643705763605659, abs=1e-4 )
    assert float( pairs[0]['mag_fit_err'] ) == pytest.approx( 0.15218841286411408, abs=1e-4 )

    assert float( pairs[1]['aperture_sum'] ) == pytest.approx( 4111.8253286559275, rel=1e-5 )
    assert float( pairs[1]['flux_fit'] ) == pytest.approx( 719.6872116027234, rel=1e-5 )
    assert float( pairs[1]['flux_fit_err'] ) == pytest.approx( 27.85733758886825, rel=1e-5 )
    assert float( pairs[1]['mag_fit'] ) == pytest.approx( -7.1428594640283265, abs=1e-4 )
    assert float( pairs[1]['mag_fit_err'] ) == pytest.approx( 0.04202620180098437, abs=1e-4 )

    # TODO : cleanup output directories!  This is scary if you're using the same
    #   directories for tests and for running... so don't do that... but the
    #   way we're set up right now, you probably are.
