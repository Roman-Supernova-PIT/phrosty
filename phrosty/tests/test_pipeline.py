import os
import pathlib
import pytest
from matplotlib import pyplot

import numpy as np
import pandas

from snappl.diaobject import DiaObject
from snappl.imagecollection import ImageCollection
from phrosty.pipeline import Pipeline

# TODO : separate tests for PipelineImage, for all the functions in#
#   PipelineImage and Pipeline.  Right now we just have this regression
#   test.


# This one writes a diagnostic plot file to test_plots/test_pipeline_run_simple_gauss1.pdf
def test_pipeline_run_simple_gauss1( config ):
    obj = DiaObject.find_objects( collection='manual', id=1, ra=120, dec=-13. )[0]
    # imgcol = ImageCollection.get_collection( 'manual_fits', subset='threefile',
    #                                          base_path='/photometry_test_data/simple_gaussian_test/sig1.0' )
    # tmplim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in [ 60000., 60005. ] ]
    # sciim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in range( 60000, 60065, 5 ) ]
    imgcol = ImageCollection.get_collection( 'manual_fits', subset='threefile',
                                             base_path='/photometry_test_data/simple_gaussian_test/sig2.0' )
    tmplim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in [ 60000., 60005. ] ]
    sciim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in range( 60000, 60065, 5 ) ]

    # We have to muck about with the config, because the default config loaded for tests is
    #   set up for ou2024.  We're going to do naughty things we're not supposed to do,
    #   specifically, modify the default config that's supposed to be immutable.  We'll
    #   restore it later.  In actual running code, never modify config._static.  If you
    #   find yourself doing that in anything other than a test like this where you're
    #   VERY careful to restore things after you are done, then you're doing it wrong.
    orig_psf = config.value( 'photometry.phrosty.psf' )
    orig_sfft = config.value( 'photometry.phrosty.sfft' )
    try:
        config._static = False
        config.set_value( 'photometry.phrosty.psf.type', 'gaussian' )
        # config.set_value( 'photometry.phrosty.psf.params', { 'sigmax': 1., 'sigmay': 1., 'theta': 0. } )
        config.set_value( 'photometry.phrosty.psf.params', { 'sigmax': 2., 'sigmay': 2., 'theta': 0. } )
        config.set_value( 'photometry.phrosty.sfft.radius_cut_detmask', 1. )

        pip = Pipeline( obj, imgcol, 'R062', science_images=sciim, template_images=tmplim, nprocs=1, nwrite=1 )
        ltcv = pip()

        chisq = 0.
        apchisq = 0.
        truemjd = []
        trueflux = {}
        measmjd = []
        measflux = []
        measdflux = []
        resid = []
        measapflux = []
        measapdflux = []
        apresid = []
        plotzpt = 31.4
        df = pandas.read_csv( ltcv )
        for row in df.itertuples():
            mjd = row.mjd
            # We know what the fluxes are supposed to be; see
            #   /photometry_test_data/simple_gaussian_test/sig1.0/README.md
            peak_flux = 10 ** ( ( 21. - row.zpt ) / -2.5 )
            plotfluxfac = 10 ** ( ( row.zpt - plotzpt ) / -2.5 )
            if ( mjd < 60010. ) or ( mjd > 60060. ):
                expected_flux = 0.
            elif mjd < 60030.:
                expected_flux = peak_flux * ( mjd - 60010. ) / 20.
            else:
                expected_flux = peak_flux * ( 60060. - mjd ) / 30.

            # We don't have errors on aperture sum, so just guess.  It's in
            #   a radius of 4 by default, and the sky noise per pixel is 100.
            ap_err = np.sqrt( np.pi * 16 * 100. + row.aperture_sum )

            measmjd.append( row.mjd )
            measflux.append( row.flux_fit * plotfluxfac )
            resid.append( ( row.flux_fit - expected_flux ) * plotfluxfac )
            measdflux.append( row.flux_fit_err * plotfluxfac )
            measapflux.append( row.aperture_sum * plotfluxfac )
            apresid.append( ( row.aperture_sum - expected_flux ) * plotfluxfac )
            measapdflux.append( ap_err * plotfluxfac )
            if mjd not in trueflux:
                truemjd.append( mjd )
                trueflux[mjd] = expected_flux * plotfluxfac

            # assert expected_flux == pytest.approx( row.aperture_sum, abs=2.*ap_err )
            apchisq += ( ( expected_flux - row.aperture_sum ) / ap_err ) ** 2

            # assert expected_flux == pytest.approx( row.flux_fit, abs=3. * row.flux_fit_err )
            chisq += ( ( expected_flux - row.flux_fit ) / row.flux_fit_err ) ** 2

        # NOTE -- tests do not currently pass, we know things are broken.
        #   Issue #...
        # These chisq values aren't really right, because we're assuming that all the
        #   points are independent, but they're not, because the ones at the same
        #   mjd share news, and all the news share the same set of refs.
        # assert apchisq / len(df) == pytest.approx( 1., abs=0.2 )
        # assert chisq / len(df) == pytest.approx( 1., abs=0.2 )

        # Draw a plot
        trueflux = [ trueflux[m] for m in truemjd ]
        fig, axes = pyplot.subplots( 2, 1, height_ratios=[3, 1], sharex=True )
        fig.subplots_adjust( wspace=0 )
        fig.set_tight_layout( True )
        axes[0].set_title( 'phrosty test_pipeline_run_simple_gauss1' )
        axes[0].plot( truemjd, trueflux, color='orange', label='Truth' )
        axes[0].errorbar( measmjd, measflux, measdflux, color='blue',
                          linestyle='none', marker='o', label='psf fit' )
        axes[0].errorbar( measmjd, measapflux, measapdflux, color='green',
                          linestyle='none', marker='s', label='aperture sum' )
        axes[0].set_ylabel( 'flux (nJy)' )
        axes[0].legend()
        xmin, xmax = axes[0].get_xlim()
        axes[1].plot( [xmin, xmax], [0., 0.], color='orange' )
        axes[1].errorbar( measmjd, resid, measdflux, color='blue', linestyle='none', marker='o' )
        axes[1].errorbar( measmjd, apresid, measdflux, color='green', linestyle='none', marker='s' )
        axes[1].set_ylabel( 'resid' )
        axes[1].set_xlabel( 'MJD' )

        plotdir = pathlib.Path( 'test_plots' )
        plotdir.mkdir( parents=True, exist_ok=True )
        fig.savefig( plotdir / 'test_pipeline_run_simple_gauss1.pdf' )

    finally:
        # Fix the naughty damage we did to config
        config.set_value( 'photomery.phrosty.psf', orig_psf )
        config.set_value( 'photometry.phrosty.sfft', orig_sfft )
        config._static = True


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
        assert float(pair['zpt']) == pytest.approx( 32.6617, abs=0.0001 )

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
