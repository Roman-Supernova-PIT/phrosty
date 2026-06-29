import numpy as np
import os
import pathlib
import pytest
from matplotlib import pyplot
from astropy.table import Table
from pathlib import Path

import astropy.units as u


from phrosty.pipeline import Pipeline
from snappl.diaobject import DiaObject
from snappl.imagecollection import ImageCollection
from snappl.lightcurve import Lightcurve
from snappl.image import FITSImageStdHeaders
import snappl.psf


# TODO : separate tests for PipelineImage, for all the functions in#
#   PipelineImage and Pipeline.  Right now we just have this regression
#   test.

# This one writes a diagnostic plot file to test_plots/test_pipeline_run_simple_gauss1.pdf
@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0), reason="SKIP_GPU_TESTS is set")
def test_pipeline_run_simple_gauss1( config ):
    obj = DiaObject.find_objects( collection='manual', name='foo', ra=120, dec=-13. )[0]
    imgcol = ImageCollection.get_collection( 'manual_fits', subset='threefile',
                                             base_path='/photometry_test_data/simple_gaussian_test/sig2.0' )

    # Use for longer test with full "lightcurve" and two templates:
    tmplim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in [ 60000., 60005. ] ]
    sciim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in range( 60010, 60065, 5 ) ]

    # Use for shorter test with only two "observations":
    # tmplim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in [ 60000 ] ]
    # sciim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in [ 60030, 60035 ] ]

    # Shortest test with only one template and one science:
    # tmplim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in [ 60000 ] ]
    # sciim = [ imgcol.get_image(path=f'test_{t:7.1f}') for t in [ 60035 ] ]

    # We have to muck about with the config, because the default config loaded for tests is
    #   set up for ou2024.  We're going to do naughty things we're not supposed to do,
    #   specifically, modify the default config that's supposed to be immutable.  We'll
    #   restore it later.  In actual running code, never modify config._static.  If you
    #   find yourself doing that in anything other than a test like this where you're
    #   VERY careful to restore things after you are done, then you're doing it wrong.
    orig_psf = config.value( 'photometry.phrosty.psf.type' )
    orig_sfft = config.value( 'photometry.phrosty.sfft.radius_cut_detmask' )
    try:
        config._static = False
        config.set_value( 'photometry.phrosty.psf.type', 'gaussian' )
        # config.set_value( 'photometry.phrosty.psf.params', { 'sigmax': 1., 'sigmay': 1., 'theta': 0. } )
        config.set_value( 'photometry.phrosty.psf.params', { 'sigmax': 2., 'sigmay': 2., 'theta': 0. } )
        config.set_value( 'photometry.phrosty.sfft.radius_cut_detmask', 1. )

        pip = Pipeline( obj,
                        imgcol,
                        'R062',
                        science_images=sciim,
                        template_images=tmplim,
                        nprocs=1,
                        nwrite=1,
                        catchfailures=False )
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
        lc_obj = Lightcurve( filepath=ltcv )
        for row in lc_obj.lightcurve:
            mjd = row['mjd']
            # We know what the fluxes are supposed to be; see
            #   /photometry_test_data/simple_gaussian_test/sig1.0/README.md

            # LA: snappl enforces astropy units in the columns. I have tried to preserve
            # where it makes sense to have units vs. dimensionless quantities below when
            # deciding where to put .value, or where to add units back in.
            # Or, I did things to make the plot work.
            peak_flux = 10 ** ( ( 21. - row['zpt'].value ) / -2.5 ) * u.count / u.s
            plotfluxfac = 10 ** ( ( row['zpt'].value - plotzpt ) / -2.5 ) * u.count / u.s
            if ( mjd.value < 60010. ) or ( mjd.value > 60060. ):
                expected_flux = 0.
            elif mjd.value < 60030.:
                expected_flux = peak_flux * ( mjd.value - 60010. ) / 20.
            else:
                expected_flux = peak_flux * ( 60060. - mjd.value ) / 30.

            # We don't have errors on aperture sum, so just guess.  It's in
            #   a radius of 4 by default, and the sky noise per pixel is 100.
            ap_err = np.sqrt( np.pi * 16 * 100. + row['aperture_sum'] )

            # LA: I stripped the units off of everything below to make the plot work.
            measmjd.append( row['mjd'].value )
            measflux.append( row['flux'].value * plotfluxfac.value )
            resid.append( ( row['flux'].value - expected_flux.value ) * plotfluxfac.value )
            measdflux.append( row['flux_err'].value * plotfluxfac.value )
            measapflux.append( row['aperture_sum'] * plotfluxfac.value )
            apresid.append( ( row['aperture_sum'] - expected_flux.value ) * plotfluxfac.value )
            measapdflux.append( ap_err * plotfluxfac.value )
            if mjd not in trueflux:
                truemjd.append( mjd.value )
                trueflux[mjd.value] = expected_flux.value * plotfluxfac.value

            # assert expected_flux == pytest.approx( row['aperture_sum'], abs=2.*ap_err )
            apchisq += ( ( expected_flux.value - row['aperture_sum'] ) / ap_err ) ** 2

            # assert expected_flux == pytest.approx( row['flux'], abs=3. * row['flux_err'] )
            chisq += ( ( expected_flux.value - row['flux'].value ) / row['flux_err'] ) ** 2

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

        # WARNING HARDCODED LIMITS these should go away when things are
        #   fixed and the test works better.  Right now PSF photometry is
        #   going nuts, and sometimes produces infinite error bars
        # axes[0].set_ylim( -2000, 16000 )
        # axes[1].set_ylim( -1000, 1000 )

        plotdir = pathlib.Path( 'test_plots' )
        plotdir.mkdir( parents=True, exist_ok=True )
        fig.savefig( plotdir / 'test_pipeline_run_simple_gauss1.pdf' )

    finally:
        # Fix the naughty damage we did to config
        config.set_value( 'photomery.phrosty.psf.type', orig_psf )
        config.set_value( 'photometry.phrosty.sfft.radius_cut_detmask', orig_sfft )
        config._static = True


@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0 ), reason="SKIP_GPU_TESTS is set" )
def test_pipeline_run( object_for_tests, ou2024_image_collection,
                       one_ou2024_template_image, two_ou2024_science_images ):

    pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                    science_images=two_ou2024_science_images,
                    template_images=one_ou2024_template_image,
                    nprocs=1, nwrite=1 )
    ltcv = pip()

    ifp = Table.read(ltcv, format='parquet')
    hdrline = tuple(ifp.columns)
    assert hdrline == ( 'mjd', 'flux', 'flux_err', 'zpt', 'NEA', 'sky_rms', 'observation_id', 'sca',
                        'pix_x', 'pix_y', 'science_name', 'template_name', 'science_id', 'template_id',
                        'template_observation_id', 'template_sca', 'aperture_sum', 'mag', 'mag_err', 'success' )
    assert len(ifp) == 2

    pairs = []
    for row in ifp:
        pairs.append( { kw: val for kw, val in zip( hdrline, tuple(row) ) } )

    for pair, img in zip( pairs, two_ou2024_science_images ):
        assert float(pair['mjd']) == pytest.approx( img.mjd, abs=1e-3 )
        assert int(pair['observation_id']) == int(img.observation_id)
        assert int(pair['sca']) == int(img.sca)
        assert int(pair['template_observation_id']) == int(one_ou2024_template_image.observation_id)
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

    dflux = float( pairs[0]['flux_err'] )
    assert dflux == pytest.approx( 540., rel=0.3 )
    dmag = float( pairs[0]['mag_err'] )
    assert dmag == pytest.approx( 0.49, abs=0.1 )
    assert float( pairs[0]['aperture_sum'] ) == pytest.approx( 1006.897, abs=0.3*dflux )
    assert float( pairs[0]['flux'] ) == pytest.approx( 1193.509, abs=0.3*dflux )
    assert float( pairs[0]['mag'] ) == pytest.approx( -7.692, abs=max( 0.3*dmag, 0.01 ) )

    dflux = float( pairs[1]['flux_err'] )
    assert dflux == pytest.approx( 525.716, rel=0.3 )
    dmag = float( pairs[1]['mag_err'] )
    assert dmag == pytest.approx( 0.11, abs=0.1 )
    assert float( pairs[1]['aperture_sum'] ) == pytest.approx( 4050.392, abs=0.3*dflux )
    assert float( pairs[1]['flux'] ) == pytest.approx( 4986.560, abs=0.3*dflux )
    assert float( pairs[1]['mag'] ) == pytest.approx( -9.24, abs=max( 0.3*dmag, 0.01 ) )

    # TODO : cleanup output directories!  This is scary if you're using the same
    #   directories for tests and for running... so don't do that... but the
    #   way we're set up right now, you probably are.


@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0 ), reason="SKIP_GPU_TESTS is set" )
def test_no_failures( config, object_for_tests, ou2024_image_collection,
                      one_ou2024_template_image, two_ou2024_science_images ):
    
    # This test makes sure the pipeline does not flag any failures for images that are
    # supposed to work fine. 
    # TODO: Expand beyond OU2024 images. 
    
    nprocss = [1, 3]
    nwrites = [1, 3]

    for i in nprocss:
        for j in nwrites:
            pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                            science_images=two_ou2024_science_images[0],
                            template_images=one_ou2024_template_image,
                            nprocs=i, nwrite=j, catchfailures=True )

            lctv = pip()

            for key in pip.failures:
                assert len(pip.failures[key]) == 0

@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0 ), reason="SKIP_GPU_TESTS is set" )
def test_psf_retrieval_failures( config, object_for_tests, ou2024_image_collection, 
                                 one_ou2024_template_image, two_ou2024_science_images ):
    # This is more of a test of the failure-flagging feature than the PSF retrieval.

    possible_psfs = ['ou24PSF', 'A25ePSF', 'STPSF']

    # Save the original PSF value so we can mess with it and set it back later...
    config._static = False
    orig_psf = config.value( 'photometry.phrosty.psf.type' )

    test_image = FITSImageStdHeaders( full_filepath='/scratch/phrosty_temp/test_nan_img',
                                      data=np.full(two_ou2024_science_images[0].image_shape, np.nan),
                                      flags=np.zeros(two_ou2024_science_images[0].image_shape),
                                      std_imagenames=True
                                    )
    test_image.band = 'horsey' # Fake band value so it can't find the PSF. 

    nprocss = [1, 3]
    nwrites = [1, 3]
    for i in nprocss:
        for j in nwrites:
            for psftype in possible_psfs:
                print('nprocs', i)
                print('nwrites', j)
                print('psftype', psftype)
                config.set_value( 'photomery.phrosty.psf.type', psftype )
                pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                                science_images=[test_image, two_ou2024_science_images[1]],
                                template_images=[one_ou2024_template_image],
                                nprocs=i, nwrite=j, catchfailures=True )

                lctv = pip()

                for key in pip.failures:
                    print(key)
                    print(len(pip.failures[key]))

                assert len(pip.failures['skysub']) == 0 
                assert len(pip.failures['get_psf']) == 1
                assert len(pip.failures['align_and_preconvolve']) == 0
                assert len(pip.failures['find_decorrelation']) == 0
                assert len(pip.failures['subtract']) == 0
                assert len(pip.failures['variance']) == 0
                assert len(pip.failures['apply_decorrelation']) == 0
                assert len(pip.failures['make_stamps']) == 0


    # Put the config back...
    config.set_value( 'photomery.phrosty.psf.type', orig_psf )
    config._static = True

@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0 ), reason="SKIP_GPU_TESTS is set" )
def test_nan_handling( config, object_for_tests, ou2024_image_collection,
                               one_ou2024_template_image, two_ou2024_science_images ):

    # Inputting an image full of NaN apparently does not cause the pipeline to fail at all.
    # Fine, as long as the corresponding row in the parquet file is also full of NaNs. 

    nprocss = [1, 3]
    nwrites = [1, 3]

    # Make an image full of NaN:
    nanpaths = [Path('test_nan_img.fits'), Path('test_nan_noise.fits')]
    nan_image = FITSImageStdHeaders( full_filepath='/scratch/phrosty_temp/test_nan_image',
                                     data=np.full(two_ou2024_science_images[0].image_shape, np.nan),
                                     noise=np.full(two_ou2024_science_images[0].image_shape, np.nan),
                                     flags=np.zeros(two_ou2024_science_images[0].image_shape),
                                     std_imagenames=True
                                   )

    nan_image._wcs = two_ou2024_science_images[0].get_wcs()
    nan_image.band = 'Y106'
    nan_image.observation_id = 35198 
    nan_image.sca = 2
    nan_image.mjd = two_ou2024_science_images[0].mjd
    nan_image.save( 
                    imagepath=nanpaths[0],
                    noisepath=nanpaths[1],
                    which=['data','noise'],
                    overwrite=True
                  )

    for i in nprocss:
        for j in nwrites:
            print('############################################################')
            print('NEW LOOP')
            print('############################################################')
            print('nprocss', i)
            print('nwrites', j)
            # print("SCIENCE IMAGE IS NAN")
            # # The "science image is full of nans" case should fail because when
            # # it tries to do photometry on a stamp, it won't find a stamp at all.
            # pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
            #                 science_images=[nan_image, two_ou2024_science_images[1]],
            #                 template_images=[one_ou2024_template_image],
            #                 nprocs=i, nwrite=j, catchfailures=True )
            

            # ltcv = pip()

            # lc_obj = Lightcurve( filepath=ltcv )

            # for key in pip.failures.keys():
            #     print(key)
            #     print(pip.failures[key])

            # assert len(pip.failures['skysub']) == 0 
            # assert len(pip.failures['get_psf']) == 0
            # assert len(pip.failures['align_and_preconvolve']) == 0
            # assert len(pip.failures['find_decorrelation']) == 0
            # assert len(pip.failures['subtract']) == 0
            # assert len(pip.failures['variance']) == 0
            # assert len(pip.failures['apply_decorrelation']) == 0
            # assert len(pip.failures['make_stamps']) == 0

            # assert len(lc_obj.lightcurve) == 2

            # assert np.isnan(lc_obj.lightcurve['flux'][0].value)
            # assert not np.isnan(lc_obj.lightcurve['flux'][1].value)

            print('TEMPLATE IMAGE IS NAN')

            # LA: The "template image is full of nans" case does not currently result
            # in a row full of nan because it does find a stamp for the corresponding
            # difference image, which is full of 0. Thus, the PSF fitting results
            # in flux=0 (with nonzero flux error, which is fine because I have not
            # set the corresponding noise image to also be full of nans). I think this
            # is a nuance of how SFFT is implemented.
            pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                            science_images=two_ou2024_science_images[0],
                            template_images=nan_image,
                            nprocs=i, nwrite=j, catchfailures=True )

            ltcv = pip()
            lc_obj = Lightcurve( filepath=ltcv )

            print(lc_obj.lightcurve)

            assert len(pip.failures['skysub']) == 0 
            assert len(pip.failures['get_psf']) == 0
            assert len(pip.failures['align_and_preconvolve']) == 0
            assert len(pip.failures['find_decorrelation']) == 0
            assert len(pip.failures['subtract']) == 0
            assert len(pip.failures['variance']) == 0
            assert len(pip.failures['apply_decorrelation']) == 0
            assert len(pip.failures['make_stamps']) == 0

            assert len(lc_obj.lightcurve) == 1
            
            assert lc_obj.lightcurve['flux'][0].value == 0