import numpy as np
import pathlib
import pytest
from matplotlib import pyplot
from astropy.table import Table

import astropy.units as u


from phrosty.pipeline import Pipeline
from phrosty.utils import map_band
from snappl.diaobject import DiaObject
from snappl.imagecollection import ImageCollection
from snappl.lightcurve import Lightcurve
from snappl.image import FITSImageStdHeaders


# TODO : separate tests for PipelineImage, for all the functions in
#   PipelineImage and Pipeline.

# This one writes a diagnostic plot file to test_plots/test_pipeline_run_simple_gauss1.pdf
def test_pipeline_run_simple_gauss1( config ):
    obj = DiaObject.find_objects( collection='manual', name='foo', ra=120, dec=-13. )[0]
    imgcol = ImageCollection.get_collection( 'manual_fits', subset='threefile',
                                             base_path='../photometry_test_data/simple_gaussian_test/sig2.0' )

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

    try:
        config._static = False
        config.set_value( 'photometry.phrosty.psf.type', 'gaussian' )
        # config.set_value( 'photometry.phrosty.psf.params', { 'sigmax': 1., 'sigmay': 1., 'theta': 0. } )
        config.set_value( 'photometry.phrosty.psf.params', { 'sigmax': 2., 'sigmay': 2., 'theta': 0. } )

        pip = Pipeline( obj,
                        imgcol,
                        'R062',
                        science_images=sciim,
                        template_images=tmplim,
                        catchfailures=False,
                        backend='numpy' )
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
        lc_obj = Lightcurve( filepath=ltcv, no_base_path=True )
        for row in lc_obj.lightcurve:
            mjd = row['mjd']
            # We know what the fluxes are supposed to be; see
            #   photometry_test_data/simple_gaussian_test/sig1.0/README.md

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
        config.set_value( 'photometry.phrosty.psf.type', orig_psf )
        config._static = True


@pytest.mark.parametrize( 'nprocs,nwrite', [(1, 1), (3, 3)] )
def test_pipeline_run( config, object_for_tests, ou2024_image_collection,
                       one_ou2024_template_image, two_ou2024_science_images,
                       nprocs, nwrite ):

    config._static = False
    orig_nprocs = config.value( 'photometry.phrosty.nprocs' )
    orig_nwrite = config.value( 'photometry.phrosty.nwrite' )

    try:
        config.set_value( 'photometry.phrosty.nprocs', nprocs )
        config.set_value( 'photometry.phrosty.nwrite', nwrite )
        pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                        science_images=two_ou2024_science_images,
                        template_images=one_ou2024_template_image,
                        backend='numpy' )
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
            # NOTE: THE ZEROPOINT CHECK IS COMMENTED OUT UNTIL THE SNAPPL ZEROPOINT STUFF
            # IS MORE COMPLETE. We can't really check zeropoint at this time because it
            # doesn't make sense to check it if we aren't sure how we are handling
            # zeropoints yet.
            # assert float(pair['zpt']) == pytest.approx( 32.6617, abs=0.0001 )

        # Tests aren't exactly reproducible from one run to the next,
        #   because the OU24 classes have random numbers in them, and we aren't
        #   controlling the seed.  So, we can only test for approximately
        #   consistent results.  Going to do 0.3 times the uncertainty,
        #   because a difference by that much is not all that meaningful
        #   change, and empirically they vary by that much.  (Which is
        #   alarming, but what can you do.)

        dflux = float( pairs[0]['flux_err'] )
        assert dflux == pytest.approx( 343., rel=0.3 )
        dmag = float( pairs[0]['mag_err'] )
        assert dmag == pytest.approx( 0.5, abs=0.1 )
        assert float( pairs[0]['aperture_sum'] ) == pytest.approx( 880., abs=0.3*dflux )
        assert float( pairs[0]['flux'] ) == pytest.approx( 708., abs=0.3*dflux )
        assert float( pairs[0]['mag'] ) == pytest.approx( -7.1, abs=max( 0.3*dmag, 0.01 ) )

        dflux = float( pairs[1]['flux_err'] )
        assert dflux == pytest.approx( 335., rel=0.3 )
        dmag = float( pairs[1]['mag_err'] )
        assert dmag == pytest.approx( 0.11, abs=0.1 )
        assert float( pairs[1]['aperture_sum'] ) == pytest.approx( 4010., abs=0.3*dflux )
        assert float( pairs[1]['flux'] ) == pytest.approx( 3028., abs=0.3*dflux )
        assert float( pairs[1]['mag'] ) == pytest.approx( -8.7, abs=max( 0.3*dmag, 0.01 ) )

        # TODO : cleanup output directories!  This is scary if you're using the same
        #   directories for tests and for running... so don't do that... but the
        #   way we're set up right now, you probably are.

    finally:
        config.set_value( 'photometry.phrosty.nprocs', orig_nprocs )
        config.set_value( 'photometry.phrosty.nwrite', orig_nwrite )
        config._static = True


@pytest.mark.parametrize( 'nprocs,nwrite', [(1, 1), (3, 3)] )
def test_no_failures( config, object_for_tests, ou2024_image_collection,
                      one_ou2024_template_image, two_ou2024_science_images,
                      nprocs, nwrite ):

    # This test makes sure the pipeline does not flag any failures for images that are
    # supposed to work fine.
    # TODO: Expand beyond OU2024 images.


    # We are, once again, going to mess with the config.
    orig_nprocs = config.value( 'photometry.phrosty.nprocs' )
    orig_nwrite = config.value( 'photometry.phrosty.nwrite' )
    config._static = False
    try:
        config.set_value( 'photometry.phrosty.nprocs', nprocs )
        config.set_value( 'photometry.phrosty.nwrite', nwrite )
        pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                        science_images=two_ou2024_science_images[0],
                        template_images=one_ou2024_template_image,
                        catchfailures=True,
                        backend='numpy' )

        pip()

        for key in pip.failures:
            assert len(pip.failures[key]) == 0

    finally:
        config.set_value( 'photometry.phrosty.nprocs', orig_nprocs )
        config.set_value( 'photometry.phrosty.nwrite', orig_nwrite )
        config._static = True


@pytest.mark.parametrize( 'nprocs,nwrite', [(1, 1), (3, 3)] )
def test_psf_retrieval_failures( config, object_for_tests, ou2024_image_collection,
                                 one_ou2024_template_image, two_ou2024_science_images,
                                 nprocs, nwrite ):
    # This is more of a test of the failure-flagging feature than the PSF retrieval.
    # We do not test the gaussian PSF because it is not dependent on any metadata
    # about the images.
    #
    # Testing for the following PSF types is broken because:
    # ou24PSF: There is a line in the tds.yaml file that specifies the path:
    # input:
    #   obseq_data:
    #     file_name: /photometry_test_data/ou2024/Roman_TDS_obseq_11_6_23.fits
    # We just discussed deprecating support for OU24 images. Is this worth fixing?
    # A25ePSF: If OU24 is deprecated, then so are these because they were made
    # on OU24 images.

    # possible_psfs = ['ou24PSF', 'A25ePSF', 'STPSF']
    possible_psfs = ['STPSF']

    # Save the original PSF value so we can mess with it and set it back later...
    config._static = False
    orig_psf = config.value( 'photometry.phrosty.psf.type' )
    orig_nprocs = config.value( 'photometry.phrosty.nprocs' )
    orig_nwrite = config.value( 'photometry.phrosty.nwrite' )

    scratchdir = pathlib.Path( config.value( 'system.paths.temp_dir' ) )
    test_image_path = scratchdir / 'test_fake_band'
    test_image = FITSImageStdHeaders( full_filepath=test_image_path,
                                      data=two_ou2024_science_images[0].data,
                                      noise=two_ou2024_science_images[0].noise,
                                      flags=np.zeros(two_ou2024_science_images[0].image_shape, dtype=np.int32),
                                      std_imagenames=True
                                    )

    test_image._header = two_ou2024_science_images[0].get_fits_header().copy()
    test_image._mjd = two_ou2024_science_images[0].mjd
    test_image._observation_id = two_ou2024_science_images[0].observation_id
    test_image._sca = two_ou2024_science_images[0].sca
    test_image.band = 'horsey' # Fake band value so it can't find the PSF.

    test_image.save(overwrite=True)

    try:
        config.set_value( 'photometry.phrosty.nprocs', nprocs )
        config.set_value( 'photometry.phrosty.nwrite', nwrite )
        for psftype in possible_psfs:
            print('psftype', psftype)
            config.set_value( 'photometry.phrosty.psf.type', psftype )
            pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                            science_images=[test_image, two_ou2024_science_images[1]],
                            template_images=[one_ou2024_template_image],
                            catchfailures=True, backend='numpy' )

            pip()

            assert len(pip.failures['skysub']) == 0
            assert len(pip.failures['get_psf']) == 1
            assert len(pip.failures['align_and_preconvolve']) == 0
            assert len(pip.failures['find_decorrelation']) == 0
            assert len(pip.failures['subtract']) == 0
            assert len(pip.failures['variance']) == 0
            assert len(pip.failures['apply_decorrelation']) == 0
            assert len(pip.failures['make_stamps']) == 0

    finally:
        # Put the config back...
        config.set_value( 'photometry.phrosty.psf.type', orig_psf )
        config.set_value( 'photometry.phrosty.nprocs', orig_nprocs )
        config.set_value( 'photometry.phrosty.nwrite', orig_nwrite )
        config._static = True


@pytest.mark.parametrize( 'nprocs,nwrite', [(1, 1), (3, 3)] )
def test_nan_handling( config, object_for_tests, ou2024_image_collection,
                       one_ou2024_template_image, two_ou2024_science_images,
                       nprocs, nwrite ):

    # Make sure that inputting a science or template image full of NaNs causes
    # the pipeline to fail.

    # We will mess with the config again.
    config._static = False
    orig_nprocs = config.value( 'photometry.phrosty.nprocs' )
    orig_nwrite = config.value( 'photometry.phrosty.nwrite' )

    # Make an image full of NaN:
    scratchdir = pathlib.Path( config.value( 'system.paths.temp_dir' ) )
    test_image_path = scratchdir / 'test_nan'
    nan_image = FITSImageStdHeaders( full_filepath=test_image_path,
                                     data=np.full(two_ou2024_science_images[0].image_shape, np.nan),
                                     noise=np.full(two_ou2024_science_images[0].image_shape, np.nan),
                                     flags=np.zeros(two_ou2024_science_images[0].image_shape, dtype=np.int32),
                                     std_imagenames=True
                                   )

    nan_image._wcs = two_ou2024_science_images[0].get_wcs()
    nan_image.band = 'Y106'
    nan_image.observation_id = 35198
    nan_image.sca = 2
    nan_image.mjd = two_ou2024_science_images[0].mjd
    nan_image.save(
                    which='all',
                    overwrite=True
                  )
    try:
        config.set_value( 'photometry.phrosty.nprocs', nprocs )
        config.set_value( 'photometry.phrosty.nwrite', nwrite )
        # First, do the science image.
        pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                        science_images=[nan_image, two_ou2024_science_images[1]],
                        template_images=[one_ou2024_template_image],
                        catchfailures=True, backend='numpy' )


        ltcv = pip()

        assert len(pip.failures['skysub']) == 1
        assert len(pip.failures['get_psf']) == 0
        assert len(pip.failures['align_and_preconvolve']) == 0
        assert len(pip.failures['find_decorrelation']) == 0
        assert len(pip.failures['subtract']) == 0
        assert len(pip.failures['variance']) == 0
        assert len(pip.failures['apply_decorrelation']) == 0
        assert len(pip.failures['make_stamps']) == 0

        lc_obj = Lightcurve( filepath=ltcv, no_base_path=True )
        lc = lc_obj.lightcurve
        nanrow = lc[lc['science_name'] == nan_image.name]
        goodrow = lc[lc['science_name'] == two_ou2024_science_images[1].name]
        assert len(lc) == 2
        assert np.isnan(nanrow['flux'][0].value)
        assert not np.isnan(goodrow['flux'][0].value)

        # Now do the template image.
        pip = Pipeline( object_for_tests, ou2024_image_collection, 'Y106',
                        science_images=two_ou2024_science_images[1],
                        template_images=nan_image,
                        catchfailures=True,
                        backend='numpy' )

        ltcv = pip()

        assert len(pip.failures['skysub']) == 1
        assert len(pip.failures['get_psf']) == 0
        assert len(pip.failures['align_and_preconvolve']) == 0
        assert len(pip.failures['find_decorrelation']) == 0
        assert len(pip.failures['subtract']) == 0
        assert len(pip.failures['variance']) == 0
        assert len(pip.failures['apply_decorrelation']) == 0
        assert len(pip.failures['make_stamps']) == 0

        lc_obj = Lightcurve( filepath=ltcv, no_base_path=True )
        lc = lc_obj.lightcurve
        assert len(lc) == 1
        assert np.isnan(lc[0]['flux'].value)
        assert not lc[0]['success']
    finally:
        # Put the config back...
        config.set_value( 'photometry.phrosty.nprocs', orig_nprocs )
        config.set_value( 'photometry.phrosty.nwrite', orig_nwrite )
        config._static = True


def test_map_band():
    # Just check the band mapping function, solo.

    assert map_band('F062') == 'F062'
    assert map_band('R062') == 'F062'
    assert map_band('F087') == 'F087'
    assert map_band('Z087') == 'F087'
    assert map_band('F106') == 'F106'
    assert map_band('Y106') == 'F106'
    assert map_band('F129') == 'F129'
    assert map_band('J129') == 'F129'
    assert map_band('F158') == 'F158'
    assert map_band('H158') == 'F158'
    # F184 doesn't have a double because it's
    # always "F" in either convention.
    assert map_band('F184') == 'F184'
    assert map_band('F213') == 'F213'
    assert map_band('K213') == 'F213'

    with pytest.raises(ValueError, match='Unknown band'):
        map_band('horsey')


def test_band_mapping( config, object_for_tests, ou2024_image_collection,
                       template_csv, two_science_csv ):

    # This test makes sure the band mapping is working properly, and that
    # a fit is not thrown whether you input "F106" or "Y106".
    # NOTE: This only tests F106/Y106 because that's what is in
    # photometry_test_data.

    for band in ['Y106', 'F106']:
        print(band)
        pip = Pipeline( object_for_tests, ou2024_image_collection, band,
                        science_csv=two_science_csv, template_csv=template_csv,
                        catchfailures=True)

        assert pip.band == 'F106'
