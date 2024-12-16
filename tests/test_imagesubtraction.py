import os
import pathlib
import pytest
import time

from astropy.io import fits
from astropy.wcs import WCS

import phrosty
import phrosty.imagesubtraction

import numpy as np
import numpy.random as random

def test_gz_and_ext( dia_out_dir, compressed_template_image_path ):
    out_path = dia_out_dir / "gunzipped.fits"
    try:
        assert out_path == phrosty.imagesubtraction.gz_and_ext( compressed_template_image_path, out_path )

        with fits.open( out_path ) as hdu:
            assert len(hdu) == 1
            assert hdu[0].data.shape == (4088, 4088)

    finally:
        out_path.unlink( missing_ok=True )


def test_sky_subtract( dia_out_dir ):
    in_path = dia_out_dir / "random.fits"
    skysub_path = dia_out_dir / "skysub.fits"
    detmask_path = dia_out_dir / "detmask.fits"

    try:
        rng = random.default_rng( 42 )
        imdata = rng.normal( 100., 10., ( 512, 512 ) )
        fits.writeto( in_path, imdata )

        skymedrms = phrosty.imagesubtraction.sky_subtract( in_path, skysub_path, detmask_path, temp_dir=dia_out_dir )
        # Surprised it's not closer to 10 given the number of pixels, but whatevs
        assert skymedrms == pytest.approx( 10., abs=0.2 )
        with fits.open( skysub_path ) as hdul:
            data = hdul[0].data
            assert data.mean() == pytest.approx( 0., abs=5 * 10. / 512. )   # 3σ
            assert data.std() == pytest.approx( 10., rel=0.05 )
            assert hdul[0].header['SKYRMS'] == pytest.approx( 10., abs=0.2 )

    finally:
        for f in ( in_path, skysub_path, detmask_path ):
            f.unlink( missing_ok=True )


def _check_resampled_image( templ, resamp ):
    with fits.open( templ ) as t, fits.open( resamp ) as r:
        # Do some quick checks to see that values aren't absurd.
        # Visually picked a spot on the image that has a cute
        # little galaxy in it.  (It should make you feel very small
        # to realize that an astronomer at z=0.4 just called the Milky
        # Way a cute little galaxy.)
        tx = [ 2495, 2495, 2564, 2564 ]
        ty = [ 3148, 3185, 3148, 3185 ]
        twcs = WCS( t[0].header )
        sc = twcs.pixel_to_world( tx, ty )
        rwcs = WCS( r[0].header )
        rx, ry = rwcs.world_to_pixel( sc )

        # They're not perfectly aligned rectangles, so try to chop out
        #   something that's close
        rx.sort()
        ry.sort()
        rxmin = int( ( rx[0] + rx[1] ) / 2. )
        rxmax = rxmin + tx[2] - tx[0]
        rymin = int( ( ry[0] + ry[1] ) / 2. )
        rymax = rymin + ty[1] - ty[0]

        tdata = t[0].data[ ty[0]:ty[1] , tx[0]:tx[2] ]
        rdata = r[0].data[ rymin:rymax , rxmin:rxmax ]
        assert tdata.shape == rdata.shape

        assert np.median( tdata ) == pytest.approx( np.median( rdata ), rel=0.02 )
        assert tdata.std() == pytest.approx( rdata.std(), rel=0.04 )

def test_run_resample( dia_out_dir, template_image_path, one_science_image_path ):
    sci = one_science_image_path
    templ = template_image_path
    resamp = dia_out_dir / "resamp.fits"

    try:
        phrosty.imagesubtraction.run_resample( templ, sci, resamp )
        _check_resampled_image( templ, resamp )
    finally:
        resamp.unlink( missing_ok=True )


# def test_imalign(  dia_out_dir, template_image_path, one_science_image_path ):
#     sci = dia_out_dir / "sci.fits"
#     templ = dia_out_dir / "templ.fits"
#     resamp = dia_out_dir / "align/resamp.fits"
#     phrosty.imagesubtraction.gz_and_ext( one_science_image_path, sci )
#     phrosty.imagesubtraction.gz_and_ext( template_image_path, templ )

#     # BROKEN, FIX PATH PASSING

#     try:
#         t0 = time.perf_counter()
#         import pdb; pdb.set_trace()
#         phrosty.imagesubtraction.imalign( templ, sci, out_path=resamp )
#         firsttime = time.perf_counter() - t0
#         _check_resampled_image( templ, resamp )

#         # Rerun again with force false, it should be faster
#         # This is a bit scary, because we're depending on the
#         # filesystem being fast, and that may be bad, because
#         # GPUs are fast too.
#         t0 = time.perf_counter()
#         phrosty.imagesubtraction.imalign( templ, sci, out_path=resamp, force=False )
#         forcefalsetime = time.perf_counter() - t0
#         assert forcefalsetime < firsttime / 10.
#         _check_resampled_image( templ, resamp )

#         # Now run with force, it should be slow again
#         t0 = time.perf_counter()
#         phrosty.imagesubtraction.imalign( templ, sci, out_path=resamp, force=True )
#         forcetruetime = time.perf_counter() - t0
#         assert forcetruetime / firsttime == pytest.approx( 1., rel=0.2 )
#         _check_resampled_image( templ, resamp )

#     finally:
#         resamp.unlink( missing_ok=True )
#         sci.unlink( missing_ok=True )
#         templ.unlink( missing_ok=True )


def test_get_imsim_psf( sims_dir, sn_info_dir, dia_out_dir,  test_dia_image, test_sn, one_science_image_path ):
    impath = one_science_image_path
    ra = test_sn[ 'ra' ]
    dec = test_sn[ 'dec' ]
    band = test_dia_image[ 'band' ]
    pointing = test_dia_image[ 'pointing' ]
    sca = test_dia_image[ 'sca' ]
    config_yaml = sn_info_dir / "tds.yaml"
    psf_path = dia_out_dir / "psf.fits"
    size = 201

    phrosty.imagesubtraction.get_imsim_psf( impath, ra, dec, band, pointing, sca, size=size,
                                            config_yaml_file=config_yaml, psf_path=psf_path )
    with fits.open( psf_path ) as psf:
        ctr = size // 2
        # Center pixel should be way brighter, it's undersampled
        for dx, dy in zip( [ -1, 1, 0, 0 ], [ 0, 0, -1, 1 ] ):
            assert psf[0].data[ ctr, ctr ] > 10.* ( psf[0].data[ ctr+dx, ctr+dy ] )
            # OMG the PSF isn't normalized, I hope we deal with this right
            assert psf[0].data.sum() == pytest.approx( 2.33, abs=0.01 )

    # TODO test force and all that


def test_get_imsim_psf_photonOps( sims_dir, sn_info_dir, dia_out_dir,
                                  test_dia_image, test_sn, one_science_image_path ):
    impath = one_science_image_path
    ra = test_sn[ 'ra' ]
    dec = test_sn[ 'dec' ]
    band = test_dia_image[ 'band' ]
    pointing = test_dia_image[ 'pointing' ]
    sca = test_dia_image[ 'sca' ]
    config_yaml = sn_info_dir / "tds.yaml"
    psf_path = dia_out_dir / "psf.fits"
    size = 201
    photonOps = True
    n_phot = 1e6
    oversampling_factor = 1

    phrosty.imagesubtraction.get_imsim_psf( impath, ra, dec, band, pointing, sca, size=size,
                                            config_yaml_file=config_yaml, psf_path=psf_path,
                                            include_photonOps=photonOps, n_phot=n_phot,
                                            oversampling_factor=oversampling_factor )
    with fits.open( psf_path ) as psf:
        ctr = size // 2
        # Center pixel should be way brighter, it's undersampled
        # Photon shooting seems to blur out the PSF a bit, so the
        # condition is center pixel is >8× neighbors, whereas it
        # was 10x in the previous test.
        for dx, dy in zip( [ -1, 1, 0, 0 ], [ 0, 0, -1, 1 ] ):
            assert psf[0].data[ ctr, ctr ] > 3.* ( psf[0].data[ ctr+dx, ctr+dy ] )
            # ...but it looks like it's normalized now!
            assert psf[0].data.sum() == pytest.approx( 1.000, abs=0.002 )

    # TODO test force and all that


def test_stampmaker( dia_out_dir, test_dia_image, test_sn, one_science_image_path ):
    savedir = dia_out_dir
    savename = 'stamp.fits'
    ra = test_sn[ 'ra' ]
    dec = test_sn[ 'dec' ]
    shape = np.array( [ 100, 100 ] )

    savepath = None
    try:
        savepath = pathlib.Path( phrosty.imagesubtraction.stampmaker( ra, dec, shape, one_science_image_path,
                                                                      savedir=savedir, savename=savename ) )
        with fits.open( savepath ) as stamp:
            assert stamp[0].data.shape == tuple( shape )
            skylevel = np.median( stamp[0].data[ 6:40, 6:40 ] )
            skysig = stamp[0].data[ 6:40, 6:40 ].std()
            # Make sure skysig is what we know it is from having looked at it before
            assert skysig == pytest.approx( 11.5, abs=0.2 )
            # Look at the supernova, make sure our square 3x3 aperture is at least 10σ
            assert ( stamp[0].data[ 50:53, 50:53 ].sum() - 9 * skylevel ) > 10. * skysig * 3.
    finally:
        if savepath is not None:
            savepath.unlink( missing_ok = True )



