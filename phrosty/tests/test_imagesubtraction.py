import os
import pathlib
import pytest

import numpy as np
import numpy.random as random

from astropy.io import fits
from astropy.wcs import WCS

from snappl.image import FITSImageOnDisk
import phrosty.imagesubtraction



@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0), reason="SKIP_GPU_TESTS is set")
def test_sky_subtract( dia_out_dir ):
    in_path = dia_out_dir / "random.fits"
    skysub_path = dia_out_dir / "skysub.fits"
    detmask_path = dia_out_dir / "detmask.fits"

    try:
        rng = random.default_rng( 42 )
        imdata = rng.normal( 100., 10., ( 512, 512 ) )
        hdr = fits.header.Header()
        fits.writeto( in_path, imdata, header=hdr )
        img = FITSImageOnDisk( path=in_path )

        subim, _detmask, skymedrms = phrosty.imagesubtraction.sky_subtract( img, temp_dir=dia_out_dir )
        assert skymedrms == pytest.approx( 10., abs=0.2 )
        assert subim.data.mean() == pytest.approx( 0., abs=5 * 10. / 512. )   # 3σ
        assert subim.data.std() == pytest.approx( 10., rel=0.05 )

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
        tx = [ 3370, 3370, 3404, 3404 ]
        ty = [  396,  422,  396,  422 ]
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

        assert np.median( tdata ) == pytest.approx( np.median( rdata ), rel=0.01 )
        assert np.sum( tdata ) == pytest.approx( np.sum( rdata ), rel=0.01 )
        # The standard deviation goes down quite a bit... correlated pixel errors from the warping doing that.
        assert tdata.std() == pytest.approx( rdata.std(), rel=0.5 )


@pytest.mark.skipif( os.getenv("SKIP_GPU_TESTS", 0), reason="SKIP_GPU_TESTS is set")
def test_stampmaker( object_for_tests, one_science_image, dia_out_dir ):
    savedir = dia_out_dir
    savename = 'stamp.fits'
    ra = object_for_tests.ra
    dec = object_for_tests.dec
    shape = np.array( [ 100, 100 ] )

    savepath = None
    try:
        savepath = pathlib.Path( phrosty.imagesubtraction.stampmaker( ra, dec, shape, one_science_image,
                                                                      savedir=savedir, savename=savename ) )
        with fits.open( savepath ) as stamp:
            assert stamp[0].data.shape == tuple( shape )
            skylevel = np.median( stamp[0].data[ 6:40, 6:40 ] )
            skysig = stamp[0].data[ 6:40, 6:40 ].std()
            # Make sure skysig is what we know it is from having looked at it before
            assert skysig == pytest.approx( 18.7, abs=0.1 )
            # Look at the supernova, make sure our square 3x3 aperture is at least 10σ
            assert ( stamp[0].data[ 50:53, 50:53 ].sum() - 9 * skylevel ) > 10. * skysig * 3.
    finally:
        if savepath is not None:
            savepath.unlink( missing_ok = True )
