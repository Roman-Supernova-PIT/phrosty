import os
import pathlib
import pytest

from astropy.io import fits
from astropy.wcs import WCS

import phrosty.imagesubtraction

import numpy as np
import numpy.random as random

def test_gz_and_ext( dia_out_dir, template_image_path ):
    out_path = dia_out_dir / "gunzipped.fits"
    try:
        assert out_path == phrosty.imagesubtraction.gz_and_ext( template_image_path, out_path )

        with fits.open( out_path ) as hdu:
            assert len(hdu) == 1
            assert hdu[0].data.shape == (4088, 4088)

    finally:
        out_path.unlink( missing_ok=True )


def test_sky_subtract( dia_out_dir, template_image_path ):
    in_path = dia_out_dir / "indata.fits"
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
            

def test_run_resample( dia_out_dir, template_image_path, one_science_image_path ):
    sci = dia_out_dir / "sci.fits"
    templ = dia_out_dir / "templ.fits"
    resamp = dia_out_dir / "resamp.fits"
    phrosty.imagesubtraction.gz_and_ext( one_science_image_path, sci )
    phrosty.imagesubtraction.gz_and_ext( template_image_path, templ )

    try:
        phrosty.imagesubtraction.run_resample( templ, sci, resamp )
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
    finally:
        resamp.unlink( missing_ok=True )
    
