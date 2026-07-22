__all__ = [ 'sky_subtract', 'stampmaker' ]

# IMPORTS Standard:
from astropy.io import fits
from astropy.stats import SigmaClip
import numpy as np
import pathlib
from photutils.background import Background2D
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
import random
from scipy.signal import convolve2d
from scipy.interpolate import griddata
from roman_datamodels import dqflags

# IMPORTS SFFT:
from sfft.utils.StampGenerator import Stamp_Generator

# IMPORTS internal
from snappl.config import Config
import snappl.image
from snappl.logger import SNLogger


def interpolate_over_bad_pixels(data, dataflags, fill_value=0):
    """Interpolate over bad pixels in an image.

    NOTE: This is stolen from sidecar, and is written by WMWV. Lauren modified it.
    Parameters
    ----------
    data : ndarray
        Input data array.
    dataflags : ndarray
        Input flags array.
    bad_pixel_flags : int, optional
        Bitwise combination of DQ flags used to identify bad pixels.
    fill_value : float, optional
        Value to use for filling bad pixels if interpolation is not possible. Default is 0.

    Returns
    -------
    interpolated_data : ndarray
        The input image data with bad pixels interpolated over.

    bad_pixel_mask : ndarray
        Boolean mask indicating the locations of bad pixels that were interpolated over.

    Notes
    -----
    This function identifies bad pixels based on the provided `bad_pixel_flags`
    and performs interpolation to fill in those pixels.
    The interpolation method can be simple (e.g., nearest neighbor) or more sophisticated
    (e.g., using neighboring pixel values), depending on the implementation.
    """
    # Create a mask of bad pixels based on the flags
    bad_pixel_flags = 2 ** len(dqflags.pixel) - 1 - dqflags.pixel.WARM - dqflags.pixel.LOW_QE
    bad_mask = np.isnan(data) | (dataflags & bad_pixel_flags > 0)
    good_pixels = ~bad_mask

    interpolated_data = data.copy()
    interpolated_data[bad_mask] = np.nan  # Mark bad pixels as NaN for interpolation

    # Example interpolation using nearest neighbor
    x, y = np.indices(data.shape)
    interpolated_data[bad_mask] = griddata(
        (x[good_pixels], y[good_pixels]), data[good_pixels], (x[bad_mask], y[bad_mask]), method="nearest"
    )

    return interpolated_data, bad_mask


def sky_subtract( img, temp_dir=None,
                  nonlinear_threshold=1000,
                  footprint_radius=10,
                  mask_radius=5,
                  **kwargs ):
    """Subtracts background using photutils.

    Parameters
    ----------
      img: snappl.image.Image
        Original image.

      temp_dir: Path, default None
        Already-existing directory where we can write a temporary file.
        Defaults to photometry.snappl.temp_dir from the config.

    Returns
    -------
      skysubim: snappl.image.FITSImage, detmask: snappl.image.FITSImage, skyrms: float

         skysubim is the sky-subtracted image.  Its location on disk
         will be underneath temp_dir. It's the caller's responsibility
         to clean this up.  The file will have been written, so you
         can pass skysubim.path to any thing that needs the path of a
         single-HDU FITS image.

         detmask is the detection mask.  Its location on disk will be
         underneath temp_dir.  It's the caller's responsibility to clean
         this up.  The file will have been written, so you can pass
         detmask.path to any thing that needs the path of a single-HDU
         FITS image.

         skyrms is the median of the skyrms image calculated by source-extractor

    """

    # SEx_SkySubtract.SSS requires FITS files to chew on.  At some point
    # we should refactor this so that we can pass data to it.  However,
    # for now, write a snappl.image.FITSImage so we have something
    # to give to it.
    temp_dir = pathlib.Path( temp_dir if temp_dir is not None else Config.get().value( 'photometry.snappl.temp_dir' ) )
    barf = "".join( random.choices( "0123456789abcdef", k=10 ) )
    tmpfitspath = temp_dir / f"{barf}_fits_from_asdf.fits"
    tmpimpath = temp_dir / f"{barf}_sub.fits"
    tmpsubpath = temp_dir / f"{barf}_sub.fits"
    tmpdetmaskpath = temp_dir / f"{barf}_detmask.fits"
    origimg = img
    if isinstance( origimg, snappl.image.CompressedFITSImage ):
        img = origimg.uncompressed_version( include=['data'] )
        hdr = img.get_fits_header()
    elif isinstance( origimg, snappl.image.RomanDatamodelImage ):
        # The next few lines that make the header are stolen from sidecar.
        hdr = origimg.get_wcs().get_astropy_wcs().to_header(relax=True)
        hdr.insert(0, ("NAXIS", 2))
        hdr.insert("NAXIS", ("NAXIS1", img.data.shape[1]), after=True)
        hdr.insert("NAXIS1", ("NAXIS2", img.data.shape[0]), after=True)
        img = snappl.image.FITSImage(
                                      full_filepath=tmpfitspath,
                                      data=origimg.get_data(which='data')[0],
                                      header=hdr
                                    )
        img.save()
    else:
        # TODO, MAYBE MAKE THIS BETTER WHEN SNAPPL SUPPORTS MORE THINGS
        # We need to exract just the image data (not the noise or flags)
        # to send to image subtraction.
        # Lauren, make an issue about this, mauybe also a snappl image
        # that says that we need a way of making imgaes from other
        # images including only data... compressedfitsimage supports
        # that right now but not fitsimage).
        # Can take out header arg when snappl issue #77 is resolved.
        img = snappl.image.FITSImage( path=tmpimpath, header=fits.header.Header() )
        img.data = origimg.data
        hdr = img.get_fits_header()
        img.save( which='data' )

    SNLogger.debug( "Interpolate over bad pixels...")
    # NOTE: Make interp_mask do something at a later time.
    interp_data, _ = interpolate_over_bad_pixels(img.data, origimg.flags)

    SNLogger.debug( "Beginning sky subtraction..." )
    # We should talk about where radius_cut_detmask should go now that I have gotten
    # rid of Source Extractor.
    # radius_cut_detmask = Config.get().value( 'photometry.phrosty.sfft.radius_cut_detmask' )

    bkg = Background2D(interp_data, box_size=64)

    sky_subtracted_data = interp_data - bkg.background
    rms = bkg.background_rms_median

    # Based on the photutils.background documentation
    sigma_clip = SigmaClip(sigma=2.0, maxiters=10)
    threshold = detect_threshold(sky_subtracted_data, n_sigma=20.0, sigma_clip=sigma_clip)

    # Build a mask of pixels that are in the non-linear retime
    mask = np.abs(sky_subtracted_data) > nonlinear_threshold
    # Grow individual pixels by mask_radius
    mask_footprint = circular_footprint(radius=mask_radius)
    mask = convolve2d(mask, mask_footprint, fillvalue=0, mode="same")

    segment_img = detect_sources(sky_subtracted_data, threshold, n_pixels=10, mask=mask)
    detection_footprint = circular_footprint(radius=10)

    # convert boolean into float 1, and 0 because data must be float (not bool or int).
    detmask_data = np.asarray(segment_img.make_source_mask(footprint=detection_footprint), dtype="float")

    SNLogger.debug( "...back from sky subtraction." )

    subim = snappl.image.FITSImage(
                                    full_filepath=tmpsubpath,
                                    data=sky_subtracted_data,
                                    header=hdr
                                  )
    subim.save()

    detmaskim = snappl.image.FITSImage(
                                        full_filepath=tmpdetmaskpath,
                                        data=detmask_data,
                                        header=hdr
                                      )
    detmaskim.save()

    return subim, detmaskim, rms


def stampmaker(ra, dec, shape, img, savedir=None, savename=None, data_prop='data'):
    """Make stamps.

    TODO : pass an array of ra and dec to make this more efficient;
    otherwise, we may be writing and deleting a FITS image over and over
    again!


    Parameters
    ----------
      ra: float
        RA of center of stamp in degrees.

      dec: float
        Dec of center of stamp in degrees.

      shape: np.array
        Shape of stamp. must be a numpy array. e.g. np.array([100,100])

      img: snappl.image.Image
        Image from whose data the stamp will be extracted.

      savedir: Path, default None
        Directory stamp will be saved to.  Defaults to "stamps" underneath
        phtometry.phrosty.paths.dia_out_dir from the config.

      savename: Path, default None
        Base name of stamp filepath; defaults to the name from img's path.

    Returns
    -------
      savepath: Path
        Full savepath of stamp. Combined inputs "savedir/savename".

    """

    if savedir is None:
        cfg = Config.get()
        savedir = pathlib.Path( cfg.value( 'photometry.phrosty.paths.dia_out_dir' ) ) / "stamps"
    else:
        savedir = pathlib.Path( savedir )
    savedir.mkdir( parents=True, exist_ok=True )

    savepath = savedir / ( savename if savename is not None else f'stamp_{img.path.stem}.fits' )

    x, y = img.get_wcs().world_to_pixel( ra, dec )
    pxradec = np.array([[x, y]])

    # Give Stamp_Generator.SG a FITS image to chew on
    barf = "".join( random.choices( "0123456789abcdef", k=10 ) )
    tmpfitspath = savedir / f"{barf}.fits"
    origimg = img
    try:
        if isinstance( origimg, snappl.image.CompressedFITSImage ):
            img = origimg.uncompressed_version( include=[data_prop] )
        elif isinstance( origimg, snappl.image.RomanDatamodelImage ):
            # The next few lines that make the header are stolen from sidecar.
            hdr = origimg.get_wcs().get_astropy_wcs().to_header(relax=True)
            hdr.insert(0, ("NAXIS", 2))
            hdr.insert("NAXIS", ("NAXIS1", img.data.shape[1]), after=True)
            hdr.insert("NAXIS1", ("NAXIS2", img.data.shape[0]), after=True)
            img = snappl.image.FITSImage(
                                          full_filepath=tmpfitspath,
                                          data=origimg.get_data(which='data')[0],
                                          header=hdr
                                        )
            img.save()

        else:
          # NOTE : this next line will break if not using an image type that
          #   can return a FITS header!  The real solution is to fix SFFT
          #   so that it's not dependent on FITS images; just pass what's
          #   needed to Stamp_Generator.SG instead of assuming it will read
          #   all the right things out of the header.
          # See issue 177: https://github.com/Roman-Supernova-PIT/phrosty/issues/177
            img = snappl.image.FITSImage( path=savedir / f"{barf}.fits", header=origimg.get_fits_header() )
            img.data = origimg.data
            img.save( which='data' )

        # TODO : if Stamp_Generator.SG can take a Path in FITS_StpLst, remove the str()
        Stamp_Generator.SG(FITS_obj=img.path, COORD=pxradec, COORD_TYPE='IMAGE',
                           STAMP_IMGSIZE=shape, FILL_VALUE=np.nan, FITS_StpLst=str(savepath))

        return savepath

    finally:
        # Clean up temp file if necessary
        if img.path != origimg.path:
            img.path.unlink( missing_ok=True )
