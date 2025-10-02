__all__ = [ 'sky_subtract', 'stampmaker' ]

# IMPORTS Standard:
from astropy.io import fits
import numpy as np
import pathlib
import random

# IMPORTS SFFT:
from sfft.utils.SExSkySubtract import SEx_SkySubtract
from sfft.utils.StampGenerator import Stamp_Generator

# IMPORTS internal
import snappl.image
from snpit_utils.logger import SNLogger
from snpit_utils.config import Config


def sky_subtract( img, temp_dir=None ):
    """Subtracts background, found with Source Extractor.

    Parameters
    ----------
      img: snappl.image.Image
        Original image.

      temp_dir: Path, default None
        Already-existing directory where we can write a temporary file.
        Defaults to photometry.snappl.temp_dir from the config.

    Returns
    -------
      skysubim: snappl.image.FITSImageOnDisk, detmask: snappl.image.FITSImageOnDisk, skyrms: float

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
    # for now, write a snappl.image.FITSImageOnDisk so we have something
    # to give to it.

    temp_dir = pathlib.Path( temp_dir if temp_dir is not None else Config.get().value( 'photometry.snappl.temp_dir' ) )
    barf = "".join( random.choices( "0123456789abcdef", k=10 ) )
    tmpimpath = temp_dir / f"{barf}_sub.fits"
    tmpsubpath = temp_dir / f"{barf}_sub.fits"
    tmpdetmaskpath = temp_dir / f"{barf}_detmask.fits"

    origimg = img
    try:
        if isinstance( origimg, snappl.image.FITSImageOnDisk ):
            img = origimg.uncompressed_version( include=['data'] )
        else:
            img = snappl.image.FITSImage( path=tmpimpath, header=fits.header.Header() )
            img.data = origimg.data
            img.save( which='data' )

        SNLogger.debug( "Calling SEx_SkySubtract.SSS..." )
        ( _SKYDIP, _SKYPEAK, _PixA_skysub,
          _PixA_sky, PixA_skyrms ) = SEx_SkySubtract.SSS(FITS_obj=img.path,
                                                         FITS_skysub=tmpsubpath,
                                                         FITS_detmask=tmpdetmaskpath,
                                                         FITS_sky=None, FITS_skyrms=None,
                                                         ESATUR_KEY='ESATUR',
                                                         BACK_SIZE=64, BACK_FILTERSIZE=3,
                                                         DETECT_THRESH=1.5, DETECT_MINAREA=5,
                                                         DETECT_MAXAREA=0,
                                                         VERBOSE_LEVEL=2, MDIR=None)
        SNLogger.debug( "...back from SEx_SkySubtract.SSS" )

        subim = snappl.image.FITSImage( path=tmpsubpath )
        detmaskim = snappl.image.FITSImage( path=tmpdetmaskpath )
        skyrms = np.median( PixA_skyrms )
        return subim, detmaskim, skyrms

    finally:
        # Clean up the image temp file if necessary
        if img.path != origimg.path:
            img.path.unlink( missing_ok=True )


def stampmaker(ra, dec, shape, img, savedir=None, savename=None):
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
    origimg = img
    try:
        if isinstance( origimg, snappl.image.FITSImageOnDisk ):
            img = origimg.uncompressed_version( include=['data'] )
        else:
            barf = "".join( random.choices( "0123456789abcdef:", k=10 ) )
            img = snappl.image.FITSImage( path=savedir / f"{barf}.fits", header=fits.header.Header() )
            img.data = origimg.data
            img.save_data( which='data' )

        # TODO : if Stamp_Generator.SG can take a Path in FITS_StpLst, remove the str()
        Stamp_Generator.SG(FITS_obj=img.path, EXTINDEX=0, COORD=pxradec, COORD_TYPE='IMAGE',
                           STAMP_IMGSIZE=shape, FILL_VALUE=np.nan, FITS_StpLst=str(savepath))

        return savepath

    finally:
        # Clean up temp file if necessary
        if img.path != origimg.path:
            img.path.unlink( missing_ok=True )
