__all__ = [ 'gz_and_ext', 'sky_subtract', 'stampmaker' ]

# IMPORTS Standard:
import os
import io
import numpy as np
import gzip
import multiprocessing
import shutil
import pathlib

# IMPORTS Astro:
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u
from astropy.wcs import WCS

# IMPORTS SFFT:
from sfft.utils.SExSkySubtract import SEx_SkySubtract
from sfft.utils.StampGenerator import Stamp_Generator

# IMPORTS internal
from snpit_utils.logger import SNLogger
from snpit_utils.config import Config


def gz_and_ext(in_path, out_path):
    """Utility function that unzips the original file and turns
       it into a single-extension FITS file.

    Parameters
    ----------
      in_path: Path
        The input path of the file to unzip and flatten.

      out_path: Path
        The desired output path of the file that is unzipped and flattened, once it is saved.

    Returns
    -------
      out_path: Path
        The output path of the file that has been unzipped and flattened.

    """

    bio = io.BytesIO()
    with gzip.open(in_path, 'rb') as f_in:
        shutil.copyfileobj(f_in, bio)
    bio.seek(0)

    with fits.open(bio) as hdu:
        newhdu = fits.HDUList([fits.PrimaryHDU(data=hdu[1].data, header=hdu[0].header)])
        SNLogger.debug( f"Process {multiprocessing.current_process().pid} Writing extracted HDU 1 to {out_path}..." )
        newhdu.writeto(out_path, overwrite=True)
        SNLogger.debug( f"...process {multiprocessing.current_process().pid} done writing "
                        f"extracted HDU 1 to {out_path}." )

    return out_path


def sky_subtract( inpath, skysubpath, detmaskpath, temp_dir=pathlib.Path("/tmp"), force=False ):
    """Subtracts background, found with Source Extractor.

    Parameters
    ----------
      inpath: Path
        Original FITS image

      skysubpath: Path
        Sky-subtracted FITS image

      detmaskpath: Path
        Detection Mask FITS Image.  (Will be uint8, I think.)

      temp_dir: Path
        Already-existing directory where we can write a temporary file.
        (If the image is .gz compressed, source-extractor can't handle
        that, so we have to write a decompressed version.)

      force: bool, default False
        If False, and outpath already exists, do nothing.  If True,
        clobber the existing file and recalculate it.

    Returns
    -------
      skyrms: float
        Median of the skyrms image calculated by source-extractor

    """

    if ( not force ) and ( skysubpath.is_file() ) and ( detmaskpath.is_file() ):
        with fits.open( skysubpath ) as hdul:
            skyrms = hdul[0].header['SKYRMS']
        return skyrms

    # do_skysub = force or ( ( not force ) and ( not skysubpath.is_file() ) and ( detmaskpath.is_file() ) )

    if inpath.name[-3:] == '.gz':
        decompressed_path = temp_dir / inpath.name[:-3]
        gz_and_ext( inpath, decompressed_path )
    else:
        decompressed_path = inpath


    SNLogger.debug( "Calling SEx_SkySubtract.SSS..." )
    ( _SKYDIP, _SKYPEAK, _PixA_skysub,
      _PixA_sky, PixA_skyrms ) = SEx_SkySubtract.SSS(FITS_obj=decompressed_path,
                                                     FITS_skysub=skysubpath,
                                                     FITS_detmask=detmaskpath,
                                                     FITS_sky=None, FITS_skyrms=None,
                                                     ESATUR_KEY='ESATUR',
                                                     BACK_SIZE=64, BACK_FILTERSIZE=3,
                                                     DETECT_THRESH=1.5, DETECT_MINAREA=5,
                                                     DETECT_MAXAREA=0,
                                                     VERBOSE_LEVEL=2, MDIR=None)
    SNLogger.debug( "...back from SEx_SkySubtract.SSS" )

    return np.median( PixA_skyrms )


def stampmaker(ra, dec, shape, imgpath, savedir=None, savename=None):
    """Make stamps.

    Parameters
    ----------
      ra: float
        RA of center of stamp in degrees.

      dec: float
        Dec of center of stamp in degrees.

      shape: np.array
        Shape of stamp. must be a numpy array. e.g. np.array([100,100])

      imgpath: Path
        Path to image that stamp will be extracted from.

      savedir: Path
        Directory stamp will be saved to.

      savename: Path
        Base name of stamp filepath.

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

    if savename is None:
        savename = f'stamp_{os.path.basename(imgpath)}.fits'

    savepath = savedir / savename

    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    with fits.open(imgpath) as hdu:
        hdun = 1 if str(imgpath)[-3:] in ('.fz', '.gz') else 0
        wcs = WCS(hdu[hdun].header)

    x, y = skycoord_to_pixel(coord, wcs)
    pxradec = np.array([[x, y]])

    # TODO : if Stamp_Generator.SG can take a Path in FITS_StpLst, remove the str()
    Stamp_Generator.SG(FITS_obj=imgpath, EXTINDEX=hdun, COORD=pxradec, COORD_TYPE='IMAGE',
                       STAMP_IMGSIZE=shape, FILL_VALUE=np.nan, FITS_StpLst=str(savepath))

    return savepath
