# IMPORTS Standard:
import os
import numpy as np
import gzip
import shutil
from tempfile import mkdtemp

# IMPORTS Astro:
from astropy.io import fits 
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS

# IMPORTS SFFT:
# from roman_imsim.utils import roman_utils
from sfft.utils.SExSkySubtract import SEx_SkySubtract
from sfft.utils.StampGenerator import Stamp_Generator
from sfft.utils.pyAstroMatic.PYSWarp import PY_SWarp

# IMPORTS Internal:
from .utils import _build_filepath

output_files_rootdir = '/cwork/lna18/imsub_out/'

def gz_decompress(in_path,out_path):
    with gzip.open(in_path,'rb') as f_in, open(out_path,'wb') as f_out:
        shutil.copyfileobj(f_in,f_out)

def sky_subtract(path=None, band=None, pointing=None, sca=None, out_path=output_files_rootdir):

    """
    Subtracts background, found with Source Extractor. 
    """
    original_imgpath = _build_filepath(path=path,band=band,pointing=pointing,sca=sca,filetype='image')
    if not os.path.exists(os.path.join(out_path, 'unzip')):
        os.mkdir(os.path.join(out_path, 'unzip'))

    decompressed_path = os.path.join(output_files_rootdir,'unzip',f'{os.path.basename(original_imgpath)[:-3]}')
    gz_decompress(original_imgpath, decompressed_path)
    output_path = os.path.join(out_path, 'skysub', f'skysub_{os.path.basename(original_imgpath)}.fits')

    if not os.path.exists(os.path.join(out_path, 'skysub')):
        os.mkdir(os.path.join(out_path, 'skysub'))

    SEx_SkySubtract.SSS(FITS_obj=decompressed_path, FITS_skysub=output_path, FITS_sky=None, FITS_skyrms=None, \
                        ESATUR_KEY='ESATUR', BACK_SIZE=64, BACK_FILTERSIZE=3, DETECT_THRESH=1.5, \
                        DETECT_MINAREA=5, DETECT_MAXAREA=0, VERBOSE_LEVEL=2, FITS_DATA_EXT=1)

    return output_path

# def imalign(tab):


def write_stamp(cutout,orig_header,savepath):
    hdu = fits.PrimaryHDU(cutout.data, header=orig_header)
    hdu.header.update(cutout.wcs.to_header())
    hdr = hdu.header
    hdr['OG_XCR'], hdr['OG_YCR'] = np.array(cutout.center_original) - 2044. # Because you padded the original image
    hdr['STMP_XCR'], hdr['STMP_YCR'] = cutout.center_cutout
    hdul = fits.HDUList([hdu])
    hdul.writeto(savepath, overwrite=True)

def stampmaker(ra,dec,
                stamp_savedir,
                size=100,
                path=None,
                band=None,
                pointing=None,
                sca=None,
                ref_path=None,
                ref_band='H158',
                ref_pointing='1394',
                ref_sca='12'):

        """
        For now, you can only input images as reference and image that contain the same coordinates. 
        """

        # After everything is done, this directory gets deleted. 
        tmpdir = mkdtemp(prefix='stmpmk_')
        original_imgpath = _build_filepath(path=path,band=band,pointing=pointing,sca=sca,filetype='image')

        # Extract .fits.gz file from original RomanDESC directory to a temporary place. 
        original_decomp_path = os.path.join(tmpdir, 'original.fits')
        gz_decompress(original_imgpath,original_decomp_path)

        # Also, pad borders with 0 so the output image contains all pixels and nothing is cut off. 
        Stamp_Generator.SG(FITS_obj=original_decomp_path, STAMP_IMGSIZE=[4088*2]*2,
                                COORD=np.array([[2044.,2044.]]), FILL_VALUE=np.nan, EXTINDEX=1, FITS_StpLst=[original_decomp_path])

        # Now, we need the reference to be padded as well so that the output image contains all pixels
        # and is not cut off. 
        refpath = _build_filepath(path=ref_path,band=ref_band,pointing=ref_pointing,sca=ref_sca,filetype='image')
        ref_decomp_path = os.path.join(tmpdir, 'ref.fits')
        gz_decompress(refpath,ref_decomp_path)

        Stamp_Generator.SG(FITS_obj=ref_decomp_path, STAMP_IMGSIZE=[4088*2]*2, 
                            COORD=np.array([[2044.,2044.]]), FILL_VALUE=np.nan, EXTINDEX=1, FITS_StpLst=[ref_decomp_path])

        rotate_path = os.path.join(tmpdir, 'orig_rotated.fits')
        PY_SWarp.PS(FITS_obj=original_decomp_path, FITS_ref=ref_decomp_path, 
                    FITS_resamp=rotate_path, IMAGE_SIZE=4088*2, 
                    NAXIS1_VAL=4088, NAXIS2_VAL=4088)

        rotated_hdu = fits.open(rotate_path)
        rotated_array = rotated_hdu[0].data
        rotated_wcs = WCS(rotated_hdu[0].header, relax=True)
        orig_hdu = fits.open(original_imgpath)
        orig_header = orig_hdu[0].header
        cutout = Cutout2D(rotated_array, position=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), size=size, wcs=rotated_wcs)

        write_stamp(cutout,orig_header,stamp_savedir)

# def get_psf(band,pointing,sca):
#     """
#     1. Get PSF model at image center. 
#     """

#     imgpath = _build_filepath(path=None,band=band,pointing=pointing,sca=sca,filetype='image')
#     imginfo = roman_utils(config_file,visit=pointing,sca=sca)