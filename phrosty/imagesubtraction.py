# IMPORTS Standard:
import os, sys
import numpy as np
import gzip
import shutil
from tempfile import mkdtemp
from glob import glob

# IMPORTS Astro:
from astropy.io import fits 
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u
from astropy.wcs import WCS
from astropy.convolution import convolve_fft

# IMPORTS SFFT:
from roman_imsim.utils import roman_utils
from sfft.utils.SExSkySubtract import SEx_SkySubtract
from sfft.utils.StampGenerator import Stamp_Generator
from sfft.utils.pyAstroMatic.PYSWarp import PY_SWarp
from sfft.utils.ReadWCS import Read_WCS
from sfft.utils.ImageZoomRotate import Image_ZoomRotate

# IMPORTS Internal:
from .utils import _build_filepath

output_files_rootdir = '/work/lna18/imsub_out/'

def gz_and_ext(in_path,out_path):
    """
    Unzips the original file and turns it into a single-extension FITS file.     
    """
    with gzip.open(in_path,'rb') as f_in, open(out_path,'wb') as f_out:
        shutil.copyfileobj(f_in,f_out)
    
    with fits.open(out_path) as hdu:
        newhdu = fits.HDUList([fits.PrimaryHDU(data=hdu[1].data, header=hdu[0].header)])
        newhdu.writeto(out_path, overwrite=True)

def sky_subtract(path=None, band=None, pointing=None, sca=None, out_path=output_files_rootdir, remove_tmpdir=True):

    """
    Subtracts background, found with Source Extractor. 
    """
    original_imgpath = _build_filepath(path=path,band=band,pointing=pointing,sca=sca,filetype='image')
    if not os.path.exists(os.path.join(out_path, 'unzip')):
        os.mkdir(os.path.join(out_path, 'unzip'))

    decompressed_path = os.path.join(output_files_rootdir,'unzip',f'{os.path.basename(original_imgpath)[:-3]}')
    gz_and_ext(original_imgpath, decompressed_path)
    output_path = os.path.join(out_path, 'skysub', f'skysub_{os.path.basename(decompressed_path)}')

    if not os.path.exists(os.path.join(out_path, 'skysub')):
        os.mkdir(os.path.join(out_path, 'skysub'))

    SEx_SkySubtract.SSS(FITS_obj=decompressed_path, FITS_skysub=output_path, FITS_sky=None, FITS_skyrms=None, \
                        ESATUR_KEY='ESATUR', BACK_SIZE=64, BACK_FILTERSIZE=3, DETECT_THRESH=1.5, \
                        DETECT_MINAREA=5, DETECT_MAXAREA=0, VERBOSE_LEVEL=2, MDIR=output_files_rootdir)

    if remove_tmpdir:
        tmpdir = glob(os.path.join(output_files_rootdir,'PYSEx_*'))
        for tdir in tmpdir:
            shutil.rmtree(tdir, ignore_errors=True)

    return output_path

def imalign(template_path, sci_path, out_path=output_files_rootdir, remove_tmpdir=True):
    """
    Align images with SWarp. 
    """
    outdir = os.path.join(out_path, 'align')
    output_path = os.path.join(outdir, f'align_{os.path.basename(sci_path)}')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    PY_SWarp.PS(FITS_obj=sci_path, FITS_ref=template_path, FITS_resamp=output_path, \
                GAIN_KEY='GAIN', SATUR_KEY='SATURATE', OVERSAMPLING=1, RESAMPLING_TYPE='BILINEAR', \
                SUBTRACT_BACK='N', FILL_VALUE=np.nan, VERBOSE_TYPE='NORMAL', VERBOSE_LEVEL=1, TMPDIR_ROOT=output_files_rootdir)

    if remove_tmpdir:
        tmpdir = glob(os.path.join(output_files_rootdir,'PYSWarp_*'))
        for tdir in tmpdir:
            shutil.rmtree(tdir, ignore_errors=True)

    return output_path

def calculate_skyN_vector(wcshdr, x_start, y_start, shift_dec=1.0):
    """
    Lei Hu 2024
    """
    w = Read_WCS.RW(wcshdr, VERBOSE_LEVEL=1)
    ra_start, dec_start = w.all_pix2world(np.array([[x_start, y_start]]), 1)[0]
    ra_end, dec_end = ra_start, dec_start + shift_dec/3600.0
    x_end, y_end = w.all_world2pix(np.array([[ra_end, dec_end]]), 1)[0]
    skyN_vector = np.array([x_end - x_start, y_end - y_start])
    return skyN_vector

def calculate_rotate_angle(vector_ref, vector_obj):
    """
    Lei Hu 2024
    """
    rad = np.arctan2(np.cross(vector_ref, vector_obj), np.dot(vector_ref, vector_obj))
    rotate_angle = np.rad2deg(rad)
    if rotate_angle < 0.0: 
        rotate_angle += 360.0 
    return rotate_angle

def get_psf(ra,dec,sci_imalign,
                    sci_skysub,
                    sci_band,
                    sci_pointing,
                    sci_sca,
                    ref_band,
                    ref_pointing,
                    ref_sca):
    """
    1. Get PSF model at specified RA, dec [degrees]. 
    2. Rotate PSF model to match reference WCS. 
        2a. Calculate rotation angle during alignment
        2b. Rotate PSF to match rotated science image
    """

    ref_path = _build_filepath(path=None,band=ref_band,pointing=ref_pointing,sca=ref_sca,filetype='image')
    ref_hdu = fits.open(ref_path)
    ref_wcs = WCS(ref_hdu[0].header)
    worldcoords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    pxradec = skycoord_to_pixel(worldcoords,ref_wcs)

    # Get PSF at specified RA, dec in science image. 
    config_path = os.path.join(os.path.dirname(__file__), 'auxiliary', 'tds.yaml')
    config = roman_utils(config_path,visit=sci_pointing,sca=sci_sca)
    psf = config.getPSF_Image(501, x=pxradec[0], y=pxradec[1]).array

    # Get vector from sky-subtracted science WCS
    hdr = fits.getheader(sci_skysub, ext=0)
    _w = Read_WCS.RW(hdr, VERBOSE_LEVEL=1)
    x0, y0 = 0.5 + int(hdr['NAXIS1'])/2.0, 0.5 + int(hdr['NAXIS2'])/2.0
    ra0, dec0 = _w.all_pix2world(np.array([[x0, y0]]), 1)[0]
    skyN_vector = calculate_skyN_vector(wcshdr=hdr, x_start=x0, y_start=y0)

    # Get vector from rotated science WCS
    hdr = fits.getheader(sci_imalign, ext=0)
    _w = Read_WCS.RW(hdr, VERBOSE_LEVEL=1)
    x1, y1 = _w.all_world2pix(np.array([[ra0, dec0]]), 1)[0]
    skyN_vectorp = calculate_skyN_vector(wcshdr=hdr, x_start=x1, y_start=y1)
    PATTERN_ROTATE_ANGLE = calculate_rotate_angle(vector_ref=skyN_vector, vector_obj=skyN_vectorp)

    # Do rotation
    psf_rotated = Image_ZoomRotate.IZR(PixA_obj=psf, ZOOM_SCALE_X=1., \
                                        ZOOM_SCALE_Y=1., PATTERN_ROTATE_ANGLE=PATTERN_ROTATE_ANGLE, \
                                        RESAMPLING_TYPE='BILINEAR', FILL_VALUE=0.0, VERBOSE_LEVEL=1)[0]

    # Save rotated PSF
    psf_dir = os.path.join(output_files_rootdir,'psf')
    if not os.path.exists(psf_dir):
        os.mkdir(psf_dir)

    psf_path = os.path.join(psf_dir, f'rot_psf_{ra}_{dec}_{sci_band}_{sci_pointing}_{sci_sca}.fits')
    fits.HDUList([fits.PrimaryHDU(data=psf_rotated, header=None)]).writeto(psf_path, overwrite=True)
    
    return psf_path

def crossconvolve(sci_img_path, sci_psf_path,
                    ref_img_path, ref_psf_path):

    savedir = os.path.join(output_files_rootdir,'convolved')
    if not os.path.exists(savedir):
        os.mkdir(savedir)

    # First convolves reference PSF on science image. 
    # Then, convolves science PSF on reference image. 
    savepaths = []
    for img, psf, name in zip([sci_img_path, ref_img_path],
                        [ref_psf_path, sci_psf_path],
                        ['sci', 'ref']):
        imgdata = fits.getdata(img, ext=0).T
        psfdata = fits.getdata(psf, ext=0).T

        convolved = convolve_fft(imgdata, psfdata, boundary='fill', nan_treatment='fill', \
                                fill_value=0.0, normalize_kernel=True)

        savename = f'conv_{os.path.basename(img)}'
        savepath = os.path.join(savedir, savename)

        with fits.open(img) as hdl:
            hdl[0].data[:, :] = convolved.T
            hdl.writeto(savepath, overwrite=True)
            savepaths.append(savepath)
    return savepaths

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