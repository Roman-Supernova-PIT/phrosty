# temp import for use with nsys
# (see "with nvtx.annotate" blocks below)
import nvtx

# IMPORTS Standard:
import os
import sys
import numpy as np
import cupy
import cupyx.scipy
import gzip
import shutil
import logging
import matplotlib.pyplot as plt
import tracemalloc
import json

from numba import cuda

# IMPORTS Astro:
from astropy.io import fits 
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u
from astropy.wcs import WCS
from astropy.convolution import convolve_fft
from astropy.visualization import ZScaleInterval
from galsim import PositionD
# from photutils.psf import FittableImageModel

# IMPORTS SFFT:
from roman_imsim.utils import roman_utils
from sfft.utils.SExSkySubtract import SEx_SkySubtract
from sfft.utils.StampGenerator import Stamp_Generator
from sfft.utils.pyAstroMatic.PYSWarp import PY_SWarp
from sfft.utils.ReadWCS import Read_WCS
from sfft.utils.ImageZoomRotate import Image_ZoomRotate
from sfft.utils.CudaResampling import Cuda_Resampling
from sfft.utils.pyAstroMatic.PYSEx import PY_SEx
from sfft.CustomizedPacket import Customized_Packet
from sfft.utils.SkyLevelEstimator import SkyLevel_Estimator
from sfft.utils.SFFTSolutionReader import Realize_MatchingKernel
from sfft.utils.DeCorrelationCalculator import DeCorrelation_Calculator
# from sfft.utils.meta.MultiProc import Multi_Proc

# IMPORTS Internal:
from .utils import _build_filepath, get_transient_radec, get_transient_mjd, get_fitsobj

# Configure logger (Rob)
_logger = logging.getLogger(f'phrosty')
if not _logger.hasHandlers():
    log_out = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(f'[%(asctime)s - %(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    log_out.setFormatter(formatter)
    _logger.addHandler(log_out)
    _logger.setLevel(logging.DEBUG) # ERROR, WARNING, INFO, or DEBUG (in that order by increasing detail)


"""
This module was written with significant contributions from 
Dr. Lei Hu (https://github.com/thomasvrussell/), and relies on his 
SFFT image subtraction package (https://github.com/thomasvrussell/sfft).
"""

output_files_rootdir = os.getenv('DIA_OUT_DIR', None)
assert output_files_rootdir is not None, 'You need to set DIA_OUT_DIR as an environment variable.'

def check_and_mkdir(dirname):
    """
    Utility function for checking if a directory exists, and if not,
    makes that directory. 
    """
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def gz_and_ext(in_path,out_path):
    """
    Utility function that unzips the original file and turns it into a single-extension FITS file.     
    """
    with gzip.open(in_path,'rb') as f_in, open(out_path,'wb') as f_out:
        shutil.copyfileobj(f_in,f_out)
    
    with fits.open(out_path) as hdu:
        newhdu = fits.HDUList([fits.PrimaryHDU(data=hdu[1].data, header=hdu[0].header)])
        newhdu.writeto(out_path, overwrite=True)

    return out_path

def sky_subtract(path=None, band=None, pointing=None, sca=None, out_path=output_files_rootdir, force=False, verbose=False):

    """
    Subtracts background, found with Source Extractor. 
    """
    original_imgpath = _build_filepath(path=path,band=band,pointing=pointing,sca=sca,filetype='image')
    zip_savedir = os.path.join(out_path, 'unzip')
    check_and_mkdir(zip_savedir)

    if f'{os.path.basename(original_imgpath)[-3:]}' == '.gz':
        decompressed_path = os.path.join(output_files_rootdir,'unzip',f'{os.path.basename(original_imgpath)[:-3]}')
        gz_and_ext(original_imgpath, decompressed_path)
    else:
        decompressed_path = os.path.join(zip_savedir,f'{os.path.basename(original_imgpath)}')
        
    output_path = os.path.join(out_path, 'skysub', f'skysub_{os.path.basename(decompressed_path)}')

    do_skysub = (force is True) or (force is False and not os.path.exists(output_path))
    skip_skysub = (not force) and os.path.exists(output_path)
    if do_skysub:
        sub_savedir = os.path.join(out_path, 'skysub')
        check_and_mkdir(sub_savedir)

        ( SKYDIP, SKYPEAK, PIxA_skysub,
          PixA_sky, PixA_skyrms, DETECT_MASK ) = SEx_SkySubtract.SSS(FITS_obj=decompressed_path, FITS_skysub=output_path, 
                                                                     FITS_sky=None, FITS_skyrms=None,
                                                                     ESATUR_KEY='ESATUR', BACK_SIZE=64, BACK_FILTERSIZE=3, 
                                                                     DETECT_THRESH=1.5, DETECT_MINAREA=5, DETECT_MAXAREA=0, 
                                                                     VERBOSE_LEVEL=2, MDIR=None)
    elif skip_skysub and verbose:
        print(output_path, 'already exists. Skipping sky subtraction.')

    return output_path, PixA_sky, PixA_skyrms, DETECT_MASK


def run_resample(FITS_obj, FITS_targ, FITS_resamp):
    """run resampling using CUDA"""
    hdr_obj = fits.getheader(FITS_obj, ext=0)
    hdr_targ = fits.getheader(FITS_targ, ext=0)

    PixA_obj = fits.getdata(FITS_obj, ext=0).T
    PixA_targ = fits.getdata(FITS_targ, ext=0).T

    if not PixA_obj.flags['C_CONTIGUOUS']:
        PixA_obj = np.ascontiguousarray(PixA_obj, np.float64)
        PixA_obj_GPU = cupy.array(PixA_obj)
    else: PixA_obj_GPU = cupy.array(PixA_obj.astype(np.float64))

    if not PixA_targ.flags['C_CONTIGUOUS']:
        PixA_targ = np.ascontiguousarray(PixA_targ, np.float64)
        PixA_targ_GPU = cupy.array(PixA_targ)
    else: PixA_targ_GPU = cupy.array(PixA_targ.astype(np.float64))

    CR = Cuda_Resampling(RESAMP_METHOD='BILINEAR', VERBOSE_LEVEL=1)
    XX_proj_GPU, YY_proj_GPU = CR.projection_sip(hdr_obj, hdr_targ, Nsamp=1024, RANDOM_SEED=10086)
    PixA_Eobj_GPU, EProjDict = CR.frame_extension(XX_proj_GPU=XX_proj_GPU, YY_proj_GPU=YY_proj_GPU, PixA_obj_GPU=PixA_obj_GPU)
    PixA_resamp = cupy.asnumpy(CR.resampling(PixA_Eobj_GPU=PixA_Eobj_GPU, EProjDict=EProjDict))

    with fits.open(FITS_targ) as hdl:
        PixA_resamp[PixA_resamp == 0.] = np.nan
        hdl[0].data[:, :] = PixA_resamp.T
        hdl.writeto(FITS_resamp, overwrite=True)
    return None


def imalign(template_path, sci_path, out_path=output_files_rootdir,savename=None,force=False, verbose=False):
    """
    Align images. 
    """

    outdir = os.path.join(out_path, 'align')
    if savename is None:
        savename = os.path.basename(sci_path)
    output_path = os.path.join(outdir, f'align_{savename}')


    do_align = (force is True) or (force is False and not os.path.exists(output_path))
    skip_align = (not force) and os.path.exists(output_path)
    if do_align:
        check_and_mkdir(outdir)

        _logger.debug( "Using Cuda_Resampling.CR to resample image" )
        
        # Cuda_Resampling.CR( sci_path, template_path, output_path, METHOD="BILINEAR" )
        run_resample( sci_path, template_path, output_path )
        
    elif skip_align and verbose:
        print(output_path, 'already exists. Skipping alignment.')

    return output_path

# def cuda_imalign(template_path, sci_path,
#                  template_psf_path, sci_psf_path,
#                  template_detection_mask_path, sci_detection_mask_path,
#                  out_path=output_files_rootdir,savename=None,
#                  force=False, verbose=False
#                 ):

#     outdir = os.path.join(out_path, 'align')
#     if savename is None:
#         savename = os.path.basename(sci_path)
#     output_path = os.path.join(outdir, f'align_{savename}')

#     do_align = (force is True) or (force is False and not os.path.exists(output_path))
#     skip_align = (not force) and os.path.exists(output_path)
#     if do_align:
#         check_and_mkdir(outdir)

#         _logger.debug( "Using SpaceSFFTCupyFlow to resample image" )

#         with nvtx.annotate("toGPU", color="#714CF9"):

#             FITS_REF = template_path
#             FITS_oSCI = sci_path

#             FITS_PSF_REF = template_psf_path
#             FITS_PSF_oSCI = sci_psf_path

#             FITS_REF_DMASK = template_detection_mask_path
#             FITS_oSCI_DMASK = sci_detection_mask_path

#             hdr_REF = fits.getheader(FITS_REF, ext=0)
#             hdr_oSCI = fits.getheader(FITS_oSCI, ext=0)

#             PixA_REF = fits.getdata(FITS_REF, ext=0).T
#             PixA_oSCI = fits.getdata(FITS_oSCI, ext=0).T

#             PixA_REF_DMASK = fits.getdata(FITS_REF_DMASK, ext=0).T
#             PixA_oSCI_DMASK = fits.getdata(FITS_oSCI_DMASK, ext=0).T

#             PSF_REF = fits.getdata(FITS_PSF_REF, ext=0).T
#             PSF_oSCI = fits.getdata(FITS_PSF_oSCI, ext=0).T

#             if not PixA_REF.flags['C_CONTIGUOUS']:
#                 PixA_REF = np.ascontiguousarray(PixA_REF, np.float64)
#                 PixA_REF_GPU = cp.array(PixA_REF)
#             else: PixA_REF_GPU = cp.array(PixA_REF.astype(np.float64))

#             if not PixA_oSCI.flags['C_CONTIGUOUS']:
#                 PixA_oSCI = np.ascontiguousarray(PixA_oSCI, np.float64)
#                 PixA_oSCI_GPU = cp.array(PixA_oSCI)
#             else: PixA_oSCI_GPU = cp.array(PixA_oSCI.astype(np.float64))

#             if not PixA_REF_DMASK.flags['C_CONTIGUOUS']:
#                 PixA_REF_DMASK = np.ascontiguousarray(PixA_REF_DMASK, np.float64)
#                 PixA_REF_DMASK_GPU = cp.array(PixA_REF_DMASK)
#             else: PixA_REF_DMASK_GPU = cp.array(PixA_REF_DMASK.astype(np.float64))

#             if not PixA_oSCI_DMASK.flags['C_CONTIGUOUS']:
#                 PixA_oSCI_DMASK = np.ascontiguousarray(PixA_oSCI_DMASK, np.float64)
#                 PixA_oSCI_DMASK_GPU = cp.array(PixA_oSCI_DMASK)
#             else: PixA_oSCI_DMASK_GPU = cp.array(PixA_oSCI_DMASK.astype(np.float64))

#             if not PSF_REF.flags['C_CONTIGUOUS']:
#                 PSF_REF = np.ascontiguousarray(PSF_REF, np.float64)
#                 PSF_REF_GPU = cp.array(PSF_REF)
#             else: PSF_REF_GPU = cp.array(PSF_REF.astype(np.float64))

#             if not PSF_oSCI.flags['C_CONTIGUOUS']:
#                 PSF_oSCI = np.ascontiguousarray(PSF_oSCI, np.float64)
#                 PSF_oSCI_GPU = cp.array(PSF_oSCI)
#             else: PSF_oSCI_GPU = cp.array(PSF_oSCI.astype(np.float64))

#         PixA_DIFF_GPU, PixA_DCDIFF_GPU, PixA_DSNR_GPU = SpaceSFFT_CupyFlow_NVTX.SSCFN(hdr_REF=hdr_REF, hdr_oSCI=hdr_oSCI, 
#             PixA_REF_GPU=PixA_REF_GPU, PixA_oSCI_GPU=PixA_oSCI_GPU, PixA_REF_DMASK_GPU=PixA_REF_DMASK_GPU, 
#             PixA_oSCI_DMASK_GPU=PixA_oSCI_DMASK_GPU, PSF_REF_GPU=PSF_REF_GPU, PSF_oSCI_GPU=PSF_oSCI_GPU, 
#             GKerHW=9, KerPolyOrder=2, BGPolyOrder=0, ConstPhotRatio=True, CUDA_DEVICE_4SUBTRACT='0', GAIN=1.0)

#         with nvtx.annotate("toCPU", color="#FD89B4"):
#             PixA_DCDIFF = cp.asnumpy(PixA_DCDIFF_GPU)
#             PixA_DSNR = cp.asnumpy(PixA_DSNR_GPU)

        # save the final results
        # FITS_DCDIFF = output_dir + "/%s.sciE.skysub.resamp.sfftdiff.decorr.fits" % sciname
        # with fits.open(FITS_REF) as hdul:
        #     hdul[0].data = PixA_DCDIFF.T
        #     hdul.writeto(FITS_DCDIFF, overwrite=True)

        # FITS_DSNR = output_dir + "/%s.sciE.skysub.resamp.sfftdiff.decorr.snr.fits" % sciname
        # with fits.open(FITS_REF) as hdul:
        #     hdul[0].data = PixA_DSNR.T
        #     hdul.writeto(FITS_DSNR, overwrite=True)
        # print("Done!")

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

def get_imsim_psf(ra, dec, band, pointing, sca, size=201, config_yaml_file=None,
                  out_path=output_files_rootdir, force=False, logger=None):

    """
    Retrieve the PSF from roman_imsim/galsim, and transform the WCS so that CRPIX and CRVAL
    are centered on the image instead of at the corner. 

    force parameter does not currently do anything.
    """

    logger = _logger if logger is None else logger

    savedir = os.path.join(out_path,'psf')
    check_and_mkdir(savedir)
    savename = f'psf_{ra}_{dec}_{band}_{pointing}_{sca}.fits'
    savepath = os.path.join(savedir,savename)

    # Get WCS of the image you need the PSF for.
    hdu = get_fitsobj(band=band,pointing=pointing,sca=sca)
    wcs = WCS(hdu[0].header)
    coord = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
    x,y = wcs.world_to_pixel(coord)

    # Get PSF at specified ra, dec.
    assert config_yaml_file is not None, "config_yaml_file is a required argument"
    config_path = config_yaml_file
    config = roman_utils(config_path,pointing,sca)
    psf = config.getPSF_Image(size,x,y)
    psf.write(savepath)

    # Change the WCS so CRPIX and CRVAL are centered. 
    pos = PositionD(x=x,y=y)
    wcs_new = psf.wcs.affine(image_pos=pos)
    psf.wcs = wcs_new

    # Save fits object.
    psf.write(savepath)

    # Transpose the data array so it works with SFFT. 
    # Can't do this like psf.array = psf.array.T because you get an error:
    # "property 'array' of 'Image' object has no setter".
    hdu = fits.open(savepath)
    hdu[0].data = hdu[0].data.T
    hdu[0].header['CRVAL1'] = ra
    hdu[0].header['CRVAL2'] = dec
    hdu[0].header['CRPIX1'] = 0.5 + int(hdu[0].header['NAXIS1'])/2.
    hdu[0].header['CRPIX2'] = 0.5 + int(hdu[0].header['NAXIS2'])/2.
    hdu.writeto(savepath,overwrite=True)

    return savepath


def rotate_psf(ra,dec,psf,target,savename=None,force=False,verbose=False):
    """
    2. Rotate PSF model to match reference WCS. 
        2a. Calculate rotation angle during alignment
        2b. Rotate PSF to match rotated science image
    """

    # Set up filepaths.
    psf_dir = os.path.join(output_files_rootdir,'psf')
    check_and_mkdir(psf_dir)

    basename = os.path.basename(psf)
    if savename is None:
        savename = f'rot_{basename}'
    psf_path = os.path.join(psf_dir, savename)

    do_psf = (force is True) or (force is False and not os.path.exists(psf_path))
    skip_psf = (not force) and os.path.exists(psf_path)
    if do_psf:
        # Get vector from original PSF WCS (i.e., not rotated to reference)
        hdr = fits.getheader(psf, ext=0)
        _w = Read_WCS.RW(hdr, VERBOSE_LEVEL=1)
        x0, y0 = 0.5 + int(hdr['NAXIS1'])/2.0, 0.5 + int(hdr['NAXIS2'])/2.0
        ra0, dec0 = _w.all_pix2world(np.array([[x0, y0]]), 1)[0]
        skyN_vector = calculate_skyN_vector(wcshdr=hdr, x_start=x0, y_start=y0)

        # Also get the PSF image for rotation
        psfimg = fits.getdata(psf, ext=0).T # Already saved as a transposed matrix from get_imsim_psf. 

        # Get vector from target WCS (i.e., rotated)
        hdr = fits.getheader(target, ext=0)
        _w = Read_WCS.RW(hdr, VERBOSE_LEVEL=1)
        x1, y1 = _w.all_world2pix(np.array([[ra0, dec0]]), 1)[0]
        skyN_vectorp = calculate_skyN_vector(wcshdr=hdr, x_start=x1, y_start=y1)
        PATTERN_ROTATE_ANGLE = calculate_rotate_angle(vector_ref=skyN_vector, vector_obj=skyN_vectorp)

        # Do rotation
        psf_rotated = Image_ZoomRotate.IZR(PixA_obj=psfimg, ZOOM_SCALE_X=1., \
                                            ZOOM_SCALE_Y=1., PATTERN_ROTATE_ANGLE=PATTERN_ROTATE_ANGLE, \
                                            RESAMPLING_TYPE='BILINEAR', FILL_VALUE=0.0, VERBOSE_LEVEL=1)[0]

        # Save rotated PSF
        fits.HDUList([fits.PrimaryHDU(data=psf_rotated.T, header=None)]).writeto(psf_path, overwrite=True)
    elif skip_psf and verbose:
        print(psf_path, 'already exists. Skipping getting PSF.')

    return psf_path

def crossconvolve(sci_img_path, sci_psf_path,
                  ref_img_path, ref_psf_path,
                  force=False,verbose=False,
                  sci_outname=None,
                  ref_outname=None,
                  out_path=output_files_rootdir):

    savedir = os.path.join(out_path,'convolved')
    check_and_mkdir(savedir)

    # First convolves reference PSF on science image. 
    # Then, convolves science PSF on reference image. 
    savepaths = []
    for img, name in zip([sci_img_path,ref_img_path],
                    [sci_outname,ref_outname]):

        if name is None:
            savename = f'conv_{os.path.basename(img)}'
        else:
            savename = name

        savepath = os.path.join(savedir, savename)
        savepaths.append(savepath)

    do_conv = (force is True) or (force is False and not any([os.path.exists(p) for p in savepaths]))
    skip_conv = (not force) and all([os.path.exists(p) for p in savepaths])
    if do_conv:
        for img, psf, name, save in zip([sci_img_path, ref_img_path],
                            [ref_psf_path, sci_psf_path],
                            ['sci', 'ref'], savepaths):

            imgdata = fits.getdata(img, ext=0).T
            psfdata = fits.getdata(psf, ext=0).T

            convolved = convolve_fft(imgdata, psfdata, boundary='fill', nan_treatment='fill', \
                                     fill_value=0.0, normalize_kernel=True, preserve_nan=True, allow_huge=True)

            with fits.open(img) as hdl:
                hdl[0].data[:, :] = convolved.T
                hdl.writeto(save, overwrite=True)

    elif skip_conv and verbose:
        print(savepaths, 'already exist. Skipping cross convolve.')

    return savepaths

def stampmaker(ra, dec, imgpath, savedir=None, savename=None, shape=np.array([1000,1000])):

    if savedir is None:
        savedir = os.path.join(output_files_rootdir,'stamps')
        check_and_mkdir(savedir)

    if savename is None:
        savename = f'stamp_{os.path.basename(imgpath)}.fits'

    savepath = os.path.join(savedir,savename)

    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    with fits.open(imgpath) as hdu:
        wcs = WCS(hdu[0].header)
    
    x, y = skycoord_to_pixel(coord, wcs)
    pxradec = np.array([[x,y]])
    
    Stamp_Generator.SG(FITS_obj=imgpath, COORD=pxradec, COORD_TYPE='IMAGE',
                       STAMP_IMGSIZE=shape, FILL_VALUE=np.nan, FITS_StpLst=savepath)

    return savepath

def bkg_mask(imgpath):
    """
    Create detection mask. Not necessary to call directly. 
    """

    source_ext_params = ['X_IMAGE', 'Y_IMAGE', 'FLUX_AUTO', 'FLUXERR_AUTO', 'MAG_AUTO', 'MAGERR_AUTO', 'FLAGS', \
    'FLUX_RADIUS', 'FWHM_IMAGE', 'A_IMAGE', 'B_IMAGE', 'KRON_RADIUS', 'THETA_IMAGE', 'SNR_WIN']

    scatalog = PY_SEx.PS(FITS_obj=imgpath, SExParam=source_ext_params, GAIN_KEY='GAIN', SATUR_KEY='SATURATE', \
                         BACK_TYPE='MANUAL', BACK_VALUE=0.0, BACK_SIZE=64, BACK_FILTERSIZE=3, DETECT_THRESH=1.5, \
                         DETECT_MINAREA=5, DETECT_MAXAREA=0, DEBLEND_MINCONT=0.001, BACKPHOTO_TYPE='LOCAL', \
                         CHECKIMAGE_TYPE='SEGMENTATION', AddRD=True, ONLY_FLAGS=None, XBoundary=0.0, YBoundary=0.0, \
                         DEFAULT_GAIN=1.0, DEFAULT_SATUR=100000, MDIR=None, VERBOSE_LEVEL=1)[1][0]

    bkg_mask = (scatalog == 0)

    # BAD IDEA : we were reading masks from different images
    #   (unaligned template, convolved neither.)
    # fname = os.path.join( output_files_rootdir, f'detect_mask/{os.path.basename(imgpath)}.npy' )
    # _logger.info( f"Trying to load detection mask from {fname}" )
    # bkg_mask = np.load( fname )
    
    return bkg_mask

def difference(scipath, refpath, 
               out_path=output_files_rootdir, savename=None, ForceConv='REF', GKerHW=9, KerPolyOrder=2, BGPolyOrder=0, 
               ConstPhotRatio=True, backend='Numpy', cudadevice='0', nCPUthreads=1, force=False, verbose=False, logger=None):

    tracemalloc.start()

    sci_basename = os.path.basename(scipath)

    sci_data = fits.getdata(scipath).T
    ref_data = fits.getdata(refpath).T

    if savename is None:
        savename = sci_basename

    savedir = os.path.join(out_path, 'subtract')
    check_and_mkdir(savedir)

    diff_savedir = os.path.join(savedir,'difference')
    soln_savedir = os.path.join(savedir, 'solution')
    masked_savedir = os.path.join(savedir, 'masked')

    for dirname in [diff_savedir, soln_savedir, masked_savedir]:
        check_and_mkdir(dirname)

    diff_savepath = os.path.join(diff_savedir, f'diff_{savename}')
    soln_savepath = os.path.join(soln_savedir, f'solution_{savename}')

    sci_masked_savepath = os.path.join(masked_savedir,f'masked_{savename}')
    ref_masked_savepath = os.path.join(masked_savedir,f'masked_{savename}')

    do_subtract = (force is True) or (force is False and not os.path.exists(diff_savepath))
    skip_subtract = (not force) and os.path.exists(diff_savepath)
    if do_subtract:
        # Make combined detection mask.
        sci_bkgmask = bkg_mask(scipath)
        ref_bkgmask = bkg_mask(refpath)
        nanmask = np.isnan(sci_data) | np.isnan(ref_data)
        _bkgmask = np.logical_and(sci_bkgmask,ref_bkgmask)
        bkgmask = np.logical_or(nanmask, _bkgmask)

        for path, msavepath in zip([refpath, scipath], \
                                    [ref_masked_savepath, sci_masked_savepath]):
            with fits.open(path) as hdu:
                hdudata = hdu[0].data.T
                hdudata[bkgmask] = 0.0
                hdu[0].data[:, :] = hdudata.T
                hdu.writeto(msavepath, overwrite=True)

        size,peak = tracemalloc.get_traced_memory()
        logger.debug(f'MEMORY IN imagesubtraction.difference() BEFORE Customized_Packet.CP: size = {size}, peak = {peak}')
        tracemalloc.reset_peak()

        # Do SFFT subtraction
        with nvtx.annotate( "Customized_Packet.CP", color=0x44ff44 ):
            Customized_Packet.CP(FITS_REF=refpath, FITS_SCI=scipath, FITS_mREF=ref_masked_savepath, FITS_mSCI=sci_masked_savepath, \
                                 ForceConv=ForceConv, GKerHW=GKerHW, FITS_DIFF=diff_savepath, FITS_Solution=soln_savepath, \
                                 KerPolyOrder=KerPolyOrder, BGPolyOrder=BGPolyOrder, ConstPhotRatio=ConstPhotRatio, \
                                 BACKEND_4SUBTRACT=backend, CUDA_DEVICE_4SUBTRACT=cudadevice, \
                                 NUM_CPU_THREADS_4SUBTRACT=nCPUthreads,logger=logger)

        size,peak = tracemalloc.get_traced_memory()
        logger.debug(f'MEMORY IN imagesubtraction.difference() AFTER Customized_Packet.CP: size = {size}, peak = {peak}')
        tracemalloc.reset_peak()

    elif skip_subtract and verbose:
        print(diff_savepath, 'already exists. Skipping image subtraction.')

    return diff_savepath, soln_savepath

def decorr_kernel(scipath, refpath, 
                  scipsfpath, refpsfpath,
                  diffpath, solnpath, out_path=output_files_rootdir, savename=None):

    savedir = os.path.join(out_path, 'dcker')
    check_and_mkdir(savedir)

    if savename is None:
        basename = os.path.basename(scipath)
        savename = F'DCKer_{basename}'

    decorr_savepath = os.path.join(savedir, savename)

    imgdatas = []
    psfdatas = []
    bkgsigs = []
    with nvtx.annotate( "get_data_and_SkyLevel", color="#ff44ff" ):
        for img, psf in zip([scipath, refpath], [scipsfpath, refpsfpath]):
            imgdata = fits.getdata(img, ext=0).T
            psfdatas.append(fits.getdata(psf, ext=0).T)
            imgdatas.append(imgdata)
            with nvtx.annotate( "SkyLevel_Estimator", color="#ff22ff" ):
                # bkgsigs.append(SkyLevel_Estimator.SLE(PixA_obj=imgdata)[1])
                # TODO : clean up interface.  This file was written in
                #   diff-img preprocess.py
                skyrmspath = os.path.join(output_files_rootdir, f'skyrms/{os.path.basename(img)}.json')
                bkgsigs.append( json.load( open( skyrmspath ) )['skyrms'] )
                

    sci_img, ref_img = imgdatas
    sci_psf, ref_psf = psfdatas
    sci_bkg, ref_bkg = bkgsigs

    with nvtx.annotate( "Realize_MatchingKernel", color="#ff88ff" ):
        N0, N1 = sci_img.shape
        XY_q = np.array([[N0/2. + 0.5, N1/2. + 0.5]])
        MKerStack = Realize_MatchingKernel(XY_q).FromFITS(FITS_Solution=solnpath)
        MK_Fin = MKerStack[0]

    with nvtx.annotate( "DeCorrelation_Calculator", color="#ffaaff" ):
        DCKer = DeCorrelation_Calculator.DCC(MK_JLst=[ref_psf], SkySig_JLst=[sci_bkg],
                                             MK_ILst=[sci_psf], SkySig_ILst=[ref_bkg], 
                                             MK_Fin=MK_Fin, KERatio=2.0, VERBOSE_LEVEL=2)

    with fits.open(scipath) as hdu:
        hdu[0].data = DCKer.T
        hdu.writeto(decorr_savepath, overwrite=True)

    return decorr_savepath

def decorr_img(imgpath, dckerpath, out_path=output_files_rootdir, savename=None):

    savedir = os.path.join(out_path,'decorr')
    check_and_mkdir(savedir)

    if savename is None:
        basename = os.path.basename(imgpath)
        savename = f'decorr_{basename}'

    decorr_savepath = os.path.join(savedir,savename)

    img_data = fits.getdata(imgpath, ext=0).T
    DCKer = fits.getdata(dckerpath)

    # Final decorrelated difference image:
    #  TODO : right now, we're acting as if
    #     the fftn= and ifftn= arguments
    #     of convolve_fft need to take and
    #     return numpy arrays, so we do all
    #     the conversions manually.  There
    #     must be a cleaner way.
    # THINK: can we just replace astropy convolve_fft
    #    with cupy and do our own post-processing?
    dcdiff = convolve_fft(img_data, DCKer, boundary='fill',
                          nan_treatment='fill', fill_value=0.0, 
                          normalize_kernel=True, preserve_nan=True,
                          fftn=lambda x: cupy.asnumpy( cupyx.scipy.fftpack.fftn( cupy.array( x ) ) ),
                          ifftn=lambda x: cupy.asnumpy( cupyx.scipy.fftpack.ifftn( cupy.array( x ) ) )
                          )

    with fits.open(imgpath) as hdu:
        hdu[0].data[:, :] = dcdiff.T
        hdu.writeto(decorr_savepath, overwrite=True)

    return decorr_savepath

def swarp_coadd(imgpath_list,refpath,out_name,out_path=output_files_rootdir,subdir='coadd',**kwargs):
    """Coadd images using SWarp. 

    kwargs: see sfft.utils.pyAstroMatic.PYSWarp.PY_SWarp.Mk_ConfigDict

    :param imgpath_list: Paths to images that will be coadded.
    :type imgpath_list: list
    :param refpath: Path to image to use as WCS reference.
    :type refpath: str
    :param savepath: Path to save coadded image.
    :type savepath: str
    :return: 
    :rtype: str
    """
    cd = PY_SWarp.Mk_ConfigDict(GAIN_DEFAULT=1., SATLEV_DEFAULT=100000., 
                                RESAMPLING_TYPE='BILINEAR', WEIGHT_TYPE='NONE', 
                                RESCALE_WEIGHTS='N', NTHREADS=1, **kwargs)

    imgpaths = []
    for p in imgpath_list:
        if p[-3:] == '.gz':
            zip_savedir = os.path.join(out_path, 'unzip')
            check_and_mkdir(zip_savedir)

            decompressed_path = os.path.join(zip_savedir,f'{os.path.basename(p)[:-3]}')
            dc_path = gz_and_ext(p, decompressed_path)

            imgpaths.append(dc_path)
        else:
            imgpaths.append(p)

    coadd_savedir = os.path.join(out_path,subdir)
    check_and_mkdir(coadd_savedir)
    coadd_savepath = os.path.join(coadd_savedir,out_name)
    
    coadd = PY_SWarp.Coadd(FITS_obj=imgpaths, FITS_ref=refpath, ConfigDict=cd,
                           OUT_path=coadd_savepath, FILL_VALUE=np.nan)

    return coadd_savepath, imgpaths

class imsub():
    """
    Class that wraps the functions in imagesubtraction.py and includes 
    path tracking. 

    USAGE EXAMPLE:
    from phrosty.imagesubtraction import imsub

    refpath = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/images/simple_model/R062/6/Roman_TDS_simple_model_R062_6_17.fits.gz'
    scipath = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/images/simple_model/R062/35083/Roman_TDS_simple_model_R062_35083_8.fits.gz'

    s = imsub(20172782,'R062','/work/lna18/imsub_out')
    s.set_template(refpath)
    s.set_science(scipath)
    s.preprocess_template()
    s.preprocess_sci()
    s.set_template_psf()
    s.set_science_psf()
    s.crossconv()
    s.diff()
    """
    def __init__(self,band,oid=None,ra=None,dec=None,
                    science_path=None,template_path=None,
                    sci_pointing=None,sci_sca=None,
                    temp_pointing=None,temp_sca=None,
                    sci_skysub_path=None,sci_align_path=None,
                    temp_psf_path=None,sci_psf_path=None,
                    cc_sci_path=None,cc_temp_path=None,
                    diff_path=None,soln_path=None,
                    decorr_kernel_path=None):
        self.oid = oid
        self.ra, self.dec = ra, dec
        self.start, self.end = None, None
        self.coord = None

        self.science_path = science_path
        self.template_path = template_path

        self.band = band
        self.sci_pointing = sci_pointing
        self.sci_sca = sci_sca
        self.temp_pointing = temp_pointing
        self.temp_sca = temp_sca

        self.sci_skysub_path = sci_skysub_path
        self.sci_align_path = sci_align_path

        self.temp_psf_path = temp_psf_path
        self.sci_psf_path = sci_psf_path

        self.cc_sci_path = cc_sci_path
        self.cc_temp_path = cc_temp_path

        self.diff_path = diff_path
        self.soln_path = soln_path

        self.decorr_kernel_path = decorr_kernel_path

    @property
    def get_radec(self):
        """Get the RA and Dec of the SN if object ID is provided."""
        self.ra, self.dec = get_transient_radec(self.oid)
        return self.ra, self.dec

    @property
    def sncoord(self):
        """Make a SkyCoord object with self.ra, self.dec."""
        assert (self.ra is not None or self.dec is not None), 'self.ra and/or self.dec are/is missing. \
                                                                Input RA and dec when initializing object, \
                                                                or use self.get_radec and provide object ID.'
        self.coord = SkyCoord(ra=self.ra*u.deg,dec=self.dec*u.deg)
        return self.coord
    
    @property
    def get_mjd(self):
        """Get the start and end MJD of the transient if object ID is provided."""
        self.start, self.end = get_transient_mjd(self.oid)
        return self.start, self.end

    def set_template(self,path):
        """Set the path to the template image."""
        assert os.path.exists(path), 'The input path does not exist.'

        if path[-3:] == '.gz':
            print('The input path needs to be .fits, not .fits.gz. Fixing this now.')
            zip_savedir = os.path.join(output_files_rootdir,'unzip')
            check_and_mkdir(zip_savedir)
            out_path = os.path.join(zip_savedir,f'{os.path.basename(path)[:-3]}')
            if not os.path.exists(out_path):
                gz_and_ext(path,out_path)
            else:
                print('Luckily,', out_path, 'already exists!')
            
            self.template_path = out_path
        else:
            self.template_path = path

        splitpath = path.split('_')
        self.temp_pointing = int(splitpath[-2])
        self.temp_sca = int(splitpath[-1].split('.')[0])

        return self.template_path

    def set_science(self,path):
        """Set the path to the science image."""
        assert os.path.exists(path), 'The input path does not exist.'

        if path[-3:] == '.gz':
            print('The input path needs to be .fits, not .fits.gz. Fixing this now.')
            zip_savedir = os.path.join(output_files_rootdir,'unzip')
            check_and_mkdir(zip_savedir)
            out_path = os.path.join(zip_savedir,f'{os.path.basename(path)[:-3]}')
            if not os.path.exists(out_path):
                gz_and_ext(path,out_path)
            else:
                print('Luckily,', out_path, 'already exists!')
            
            self.science_path = out_path
        else:
            self.science_path = path

        splitpath = path.split('_')
        self.sci_pointing = int(splitpath[-2])
        self.sci_sca = int(splitpath[-1].split('.')[0])

        return self.science_path

    def preprocess_template(self,skysub_path=None,align_path=None,force=False):
        """
        Sky subtract and 'align' the template image (to itself).
        Updates self.template_path to the sky-subtracted and 
        'aligned' template image rather than the original input.
        Optional step. 
        """
        if skysub_path is None:
            skysub_path = sky_subtract(path=self.template_path,force=force)

        if align_path is None:
            align_path = imalign(template_path=skysub_path,sci_path=skysub_path,force=force)

        self.template_path = align_path

        return self.template_path

    def preprocess_sci(self,skysub_path=None,align_path=None,force=False):
        assert self.template_path is not None, 'You need to set a template path. Go back and use self.set_template.'
        
        if skysub_path is None:
            print('Path to sky subtracted image not provided. Doing this now.')
            skysub_path = sky_subtract(path=self.science_path,force=force)
        if align_path is None:
            print('Path to aligned science image not provided. Doing this now. ')
            align_path = imalign(template_path=self.template_path,sci_path=skysub_path,force=force)

        self.sci_skysub_path = skysub_path
        self.sci_align_path = align_path

        return self.sci_skysub_path, self.sci_align_path

    def set_template_psf(self,path=None,force=False):
        """Set path to template psf."""

        if path is None:
            print('Path to template PSF not provided. Getting this now.')
            psf_path = rotate_psf(self.ra,self.dec,
                                  self.template_path,
                                  self.template_path,
                                  self.band,
                                  self.temp_pointing,
                                  self.temp_sca,
                                  ref_path=self.template_path,
                                  force=force)

            self.temp_psf_path = psf_path
        else:
            assert os.path.exists(path), f'Input path {path} does not exist.'
            self.temp_psf_path = path
        return self.temp_psf_path

    def set_science_psf(self,path=None,force=False):
        """Set path to science image PSF."""
        if path is None:
            print('Path to science PSF not provided. Getting this now.')
            psf_path = rotate_psf(self.ra,self.dec,
                                  self.sci_skysub_path,
                                  self.sci_align_path,
                                  self.band,
                                  self.sci_pointing,
                                  self.sci_sca,
                                  ref_path=self.template_path,
                                  force=force)

            self.sci_psf_path = psf_path
        else:
            assert os.path.exists(path), f'Input path {path} does not exist.'
            self.sci_psf_path = path
        return self.sci_psf_path

    def crossconv(self,force=False):
        """"Cross convolve the images with the PSFs.
        
        NOTE: Need to make 'force' arg actually do something. 
        """
        assert self.sci_align_path is not None, 'The path to the aligned science image does not exist. Go back and use self.preprocess_sci.'
        assert self.sci_psf_path is not None, 'The path to the science PSF does not exist. Go back and use self.set_science_psf.'
        assert self.template_path is not None, 'The path to the template image does not exist. Go back and use self.set_template.'
        assert self.temp_psf_path is not None, 'The path to the template PSF does not exist. Go back and use self.set_template_psf.'
        
        sci,ref = crossconvolve(self.sci_align_path, self.sci_psf_path,
                                self.template_path, self.temp_psf_path)

        self.cc_sci_path,self.cc_ref_path = sci,ref

        return self.cc_sci_path, self.cc_ref_path

    def diff(self,diffpath=None,solnpath=None,force=False):
        """Do the difference imaging."""
        
        if diffpath is None or solnpath is None:
            print('Path to difference image and/or solution not provided. Doing this now.')
            diffpath, solnpath = difference(self.cc_sci_path, self.cc_ref_path, 
                                    self.sci_psf_path, self.temp_psf_path, 
                                    backend='Numpy')
            self.diff_path,self.soln_path = diffpath,solnpath

        return self.diff_path,self.soln_path

    def show_diff(self,savepath=None):
        """Quick display and/or save the difference image."""
        assert self.diff_path is not None, 'The path to the difference image does not exist. Go back and use self.diff.'

        with fits.open(self.diff_path) as hdu:
            img = hdu[0].data
        vmin,vmax = ZScaleInterval().get_limits(img)
        plt.imshow(img,vmin=vmin,vmax=vmax)
        plt.colorbar()

        if savepath is not None:
            plt.savefig(savepath,bbox_inches='tight',dpi=300)
        else:
            plt.show()

    def set_dcker(self,dckerpath=None):
        """Make the decorrelation kernel."""

        if dckerpath is not None:
            assert os.path.exists(dckerpath), 'The path to the decorrelation kernel does not exist.'

        elif dckerpath is None:
            inputs = [self.cc_sci_path,self.cc_temp_path,
                                        self.sci_psf_path,self.temp_psf_path,
                                        self.diff_path,self.soln_path]
            
            assert None not in inputs, 'You need to define the paths to one or more of the following: \
                                    cross-convolved science image, cross-convolved template image, \
                                    science PSF, template PSF, difference, solution.'
            dckerpath = decorr_kernel(self.cc_sci_path,self.cc_temp_path,
                                        self.sci_psf_path,self.temp_psf_path,
                                        self.diff_path,self.soln_path)

        self.decorr_kernel_path = dckerpath
        return self.decorr_kernel_path

    def decorr(self,imgpath,dckerpath):
        """Decorrelate an image."""
        assert os.path.exists(dckerpath), 'The path to the decorrelation kernel does not exist.'
        assert os.path.exists(imgpath), 'The path to the image to be decorrelated does not exist.'
    
        dc = decorr_img(imgpath,dckerpath)
        return dc
