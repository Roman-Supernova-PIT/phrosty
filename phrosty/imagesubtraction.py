# IMPORTS Standard:
import os, sys
import numpy as np
import gzip
import shutil
from tempfile import mkdtemp
from glob import glob

# IMPORTS Astro:
from astropy.io import fits 
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
from sfft.utils.pyAstroMatic.PYSEx import PY_SEx
from sfft.CustomizedPacket import Customized_Packet
from sfft.utils.SkyLevelEstimator import SkyLevel_Estimator
from sfft.utils.SFFTSolutionReader import Realize_MatchingKernel
from sfft.utils.DeCorrelationCalculator import DeCorrelation_Calculator

# IMPORTS Internal:
from .utils import _build_filepath

"""
This module was written with significant contributions from 
Dr. Lei Hu (https://github.com/thomasvrussell/), and relies on his 
SFFT image subtraction package (https://github.com/thomasvrussell/sfft).
"""

output_files_rootdir = '/work/lna18/imsub_out/'

def check_and_mkdir(dirname):
    """
    Utility function for checking if a directory exists, and if not,
    makes that directory. 
    """
    if not os.path.exists(dirname):
        os.mkdir(dirname)

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
    zip_savedir = os.path.join(out_path, 'unzip')
    check_and_mkdir(zip_savedir)

    decompressed_path = os.path.join(output_files_rootdir,'unzip',f'{os.path.basename(original_imgpath)[:-3]}')
    gz_and_ext(original_imgpath, decompressed_path)
    output_path = os.path.join(out_path, 'skysub', f'skysub_{os.path.basename(decompressed_path)}')

    sub_savedir = os.path.join(out_path, 'skysub')
    check_and_mkdir(sub_savedir)

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
    check_and_mkdir(outdir)

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
    1. Get PSF model at specified RA, dec [degrees] from RomanDESC sims. 
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
    check_and_mkdir(psf_dir)

    psf_path = os.path.join(psf_dir, f'rot_psf_{ra}_{dec}_{sci_band}_{sci_pointing}_{sci_sca}.fits')
    fits.HDUList([fits.PrimaryHDU(data=psf_rotated, header=None)]).writeto(psf_path, overwrite=True)
    
    return psf_path

def crossconvolve(sci_img_path, sci_psf_path,
                    ref_img_path, ref_psf_path):

    savedir = os.path.join(output_files_rootdir,'convolved')
    check_and_mkdir(savedir)

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

def stampmaker(ra, dec, imgpath, shape=np.array([1000,1000])):

    savedir = os.path.join(output_files_rootdir,'stamps')
    check_and_mkdir(savedir)

    savename = os.path.basename(imgpath)
    savepath = os.path.join(savedir,f'{ra}_{dec}_stamp_{savename}')

    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    with fits.open(imgpath) as hdu:
        wcs = WCS(hdu[0].header)
    
    x, y = skycoord_to_pixel(coord, wcs)
    pxradec = np.array([[x,y]])
    
    Stamp_Generator.SG(FITS_obj=imgpath, COORD=pxradec, COORD_TYPE='IMAGE', \
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
                            MDIR=None, VERBOSE_LEVEL=1)[1][0]

    bkg_mask = (scatalog == 0)

    return bkg_mask

def sfft(scipath, refpath, 
        scipsfpath, refpsfpath, ForceConv='REF', GKerHW=9, KerPolyOrder=3, BGPolyOrder=0, 
        ConstPhotRatio=True, backend='Numpy', cudadevice='0', nCPUthreads=8):

    sci_basename = os.path.basename(scipath)
    ref_basename = os.path.basename(refpath)

    savedir = os.path.join(output_files_rootdir, 'subtract')
    check_and_mkdir(savedir)

    diff_savedir = os.path.join(savedir,'difference')
    soln_savedir = os.path.join(savedir, 'solution')
    masked_savedir = os.path.join(savedir, 'masked')

    for dirname in [diff_savedir, soln_savedir, masked_savedir]:
        check_and_mkdir(dirname)

    diff_savepath = os.path.join(diff_savedir, f'diff_{sci_basename}')
    soln_savepath = os.path.join(soln_savedir, f'solution_{sci_basename}')

    sci_masked_savepath = os.path.join(masked_savedir,f'masked_{sci_basename}')
    ref_masked_savepath = os.path.join(masked_savedir,f'masked_{ref_basename}')

    # Make combined detection mask.
    sci_bkgmask = bkg_mask(scipath)
    ref_bkgmask = bkg_mask(refpath)
    bkgmask = np.logical_and(sci_bkgmask,ref_bkgmask)


    for path, msavepath in zip([refpath, scipath], \
                                [ref_masked_savepath, sci_masked_savepath]):
        with fits.open(path) as hdu:
            hdudata = hdu[0].data.T
            hdudata[bkgmask] = 0.0
            hdu[0].data[:, :] = hdudata.T
            hdu.writeto(msavepath, overwrite=True)

    # Do SFFT subtraction
    Customized_Packet.CP(FITS_REF=refpath, FITS_SCI=scipath, FITS_mREF=ref_masked_savepath, FITS_mSCI=sci_masked_savepath, \
                        ForceConv=ForceConv, GKerHW=GKerHW, FITS_DIFF=diff_savepath, FITS_Solution=soln_savepath, \
                        KerPolyOrder=KerPolyOrder, BGPolyOrder=BGPolyOrder, ConstPhotRatio=ConstPhotRatio, \
                        BACKEND_4SUBTRACT=backend, CUDA_DEVICE_4SUBTRACT=cudadevice, \
                        NUM_CPU_THREADS_4SUBTRACT=nCPUthreads)

    return diff_savepath, soln_savepath

def decorr(scipath, refpath, 
            scipsfpath, refpsfpath,
            diffpath, solnpath):

    sci_basename = os.path.basename(scipath)
    ref_basename = os.path.basename(refpath)

    savedir = os.path.join(output_files_rootdir, 'subtract')
    check_and_mkdir(savedir)

    decorr_savedir = os.path.join(savedir, 'decorr')
    check_and_mkdir(decorr_savedir)

    decorr_savepath = os.path.join(decorr_savedir, f'decorr_{sci_basename}')

    imgdatas = []
    psfdatas = []
    bkgsigs = []
    for img, psf in zip([scipath, refpath], [scipsfpath, refpsfpath]):
        imgdata = fits.getdata(img, ext=0).T
        psfdatas.append(fits.getdata(psf, ext=0).T)
        imgdatas.append(imgdata)
        bkgsigs.append(SkyLevel_Estimator.SLE(PixA_obj=imgdata)[1])

    sci_img, ref_img = imgdatas
    sci_psf, ref_psf = psfdatas
    sci_bkg, ref_bkg = bkgsigs

    N0, N1 = sci_img.shape
    XY_q = np.array([[N0/2. + 0.5, N1/2. + 0.5]])
    MKerStack = Realize_MatchingKernel(XY_q).FromFITS(FITS_Solution=solnpath)
    MK_Fin = MKerStack[0]

    DCKer = DeCorrelation_Calculator.DCC(MK_JLst = [ref_psf], SkySig_JLst=[sci_bkg], \
                                        MK_ILst=[sci_psf], SkySig_ILst=[ref_bkg], MK_Fin=MK_Fin, \
                                        KERatio=2.0, VERBOSE_LEVEL=2)

    diff_data = fits.getdata(diffpath, ext=0).T
    dcdiff = convolve_fft(diff_data, DCKer, boundary='fill', \
                            nan_treatment='fill', fill_value=0.0, normalize_kernel=True)

    with fits.open(diffpath) as hdu:
        hdu[0].data[:, :] = dcdiff.T
        hdu.writeto(decorr_savepath, overwrite=True)

    return decorr_savepath