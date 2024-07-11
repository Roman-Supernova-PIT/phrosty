# IMPORTS Standard:
import os
import numpy as np
import gzip
import shutil
import matplotlib.pyplot as plt

# IMPORTS Astro:
from astropy.io import fits 
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u
from astropy.wcs import WCS
from astropy.convolution import convolve_fft
from astropy.visualization import ZScaleInterval
# from photutils.psf import FittableImageModel

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
# from sfft.utils.meta.MultiProc import Multi_Proc

# IMPORTS Internal:
from .utils import _build_filepath, get_transient_radec, get_transient_mjd

"""
This module was written with significant contributions from 
Dr. Lei Hu (https://github.com/thomasvrussell/), and relies on his 
SFFT image subtraction package (https://github.com/thomasvrussell/sfft).
"""

output_files_rootdir = os.getenv('IMSUB_OUT','/work/lna18/imsub_out/')

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

def sky_subtract(path=None, band=None, pointing=None, sca=None, out_path=output_files_rootdir, force=False):

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

    if (force is True) or (force is False and not os.path.exists(output_path)):
        sub_savedir = os.path.join(out_path, 'skysub')
        check_and_mkdir(sub_savedir)

        SEx_SkySubtract.SSS(FITS_obj=decompressed_path, FITS_skysub=output_path, FITS_sky=None, FITS_skyrms=None, \
                            ESATUR_KEY='ESATUR', BACK_SIZE=64, BACK_FILTERSIZE=3, DETECT_THRESH=1.5, \
                            DETECT_MINAREA=5, DETECT_MAXAREA=0, VERBOSE_LEVEL=2, MDIR=None)
    elif not force and os.path.exists(output_path):
        print(output_path, 'already exists. Skipping sky subtraction.')

    return output_path

def imalign(template_path, sci_path, out_path=output_files_rootdir, force=False):
    """
    Align images with SWarp. 
    """
    outdir = os.path.join(out_path, 'align')
    output_path = os.path.join(outdir, f'align_{os.path.basename(sci_path)}')
    if (force is True) or (force is False and not os.path.exists(output_path)):
        check_and_mkdir(outdir)

        cd = PY_SWarp.Mk_ConfigDict(GAIN_KEY='GAIN', SATUR_KEY='SATURATE', OVERSAMPLING=1, RESAMPLING_TYPE='BILINEAR', \
                                    SUBTRACT_BACK='N', VERBOSE_TYPE='NORMAL', GAIN_DEFAULT=1., SATLEV_DEFAULT=100000.)
        PY_SWarp.PS(FITS_obj=sci_path, FITS_ref=template_path, ConfigDict=cd, FITS_resamp=output_path, \
                    FILL_VALUE=np.nan, VERBOSE_LEVEL=1, TMPDIR_ROOT=None)
    elif not force and os.path.exists(output_path):
        print(output_path, 'already exists. Skipping alignment.')

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

def get_imsim_psf(ra,dec,
                    sci_band,
                    sci_pointing,
                    sci_sca,
                    ref_band=None,
                    ref_pointing=None,
                    ref_sca=None,
                    ref_path=None):

    """
    Retrieves PSF directly from roman_imsim. If you have a reference image, retrieves it 
    according to the reference WCS. 
    """

    # Check if reference image was provided. If not, just retrieve PSF from science
    # image without changing the WCS. 
    if all(val is None for val in [ref_band,ref_pointing,ref_sca,ref_path]):
        print('Warning! No reference provided. WCS will not be rotated.')
        wcsband, wcspointing, wcssca, wcspath = sci_band, sci_pointing, sci_sca, ref_path
    else:
        wcsband, wcspointing, wcssca, wcspath = ref_band, ref_pointing, ref_sca, ref_path

    ref_path = _build_filepath(path=wcspath,band=wcsband,pointing=wcspointing,sca=wcssca,filetype='image')
    ref_hdu = fits.open(ref_path)
    ref_wcs = WCS(ref_hdu[0].header)
    worldcoords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    pxradec = skycoord_to_pixel(worldcoords,ref_wcs)

    # Get PSF at specified RA, dec in science image. 
    config_path = os.path.join(os.path.dirname(__file__), 'auxiliary', 'tds.yaml')
    config = roman_utils(config_path,visit=sci_pointing,sca=sci_sca)
    # add transpose here as Image_ZoomRotate only accept that form (LH, 2024/06/07)
    psf = config.getPSF_Image(501, x=pxradec[0], y=pxradec[1]).array.T 

    return psf

def rotate_psf(ra,dec,sci_skysub,
                    sci_imalign,
                    sci_band,
                    sci_pointing,
                    sci_sca,
                    ref_band=None,
                    ref_pointing=None,
                    ref_sca=None,
                    ref_path=None,
                    force=False):
    """
    2. Rotate PSF model to match reference WCS. 
        2a. Calculate rotation angle during alignment
        2b. Rotate PSF to match rotated science image
    """

    if all(val is None for val in [ref_band,ref_pointing,ref_sca,ref_path]):
        raise ValueError('You need to provide either [ref_band, ref_pointing, ref_sca] OR ref_path.')

    # Set up filepaths.
    psf_dir = os.path.join(output_files_rootdir,'psf')
    check_and_mkdir(psf_dir)

    psf_path = os.path.join(psf_dir, f'rot_psf_{ra}_{dec}_{sci_band}_{sci_pointing}_{sci_sca}.fits')

    if (force is True) or (force is False and not os.path.exists(psf_path)):
        psf = get_imsim_psf(ra,dec,
                            sci_band,sci_pointing,sci_sca,
                            ref_band,ref_pointing,ref_sca,ref_path)

        # Get vector from sky-subtracted science WCS (i.e., not rotated to reference)
        hdr = fits.getheader(sci_skysub, ext=0)
        _w = Read_WCS.RW(hdr, VERBOSE_LEVEL=1)
        x0, y0 = 0.5 + int(hdr['NAXIS1'])/2.0, 0.5 + int(hdr['NAXIS2'])/2.0
        ra0, dec0 = _w.all_pix2world(np.array([[x0, y0]]), 1)[0]
        skyN_vector = calculate_skyN_vector(wcshdr=hdr, x_start=x0, y_start=y0)

        # Get vector from rotated science WCS (i.e., rotated to reference)
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
        fits.HDUList([fits.PrimaryHDU(data=psf_rotated.T, header=None)]).writeto(psf_path, overwrite=True)
    elif not force and os.path.exists(psf_path):
        print(psf_path, 'already exists. Skipping getting PSF.')

    return psf_path

def crossconvolve(sci_img_path, sci_psf_path,
                    ref_img_path, ref_psf_path,
                    force=False):

    savedir = os.path.join(output_files_rootdir,'convolved')
    check_and_mkdir(savedir)

    # First convolves reference PSF on science image. 
    # Then, convolves science PSF on reference image. 
    savepaths = []
    for img in [sci_img_path,ref_img_path]:
        savename = f'conv_{os.path.basename(img)}'
        savepath = os.path.join(savedir, savename)
        savepaths.append(savepath)

    if (force is True) or (force is False and not any([os.path.exists(p) for p in savepaths])):
        for img, psf, name in zip([sci_img_path, ref_img_path],
                            [ref_psf_path, sci_psf_path],
                            ['sci', 'ref']):

            imgdata = fits.getdata(img, ext=0).T
            psfdata = fits.getdata(psf, ext=0).T

            convolved = convolve_fft(imgdata, psfdata, boundary='fill', nan_treatment='fill', \
                                    fill_value=0.0, normalize_kernel=True, preserve_nan=True, allow_huge=True)

            with fits.open(img) as hdl:
                hdl[0].data[:, :] = convolved.T
                hdl.writeto(savepath, overwrite=True)
                savepaths.append(savepath)
    elif not force and all([os.path.exists(p) for p in savepaths]):
        print(savepaths, 'already exist. Skipping cross convolve.')

    return savepaths

def stampmaker(ra, dec, imgpath, savepath=None, shape=np.array([1000,1000])):

    if savepath is None:
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

def difference(scipath, refpath, 
        scipsfpath, refpsfpath, ForceConv='REF', GKerHW=9, KerPolyOrder=3, BGPolyOrder=0, 
        ConstPhotRatio=True, backend='Numpy', cudadevice='0', nCPUthreads=8, force=False):

    sci_basename = os.path.basename(scipath)
    ref_basename = os.path.basename(refpath)

    sci_data = fits.getdata(scipath).T
    ref_data = fits.getdata(refpath).T

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


    if (force is True) or (force is False and not os.path.exists(diff_savepath)):
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

        # Do SFFT subtraction
        Customized_Packet.CP(FITS_REF=refpath, FITS_SCI=scipath, FITS_mREF=ref_masked_savepath, FITS_mSCI=sci_masked_savepath, \
                            ForceConv=ForceConv, GKerHW=GKerHW, FITS_DIFF=diff_savepath, FITS_Solution=soln_savepath, \
                            KerPolyOrder=KerPolyOrder, BGPolyOrder=BGPolyOrder, ConstPhotRatio=ConstPhotRatio, \
                            BACKEND_4SUBTRACT=backend, CUDA_DEVICE_4SUBTRACT=cudadevice, \
                            NUM_CPU_THREADS_4SUBTRACT=nCPUthreads)

    elif not force and os.path.exists(diff_savepath):
        print(diff_savepath, 'already exists. Skipping image subtraction.')

    return diff_savepath, soln_savepath

def decorr_kernel(scipath, refpath, 
            scipsfpath, refpsfpath,
            diffpath, solnpath):

    sci_basename = os.path.basename(scipath)

    savedir = os.path.join(output_files_rootdir, 'dcker')
    check_and_mkdir(savedir)

    decorr_savepath = os.path.join(savedir, f'DCKer_{sci_basename}')

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

    with fits.open(scipath) as hdu:
        hdu[0].data = DCKer.T
        hdu.writeto(decorr_savepath, overwrite=True)

    return decorr_savepath

def decorr_img(imgpath, dckerpath, imgtype='difference'):
    decorr_basename = os.path.basename(imgpath)
    if imgtype == 'difference':
        savedir = os.path.join(output_files_rootdir,'subtract','decorr')
        check_and_mkdir(savedir)
        decorr_savepath = os.path.join(savedir,f'decorr_{decorr_basename}')

    elif imgtype == 'science':
        savedir = os.path.join(output_files_rootdir,'science')
        check_and_mkdir(savedir)
        decorr_savepath = os.path.join(savedir,f'decorr_{decorr_basename}')

    img_data = fits.getdata(imgpath, ext=0).T
    DCKer = fits.getdata(dckerpath)

    # Final decorrelated difference image: 
    dcdiff = convolve_fft(img_data, DCKer, boundary='fill', \
                            nan_treatment='fill', fill_value=0.0, normalize_kernel=True,
                            preserve_nan=True)

    with fits.open(imgpath) as hdu:
        hdu[0].data[:, :] = dcdiff.T
        hdu.writeto(decorr_savepath, overwrite=True)

    return decorr_savepath

def calc_psf(scipath, refpath,
            scipsfpath, refpsfpath, 
            dckerpath,
            SUBTTAG='DCSCI', nproc=1, TILESIZE_RATIO=5, GKerHW=9):
    """
    Calculate the PSF. 
    scipath -- should be sky-subtracted, aligned, and cross-convolved with the reference PSF. 
    refpath -- should be sky-subtracted and cross-convolved with the science PSF. 
    """

    psf_basename = os.path.basename(scipath[:-5])
    savedir = os.path.join(output_files_rootdir,'psf_final')
    check_and_mkdir(savedir)
    psf_savepath = os.path.join(output_files_rootdir,'psf_final',f'{psf_basename}.sfft_{SUBTTAG}.DeCorrelated.dcPSFFStack.fits')

    # * define an image grid (use one psf size)
    _hdr = fits.getheader(scipath, ext=0) # FITS_lSCI = output_dir + '/%s.sciE.skysub.fits' %sciname
    # N0, N1 = int(_hdr['NAXIS1']), int(_hdr['NAXIS2'])

    # lab = 0
    # XY_TiC = []
    # TILESIZE_RATIO = 10
    # TiHW = round(TILESIZE_RATIO * GKerHW) 
    # TiN = 2*TiHW+1
    # AllocatedL = np.zeros((N0, N1), dtype=int)
    # for xs in np.arange(0, N0, TiN):
    #     xe = np.min([xs+TiN, N0])
    #     for ys in np.arange(0, N1, TiN):
    #         ye = np.min([ys+TiN, N1])
    #         AllocatedL[xs: xe, ys: ye] = lab
    #         x_q = 0.5 + xs + (xe - xs)/2.0   # tile-center (x)
    #         y_q = 0.5 + ys + (ye - ys)/2.0   # tile-center (y)
    #         XY_TiC.append([x_q, y_q])
    #         lab += 1
    # XY_TiC = np.array(XY_TiC)
    # NTILE = XY_TiC.shape[0]

    # XY_q = np.array([[N0/2.+0.5, N1/2.+0.5]])
    # MKerStack = Realize_MatchingKernel(XY_q).FromFITS(FITS_Solution=soln).T

    # PixA_lREF = fits.getdata(refpath, ext=0).T # use stamp
    # PixA_lSCI = fits.getdata(scipath, ext=0).T # use stamp

    # PixA_PSF_lREF = fits.getdata(refpsfpath, ext=0).T
    PixA_PSF_lSCI = fits.getdata(scipsfpath, ext=0).T

    # bkgsig_REF = SkyLevel_Estimator.SLE(PixA_obj=PixA_lREF)[1]
    # bkgsig_SCI = SkyLevel_Estimator.SLE(PixA_obj=PixA_lSCI)[1]

    # record model PSF of decorrelated image (DCMREF, DCSCI, DCDIFF) for the grid
    # FITS_dcPSFFStack = os.path.join(output_files_rootdir,'psf_final',f'{os.path.basename(scipath)[:-5]}.sfft_{SUBTTAG}.DeCorrelated.dcPSFFStack.fits')

    # XY_q = np.array([[N0/2.+0.5, N1/2.+0.5]])
    # MKerStack = Realize_MatchingKernel(XY_q).FromFITS(FITS_Solution=soln)
    # MK_Fin = MKerStack[0]

    # calculate decorrelation kernels on the grid
    # DCKer = DeCorrelation_Calculator.DCC(MK_JLst=[PixA_PSF_lREF], SkySig_JLst=[bkgsig_SCI], \
    #         MK_ILst=[PixA_PSF_lSCI], SkySig_ILst=[bkgsig_REF], MK_Fin=MK_Fin, \
    #         KERatio=2.0, VERBOSE_LEVEL=0)

    DCKer = fits.getdata(dckerpath)
    NX_DCKer, NY_DCKer = DCKer.shape

    PixA_dcPSF = convolve_fft(PixA_PSF_lSCI, DCKer, boundary='fill', \
            nan_treatment='fill', fill_value=0.0, normalize_kernel=True)

    _hdr = fits.Header()
    # _hdr['NTILE'] = NTILE
    _hdr['NX_PSF'] = PixA_dcPSF.shape[0]
    _hdr['NY_PSF'] = PixA_dcPSF.shape[1]
    _hdl = fits.HDUList([fits.PrimaryHDU(PixA_dcPSF.T, header=_hdr)])
    _hdl.writeto(psf_savepath, overwrite=True)
    print("MeLOn CheckPoint: PSF for DeCorrelated images Saved! \n # %s!" %psf_savepath)

    return psf_savepath

def swarp_coadd_img(imgpath_list,refpath,out_name,out_path=output_files_rootdir,**kwargs):
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
                                RESAMPLING_TYPE='BILINEAR', **kwargs)

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

    coadd_savedir = os.path.join(out_path,'coadd')
    check_and_mkdir(coadd_savedir)
    coadd_savepath = os.path.join(coadd_savedir,out_name)
    
    coadd = PY_SWarp.Coadd(FITS_obj=imgpaths, FITS_ref=refpath, ConfigDict=cd,
                            OUT_path=coadd_savepath, FILL_VALUE=np.nan)

    return coadd_savepath, imgpaths

def swarp_coadd_psf(ra,dec,sci_skysub_paths,sci_imalign_paths,
                    ref_table,refpath,out_name,out_path=output_files_rootdir):
    """_summary_

    :param ra: _description_
    :type ra: _type_
    :param dec: _description_
    :type dec: _type_
    :param sci_imalign: _description_
    :type sci_imalign: _type_
    :param sci_skysub: _description_
    :type sci_skysub: _type_
    :param ref_table: filter, pointing, sca astropy table w/ images in the coadded template. 
    :type ref_table: _type_
    :param refpath: Path to coadded template. For WCS info. 
    :type refpath: _type_
    :param out_name: _description_
    :type out_name: _type_
    :param out_path: _description_, defaults to output_files_rootdir
    :type out_path: _type_, optional
    :return: _description_
    :rtype: _type_
    """

    coadd_psf_savedir = os.path.join(out_path,'coadd_psf')
    check_and_mkdir(coadd_psf_savedir)
    coadd_psf_savepath = os.path.join(coadd_psf_savedir,out_name)

    psf_list = []
    # Can be parallelized: 
    for i, row in enumerate(ref_table):
        psf = rotate_psf(ra,dec,sci_skysub_paths[i],sci_imalign_paths[i],
                         sci_band=row['filter'],sci_pointing=row['pointing'],sci_sca=row['sca'],
                         ref_path=refpath)
        psf_list.append(psf)

    coadd_psf = swarp_coadd_img(psf_list,refpath,coadd_psf_savepath)

    return coadd_psf_savepath

class imsub():
    def __init__(self,oid,band,rootdir):
        self.ra, self.dec = get_transient_radec(oid)
        self.start, self.end = get_transient_mjd(oid)
        self.coord = SkyCoord(ra=self.ra*u.deg,dec=self.dec*u.deg)

        self.science_path = None
        self.template_path = None

        self.band = band
        self.sci_pointing = None
        self.sci_sca = None
        self.temp_pointing = None
        self.temp_sca = None

        self.sci_skysub_path = None
        self.sci_align_path = None

        self.temp_psf_path = None
        self.sci_psf_path = None

        self.cc_sci_path = None
        self.cc_temp_path = None

        self.diff_path = None
        self.soln_path = None

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
            align_path = imalign(template_path=self.template_path,sci_path=skysub,force=force)

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

    def show_diff(self):
        assert self.diff_path is not None, 'The path to the difference image does not exist. Go back and use self.diff.'

        with fits.open(self.diff_path) as hdu:
            img = hdu[0].data
        vmin,vmax = ZScaleInterval().get_limits(img)
        plt.imshow(img,vmin=vmin,vmax=vmax)
        plt.colorbar()
        plt.show()