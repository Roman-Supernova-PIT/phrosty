# IMPORTS Standard:
import os
import numpy as np
import gzip
import shutil
import matplotlib.pyplot as plt
from tempfile import mkdtemp

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
from sfft.utils.pyAstroMatic.PYSEx import PY_SEx
from sfft.CustomizedPacket import Customized_Packet
from sfft.utils.SkyLevelEstimator import SkyLevel_Estimator
from sfft.utils.SFFTSolutionReader import Realize_MatchingKernel
from sfft.utils.DeCorrelationCalculator import DeCorrelation_Calculator
# from sfft.utils.meta.MultiProc import Multi_Proc

# IMPORTS Internal:
from .utils import _build_filepath, get_transient_radec, get_transient_mjd, get_fitsobj

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

    do_skysub = (force is True) or (force is False and not os.path.exists(output_path))
    skip_skysub = (not force) and os.path.exists(output_path)
    if do_skysub:
        sub_savedir = os.path.join(out_path, 'skysub')
        check_and_mkdir(sub_savedir)

        SEx_SkySubtract.SSS(FITS_obj=decompressed_path, FITS_skysub=output_path, FITS_sky=None, FITS_skyrms=None, \
                            ESATUR_KEY='ESATUR', BACK_SIZE=64, BACK_FILTERSIZE=3, DETECT_THRESH=1.5, \
                            DETECT_MINAREA=5, DETECT_MAXAREA=0, VERBOSE_LEVEL=2, MDIR=None)
    elif skip_skysub:
        print(output_path, 'already exists. Skipping sky subtraction.')

    return output_path

def imalign(template_path, sci_path, out_path=output_files_rootdir, force=False):
    """
    Align images with SWarp. 
    """
    outdir = os.path.join(out_path, 'align')
    output_path = os.path.join(outdir, f'align_{os.path.basename(sci_path)}')

    do_align = (force is True) or (force is False and not os.path.exists(output_path))
    skip_align = (not force) and os.path.exists(output_path)
    if do_align:
        check_and_mkdir(outdir)

        cd = PY_SWarp.Mk_ConfigDict(GAIN_KEY='GAIN', SATUR_KEY='SATURATE', OVERSAMPLING=1, RESAMPLING_TYPE='BILINEAR', \
                                    SUBTRACT_BACK='N', VERBOSE_TYPE='NORMAL', GAIN_DEFAULT=1., SATLEV_DEFAULT=100000.)
        PY_SWarp.PS(FITS_obj=sci_path, FITS_ref=template_path, ConfigDict=cd, FITS_resamp=output_path, \
                    FILL_VALUE=np.nan, VERBOSE_LEVEL=1, TMPDIR_ROOT=None)
    elif skip_align:
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

def get_imsim_psf(ra,dec,band,pointing,sca,size=201,out_path=output_files_rootdir,force=False):

    """
    Retrieve the PSF from roman_imsim/galsim, and transform the WCS so that CRPIX and CRVAL
    are centered on the image instead of at the corner. 

    force parameter does not currently do anything.
    """

    savedir = os.path.join(output_files_rootdir,'psf')
    check_and_mkdir(savedir)
    savename = f'psf_{ra}_{dec}_{band}_{pointing}_{sca}.fits'
    savepath = os.path.join(savedir,savename)

    # Get WCS of the image you need the PSF for.
    hdu = get_fitsobj(band=band,pointing=pointing,sca=sca)
    wcs = WCS(hdu[0].header)
    coord = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
    x,y = wcs.world_to_pixel(coord)

    # Get PSF at specified ra, dec. 
    config_path = os.path.join(os.path.dirname(__file__), 'auxiliary', 'tds.yaml')
    config = roman_utils(config_path,pointing,sca)
    psf = config.getPSF_Image(size,x,y)

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

def rotate_psf(ra,dec,psf,target,force=False,verbose=False):
    """
    2. Rotate PSF model to match reference WCS. 
        2a. Calculate rotation angle during alignment
        2b. Rotate PSF to match rotated science image
    """

    # Set up filepaths.
    psf_dir = os.path.join(output_files_rootdir,'psf')
    check_and_mkdir(psf_dir)

    basename = os.path.basename(psf)
    psf_path = os.path.join(psf_dir, f'rot_{basename}.fits')

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
        psfimg = fits.getdata(psf, ext=0) # Already saved as a transposed matrix from get_imsim_psf. 

        # Get vector from target WCS (i.e., rotated)
        hdr = fits.getheader(target, ext=0)
        _w = Read_WCS.RW(hdr, VERBOSE_LEVEL=1)
        x1, y1 = _w.all_world2pix(np.array([[ra0, dec0]]), 1)[0]
        skyN_vectorp = calculate_skyN_vector(wcshdr=hdr, x_start=x1, y_start=y1)
        PATTERN_ROTATE_ANGLE = calculate_rotate_angle(vector_ref=skyN_vector, vector_obj=skyN_vectorp)

        # Do rotation
        psf_rotated, psf_rot_head, func = Image_ZoomRotate.IZR(PixA_obj=psfimg, ZOOM_SCALE_X=1., \
                                            ZOOM_SCALE_Y=1., PATTERN_ROTATE_ANGLE=PATTERN_ROTATE_ANGLE, \
                                            RESAMPLING_TYPE='BILINEAR', FILL_VALUE=0.0, VERBOSE_LEVEL=1,
                                            RETURN_IMG_ONLY=False)

        # Save rotated PSF
        fits.HDUList([fits.PrimaryHDU(data=psf_rotated.T, header=psf_rot_head)]).writeto(psf_path, overwrite=True)
    elif skip_psf and verbose:
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

    do_conv = (force is True) or (force is False and not any([os.path.exists(p) for p in savepaths]))
    skip_conv =  (not force) and all([os.path.exists(p) for p in savepaths])
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

    elif skip_conv:
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

        # Do SFFT subtraction
        Customized_Packet.CP(FITS_REF=refpath, FITS_SCI=scipath, FITS_mREF=ref_masked_savepath, FITS_mSCI=sci_masked_savepath, \
                            ForceConv=ForceConv, GKerHW=GKerHW, FITS_DIFF=diff_savepath, FITS_Solution=soln_savepath, \
                            KerPolyOrder=KerPolyOrder, BGPolyOrder=BGPolyOrder, ConstPhotRatio=ConstPhotRatio, \
                            BACKEND_4SUBTRACT=backend, CUDA_DEVICE_4SUBTRACT=cudadevice, \
                            NUM_CPU_THREADS_4SUBTRACT=nCPUthreads)

    elif skip_subtract:
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
                                RESCALE_WEIGHTS='N', **kwargs)

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