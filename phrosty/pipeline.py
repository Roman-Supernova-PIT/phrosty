__all__ = [ 'PipelineImage', 'Pipeline' ]

# Imports STANDARD
import sys
import argparse
import cupy as cp
from functools import partial
import json
import logging
from multiprocessing import Pool
import numpy as np
import nvtx
import pathlib
import pyarrow as pa
import pyarrow.parquet as pq
import re
import shutil
import tracemalloc
import uuid

# Imports ASTRO
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
from astropy.wcs.utils import skycoord_to_pixel
import fitsio

# Imports INTERNAL
import phrosty
from phrosty.imagesubtraction import sky_subtract, stampmaker
from sfft.SpaceSFFTCupyFlow import SpaceSFFT_CupyFlow
from snappl.dbclient import SNPITDBClient
from snappl.diaobject import DiaObject
from snappl.imagecollection import ImageCollection
from snappl.image import FITSImageOnDisk
from snappl.lightcurve import Lightcurve
from snappl.provenance import Provenance
from snappl.psf import PSF
from snappl.config import Config
from snappl.logger import SNLogger


class PipelineImage:
    """Holds a snappl.image.Image, with some other stuff the pipeline needs."""

    def __init__( self, image, pipeline ):
        """Create a PipelineImage

        Parameters:
        -----------
           image : snappl.image.Image
              The image we're encapsulating.  Pass either this or imagepath.

           pipeline : phrosty.pipeline.Pipeline
              The pipeline that owns this image.

        """

        self.config = Config.get()
        self.temp_dir = pipeline.temp_dir
        self.keep_intermediate = self.config.value( 'photometry.phrosty.keep_intermediate' )
        if self.keep_intermediate:
            self.save_dir = pathlib.Path( self.config.value( 'system.paths.scratch_dir' ) )
        elif not self.keep_intermediate:
            self.save_dir = self.temp_dir

        self.image = image
        if self.image.band != pipeline.band:
            raise ValueError( f"Image {self.image.path.name} has a band {self.image.band}, "
                              f"which is different from the pipeline band {pipeline.band}" )

        # Intermediate files
        if self.keep_intermediate:
            # Set to None. The path gets defined later on.
            # They have to be defined here in __init__ so that they exist
            # and are accessible in later functions.
            self.skysub_img = None
            self.detmask_img = None
            self.input_sci_psf_path = None
            self.input_templ_psf_path = None
            self.aligned_templ_img_path = None
            self.aligned_templ_var_path = None
            self.aligned_templ_psf_path = None
            self.crossconv_sci_path = None
            self.crossconv_templ_path = None
            self.diff_path = None
            self.decorr_kernel_path = None

        # Always save and output these
        self.decorr_psf_path = {}
        self.decorr_zptimg_path = {}
        self.decorr_diff_path = {}
        self.zpt_stamp_path = {}
        self.diff_var_path = {}
        self.diff_var_stamp_path = {}
        self.diff_stamp_path = {}

        # Held in memory
        self.skyrms = None
        self.psfobj = None
        self.psf_data = None

    def run_sky_subtract( self, mp=True ):
        """Run sky subtraction using Source Extractor.

        Parameters
        ----------
        mp : bool, optional
            Toggle multiprocessing, by default True

        Returns
        -------
        tuple
            Tuple containing sky subtracted image, detection
            mask array, and sky RMS value.
            Output of phrosty.imagesubtraction.sky_subtract().
        """
        try:
            return sky_subtract( self.image )
        except Exception as ex:
            SNLogger.exception( ex )
            raise

    def save_sky_subtract_info( self, info ):
        """Saves the sky-subtracted image, detection mask array,
        and sky RMS values to attributes.

        Parameters
        ----------
        info : tuple
            Output of self.run_sky_subtract(). See documentation
            for phrosty.imagesubtraction.sky_subtract().
        """

        SNLogger.debug( f"Saving sky_subtract info for path {info[0]}" )
        self.skysub_img = info[0]
        self.detmask_img = info[1]
        self.skyrms = info[2]

    def get_psf( self, ra, dec ):
        """Get the at the right spot on the image.

        Parameters
        ----------
          ra, dec : float
             The coordinates in decimal degrees where we want the PSF.

        Returns
        -------
        np.array
            PSF stamp. If this function fails, None is returned.

        """

        # TODO: right now snappl.psf.PSF.get_psf_object just
        #   passes the keyword arguments on to whatever makes
        #   the psf... and it's different for each type of
        #   PSF.  We need to fix that... somehow....

        wcs = self.image.get_wcs()
        x, y = wcs.world_to_pixel( ra, dec )

        if self.psfobj is None:
            psftype = self.config.value( 'photometry.phrosty.psf.type' )
            psfparams = self.config.value( 'photometry.phrosty.psf.params' )
            self.psfobj = PSF.get_psf_object( psftype, x=x, y=y,
                                              band=self.image.band,
                                              pointing=self.image.pointing,
                                              sca=self.image.sca,
                                              **psfparams )

        stamp = self.psfobj.get_stamp( x, y )
        if self.keep_intermediate:
            outfile = self.save_dir / f"psf_{self.image.name}.fits"
            fitsio.write( outfile, stamp, clobber=True )

        return stamp

    def keep_psf_data( self, psf_data ):
        """Save PSF data to attribute.

        Parameters
        ----------
        psf_data : np.array
            PSF stamp.

        """

        self.psf_data = psf_data

    def free( self, free_psf_data=False ):
        """Try to free memory.  More might be done here."""
        self.image.free()
        self.skysub_img.free()
        self.detmask_img.free()
        if free_psf_data:
            self.psf_data = None


class Pipeline:
    """Phrosty's top-level pipeline"""

    def __init__( self, diaobj, imgcol, band,
                  science_images=None,
                  template_images=None,
                  science_csv=None,
                  template_csv=None,
                  oid=None,
                  ltcv_prov_tag=None,
                  dbsave=False,
                  dbclient=None,
                  nprocs=1,
                  nwrite=5,
                  verbose=False,
                  memtrace=False ):

        """Create the a pipeline object.

        Parameters
        ----------
           diaobj : DiaObject
             The object we're building a lightcurve for

           imgcol : ImageCollection
             snappl.imagecollection.ImageCollection

           band: str
             One of R062, Z087, Y106, J129, H158, F184, K213

           science_images: list of snappl.image.Image
             The science images.

           template_images: list of snappl.image.Image
             The template images.

           science_csv: Path or str
             CSV file with the science images.  The first line must be::

               path pointing sca mjd band

             subsequent lines must have that information for all the
             science images.  path must be relative to ou24.images in
             config.  Pipeline will extract the images from this file
             whose band matches the band of the pipeline (and ignore
             the rest)

           template_csv: Path or str
             CSV file with template images.  Same format as science_csv.

           oid: str
             Object ID. This is probably a temporary argument. It is only used
             if diaobj is None, and is really only used to build a filepath
             without "None" in it when saving the light curve parquet file.

           ltcv_prov_tag: str
             Provenance tag for light curve. Required to use SN PIT database.

           dbsave: bool
             Are we saving to the database?
             Default False.

           dbclient: snappl.dbclient.SNPITDBClient

           nprocs: int, default 1
             Number of cpus for the CPU multiprocessing segments of the pipeline.
             (GPU segments will run a single process.)

           nwrite: int, default 5
             Number of asynchronous FITS writer processes.

           verbose: bool, default True
             Toggle verbose output.

           memtrace: bool, default False
             Toggle memory tracing.

        """

        SNLogger.setLevel( logging.DEBUG if verbose else logging.INFO )

        self.config = Config.get()
        self.imgcol = imgcol
        self.diaobj = diaobj
        self.band = band
        self.oid = oid

        self.dia_out_dir = pathlib.Path( self.config.value( 'system.paths.dia_out_dir' ) )
        self.scratch_dir = pathlib.Path( self.config.value( 'system.paths.scratch_dir' ) )
        self.temp_dir_parent = pathlib.Path( self.config.value( 'system.paths.temp_dir' ) )
        self.temp_dir = self.temp_dir_parent / str(uuid.uuid1())
        self.temp_dir.mkdir()
        self.ltcv_dir = pathlib.Path( self.config.value( 'system.paths.ltcv_dir' ) )

        if ( science_images is None) == ( science_csv is None ):
            raise ValueError( "Pass exactly one of science_images or science_csv" )
        if science_csv is not None:
            science_images = self._read_csv( science_csv )

        if ( template_images is None ) == ( template_csv is None ):
            raise ValueError( "Pass exactly one of template_images or template_csv" )
        if template_csv is not None:
            template_images = self._read_csv( template_csv )

        self.science_images = [ PipelineImage( i, self ) for i in science_images ]
        self.template_images = [ PipelineImage( i, self ) for i in template_images ]

        # All of our failures.
        # LA NOTE: I want to add keys for when this fails in sky subtraction and PSF retrieval
        # as well. But that is slightly more involved because those things happen in
        # PipelineImage, not Pipeline. So, for now, anything that fails at 'make_lightcurve'
        # may have failed at earlier steps first. Also, source extractor doesn't fail on
        # an image full of NaNs for some reason. It just reports 0 sources. I would expect
        # that it should fail, here...
        self.failures = {'align_and_preconvolve': [],
                         'find_decorrelation': [],
                         'subtract': [],
                         'variance': [],
                         'apply_decorrelation': [],
                         'make_lightcurve': [],
                         'make_stamps': []}

        self.ltcv_prov_tag = ltcv_prov_tag
        self.dbsave = dbsave
        self.dbclient = dbclient
        self.nprocs = nprocs
        self.nwrite = nwrite

        self.keep_intermediate = self.config.value( 'photometry.phrosty.keep_intermediate' )
        self.remove_temp_dir = self.config.value( 'photometry.phrosty.remove_temp_dir' )
        self.mem_trace = memtrace

    def _read_csv( self, csvfile ):
        """Reads input csv files with columns:
        'path pointing sca mjd band'.

        Parameters
        ----------
        csvfile : str
            Path to an input csv file.

        Returns
        -------
        list of snappl.image.Image

        Raises
        ------
        ValueError
            If the first line of the csv file doesn't match
            'path pointing sca mjd band', a ValueError is raised.
        """

        imlist = []
        with open( csvfile ) as ifp:
            hdrline = ifp.readline()
            if not re.search( r"^\s*path\s+pointing\s+sca\s+mjd\s+band\s*$", hdrline ):
                raise ValueError( f"First line of list file {csvfile} didn't match what was expected." )
            for line in ifp:
                path, pointing, sca, _mjd, band = line.split()
                if band == self.band:
                    # This should yell at us if the pointing or sca doesn't match what is read from the path
                    imlist.append( self.imgcol.get_image( path=path, pointing=pointing, sca=sca, band=band ) )
        return imlist


    def sky_sub_all_images( self ):
        """Sky subtracts all snappl.image.Image objects in
        self.science_images and self.template_images using
        Source Extractor.

        Contains its own error logging function, log_error().

        """

        # Currently, this writes out a bunch of FITS files.  Further refactoring needed
        #   to support more general image types.
        all_imgs = self.science_images.copy()     # shallow copy
        all_imgs.extend( self.template_images )

        def log_error( img, x ):
            SNLogger.error( f"Sky subtraction failure on {img.image.path}: {x}" )

        if self.nprocs > 1:
            with Pool( self.nprocs ) as pool:
                for img in all_imgs:
                    pool.apply_async( img.run_sky_subtract, (), {},
                                      callback=img.save_sky_subtract_info,
                                      error_callback=partial(log_error, img) )
                pool.close()
                pool.join()
        else:
            for img in all_imgs:
                img.save_sky_subtract_info( img.run_sky_subtract( mp=False ) )

    def get_psfs( self ):
        """Retrieve PSFs for all snappl.image.Image objects in
        self.science_images and self.template_images.

        Contains its own error logging function, log_error().

        """

        all_imgs = self.science_images.copy()     # shallow copy
        all_imgs.extend( self.template_images )

        def log_error( img, x ):
            SNLogger.error( f"get_psf failure on {img.image.path}: {x}" )

        if self.nprocs > 1:
            with Pool( self.nprocs ) as pool:
                for img in all_imgs:
                    # callback_partial = partial( img.save_psf_path, all_imgs )
                    pool.apply_async( img.get_psf, (self.diaobj.ra, self.diaobj.dec), {},
                                      callback=img.keep_psf_data,
                                      error_callback=partial(log_error, img) )
                pool.close()
                pool.join()
        else:
            for img in all_imgs:
                img.keep_psf_data( img.get_psf(self.diaobj.ra, self.diaobj.dec) )

    def align_and_pre_convolve(self, templ_image, sci_image ):
        """Align and pre convolve a single template/science pair.

        Parameters
        ----------
          sci_image: phrosty.PipelineImage
            The science (new) image.

          templ_image: phrosty.PipelineImage
            The template (ref) image that will be subtracted from sci_image.

        Returns
        -------
          sfftifier: SpaceSFFT_CupyFlow
            Use this object for futher SFFT work.  Be sure to
            dereference it to free the prodigious amount of memory it
            allcoates.

        """

        # SFFT needs FITS headers with a WCS and with NAXIS[12]
        hdr_sci = sci_image.image.get_wcs().get_astropy_wcs().to_header( relax=True )
        hdr_sci.insert( 0, ('NAXIS', 2) )
        hdr_sci.insert( 'NAXIS', ('NAXIS1', sci_image.image.data.shape[1] ), after=True )
        hdr_sci.insert( 'NAXIS1', ('NAXIS2', sci_image.image.data.shape[0] ), after=True )
        data_sci = cp.array( np.ascontiguousarray(sci_image.image.data.T), dtype=cp.float64 )
        noise_sci = cp.array( np.ascontiguousarray(sci_image.image.noise.T), dtype=cp.float64 )

        hdr_templ = templ_image.image.get_wcs().get_astropy_wcs().to_header( relax=True )
        hdr_templ.insert( 0, ('NAXIS', 2) )
        hdr_templ.insert( 'NAXIS', ('NAXIS1', templ_image.image.data.shape[1] ), after=True )
        hdr_templ.insert( 'NAXIS1', ('NAXIS2', templ_image.image.data.shape[0] ), after=True )
        data_templ = cp.array( np.ascontiguousarray(templ_image.image.data.T), dtype=cp.float64 )
        noise_templ = cp.array( np.ascontiguousarray(templ_image.image.noise.T), dtype=cp.float64 )

        sci_psf = cp.ascontiguousarray( cp.array( sci_image.psf_data.T, dtype=cp.float64 ) )
        templ_psf = cp.ascontiguousarray( cp.array( templ_image.psf_data.T, dtype=cp.float64 ) )

        sci_detmask = cp.array( np.ascontiguousarray( sci_image.detmask_img.data.T ) )
        templ_detmask = cp.array( np.ascontiguousarray( templ_image.detmask_img.data.T ) )

        sfftifier = SpaceSFFT_CupyFlow(
            hdr_sci, hdr_templ,
            sci_image.skyrms, templ_image.skyrms,
            data_sci, data_templ,
            noise_sci, noise_templ,
            sci_detmask, templ_detmask,
            sci_psf, templ_psf,
            KerPolyOrder=Config.get().value('photometry.phrosty.kerpolyorder')
        )

        sfftifier.resampling_image_mask_psf()
        sfftifier.cross_convolution()

        return sfftifier

    def phot_at_coords( self, img, psf, pxcoords=(50, 50), ap_r=4 ):
        """Do photometry at forced set of pixel coordinates.

        Parameters
        ----------
          img: snappl.image.Image
            The image on which to do the photometry

          psf: snappl.psf.PSF
            The PSF.

          pxcoords: tuple of (int, int)
            The position on the image to do the photometry

          ap_r: float
            Radius of aperture.

        Returns
        -------
          results: dict
            Keys and values are:
            * 'aperture_sum': flux in aperture of radius ap_r
            * 'flux_fit': flux from PSF photometry
            * 'flux_fit_err': uncertainty on flux_fit
            * 'mag_fit': instrumental magnitude (i.e. no zeropoint) from flux_fit
            * 'mag_fit_err': uncertainty on mag_fit

            All values are floats.

        """

        forcecoords = Table([[float(pxcoords[0])], [float(pxcoords[1])]], names=["x", "y"])
        init = img.ap_phot( forcecoords, ap_r=ap_r )
        init.rename_column( 'aperture_sum', 'flux_init' )
        init.rename_column( 'xcenter', 'xcentroid' )
        init.rename_column( 'ycenter', 'ycentroid' )
        final = img.psf_phot( init, psf, forced_phot=True )

        flux = final['flux_fit'][0]
        flux_err = final['flux_err'][0]
        mag = -2.5 * np.log10(final["flux_fit"][0])
        mag_err = (2.5 / np.log(10)) * np.abs(final["flux_err"][0] / final["flux_fit"][0])

        results_dict = {
                        'flux': flux,
                        'flux_err': flux_err,
                        'aperture_sum': init['flux_init'][0],  # Has to be renamed and named back because photutils.
                        'mag': mag,
                        'mag_err': mag_err
                        }

        return results_dict

    def make_phot_info_dict( self, sci_image, templ_image, ap_r=4 ):
        """"
        Do photometry on a difference image generated from sci_image
        and templ_image. Collect the output in a dictionary.

        Parmaeters
        ----------
          sci_image: PipelineImage
            science image wrapper

          temp_image: PipelineImage
            template image wrapper

          ap_r: float, default 4
             Radius of aperture to use in aperture photometry.

        Returns
        -------
          results_dict: dict
            Dictionary with keys sci_name, templ_name, success, ra,
            dec, mjd, filter, pointing, sca, template_pointing,
            template_sca, zpt, aperture_sum, flux_fit, flux_fit_err,
            mag_fit, and mag_fit_err.

        """

        # Do photometry on stamp because it will read faster.
        # (We hope.  But CFS latency will kill you at 1 byte.)
        SNLogger.debug( "...make_phot_info_dict reading stamp and psf" )
        diff_img = FITSImageOnDisk( sci_image.diff_stamp_path[ templ_image.image.name ],
                                    noisepath=sci_image.diff_var_stamp_path[ templ_image.image.name ] )
        psf_img = FITSImageOnDisk( sci_image.decorr_psf_path[ templ_image.image.name ] )

        # Required results keys
        req_results_keys = ['mjd', 'flux', 'flux_err', 'zpt', 'NEA', 'sky_rms',
                            'pointing', 'sca', 'pix_x', 'pix_y']
        phrosty_results_keys = ['science_name', 'template_name',
                                'science_id', 'template_id',
                                'template_pointing', 'template_sca',
                                'aperture_sum', 'mag', 'mag_err', 'success']

        results_keys = req_results_keys + phrosty_results_keys
        results_dict = {key: np.nan for key in results_keys}
        results_dict['mjd'] = sci_image.image.mjd
        results_dict['pointing'] = sci_image.image.pointing
        results_dict['sca'] = sci_image.image.sca

        # Additional phrosty keys
        results_dict['science_name'] = sci_image.image.name
        results_dict['template_name'] = templ_image.image.name
        results_dict['science_id'] = str(sci_image.image.id)
        results_dict['template_id'] = str(templ_image.image.id)
        results_dict['template_pointing'] = templ_image.image.pointing
        results_dict['template_sca'] = templ_image.image.sca
        results_dict['success'] = False

        try:
            # Make sure the files are there.  Has side effect of loading the header and data.
            diff_img.get_data( which='data', cache=True )
            diff_img.get_data( which='noise', cache=True )
            # The thing written to disk was actually variance, so fix that
            diff_img.noise = np.sqrt( diff_img.noise )
            psf_img.get_data( which='data', cache=True )
        except Exception:
            # TODO : change this to a more specific exception
            SNLogger.warning( f"Post-processed image files for "
                              f"{self.band}_{sci_image.image.pointing}_{sci_image.image.sca}-"
                              f"{self.band}_{templ_image.image.pointing}_{templ_image.image.sca} "
                              f"do not exist.  Skipping." )
            # results_dict['ap_zpt'] = np.nan

        SNLogger.debug( "...make_phot_info_dict getting psf" )
        coord = SkyCoord(ra=self.diaobj.ra * u.deg, dec=self.diaobj.dec * u.deg)
        pxcoords = skycoord_to_pixel( coord, diff_img.get_wcs().get_astropy_wcs() )
        psf = PSF.get_psf_object( 'OversampledImagePSF',
                                  x=diff_img.data.shape[1]/2., y=diff_img.data.shape[1]/2.,
                                  oversample_factor=1.,
                                  data=psf_img.data )
        SNLogger.debug( "...make_phot_info_dict doing photometry" )
        try:
            results_dict.update( self.phot_at_coords( diff_img, psf, pxcoords=pxcoords, ap_r=ap_r) )
            # Add additional info to the results dictionary so it can be merged into a nice file later.
            SNLogger.debug( "...make_phot_info_dict getting zeropoint" )
            results_dict['zpt'] = sci_image.image.zeropoint
            results_dict['success'] = True

        except:
            # results_dict['ap_zpt'] = np.nan
            SNLogger.debug( f"...make_phot_info_dict failed for \
                             {sci_image.image.name} - {templ_image.image.name}." )

        SNLogger.debug( "...make_phot_info_dict done." )
        return results_dict

    def add_to_results_dict( self, one_pair ):
        """Record results from self.make_phot_info_dict() to the
        aggregate dictionary for the entire light curve.

        Parameters
        ----------
        one_pair : dict
            Dictionary output from self.make_phot_info_dict().
        """

        for key, arr in self.results_dict.items():
            arr.append( one_pair[ key ] )

        if not one_pair['success']:
            self.failures['make_lightcurve'].append({'science': f"{one_pair['filter']} \
                                                                  {one_pair['pointing']} \
                                                                  {one_pair['sca']}",
                                                     'template': f"{one_pair['filter']} \
                                                                   {one_pair['template_pointing']} \
                                                                   {one_pair['template_sca']}"
                                                    })

        SNLogger.debug( "Done adding to results dict" )

    def save_stamp_paths( self, sci_image, templ_image, paths ):
        """Helper function for recording the stamp paths returned in
        self.do_stamps.

        Parameters
        ----------
        sci_image : snappl.image.Image
            Science image with supernova.

        templ_image : snappl.image.Image
            Template image without supernova.

        paths : list of pathlib.Path
            Output from self.do_stamps().

        """
        sci_image.zpt_stamp_path[ templ_image.image.name ] = paths[0]
        sci_image.diff_stamp_path[ templ_image.image.name ] = paths[1]
        sci_image.diff_var_stamp_path[ templ_image.image.name ] = paths[2]

    def do_stamps( self, sci_image, templ_image ):
        """Make stamps from the zero point image, decorrelated
        difference image, and variance image centered at the
        location of the supernova.

        Parameters
        ----------
        sci_image : snappl.image.Image
            Science image with supernova.
        templ_image : snappl.image.Image
            Template image without supernova.

        Returns
        -------
        list of pathlib.Path
            Paths to the stamps corresponding to the zero point image,
            decorrelated difference image, and variance image centered
            at the location of the supernova.
        """

        try:
            zptim = FITSImageOnDisk( sci_image.decorr_zptimg_path[ templ_image.image.name ] )
            zpt_stampname = stampmaker( self.diaobj.ra, self.diaobj.dec, np.array([100, 100]),
                                        zptim,
                                        savedir=self.dia_out_dir,
                                        savename=f"stamp_{zptim.path.name}" )
            zptim.free()

            diffim = FITSImageOnDisk( sci_image.decorr_diff_path[ templ_image.image.name ] )
            diff_stampname = stampmaker( self.diaobj.ra, self.diaobj.dec, np.array([100, 100]),
                                        diffim,
                                        savedir=self.dia_out_dir,
                                        savename=f"stamp_{diffim.path.name}" )
            diffim.free()

            diffvarim = FITSImageOnDisk( sci_image.diff_var_path[ templ_image.image.name ] )
            diffvar_stampname = stampmaker( self.diaobj.ra, self.diaobj.dec, np.array([100, 100]),
                                            diffvarim,
                                            savedir=self.dia_out_dir,
                                            savename=f"stamp_{diffvarim.path.name}" )
            diffvarim.free()

            SNLogger.info(f"Decorrelated diff stamp path: {pathlib.Path( diff_stampname )}")
            SNLogger.info(f"Zpt image stamp path: {pathlib.Path( zpt_stampname )}")
            SNLogger.info(f"Decorrelated diff variance stamp path: {pathlib.Path( diffvar_stampname )}")

            return pathlib.Path( zpt_stampname ), pathlib.Path( diff_stampname ), pathlib.Path( diffvar_stampname )

        except:
            SNLogger.error( f"do_stamps failure for {sci_image.image.pointing} \
                                                    {sci_image.image.sca} - \
                                                    {templ_image.image.pointing} \
                                                    {templ_image.image.sca}: {x} " )
            self.failures['make_stamps'].append({'science': f'{sci_image.image.band} \
                                                              {sci_image.image.pointing} \
                                                              {sci_image.image.sca}',
                                                 'template': f'{templ_image.image.band} \
                                                               {templ_image.image.pointing} \
                                                               {templ_image.image.sca}'
                                                })

    def make_lightcurve( self ):
        """Collect all results from photometry in one dictionary.
        Write the output to a csv as a table.

        Contains its own error logging function, log_error().

        Returns
        -------
        pathlib.Path
            Path to output csv file that contains a light curve.
        """
        SNLogger.info( "Making lightcurve." )

        self.metadata = {
            'provenance_id': None,
            'diaobject_id': self.diaobj.id,
            'diaobject_position_id': None, 
            'iau_name': self.diaobj.iauname,
            'band': self.band,
            'ra': self.diaobj.ra,
            'dec': self.diaobj.dec,
            'ra_err': None,
            'dec_err': None,
            'ra_dec_covar': None,
            f'local_surface_brightness_{self.band}': -999.  # phrosty does not output this value!
                                                            # add it later!
        }

        self.results_dict = {
            # Required keys in specified order
            'mjd': [],  # Days, float
            'flux': [],  # DN/s, float
            'flux_err': [],  # DN/s, float
            'zpt': [],  # AB mag of object is m = -2.5 * log(flux) + zpt, float
            'NEA': [],  # px^2, float
            'sky_rms': [],  # DN/s, float
            'pointing': [],  # int/string, temporary name
            'sca': [],  # int
            'pix_x': [],  # x-position of SN on detector w/ 0-offset, float
            'pix_y': [],  # y-position of SN on detector w/ 0-offset, float

            # Additional phrosty keys
            'science_name': [],
            'template_name': [],
            'science_id': [],
            'template_id': [],
            'template_pointing': [],
            'template_sca': [],
            'aperture_sum': [],
            'mag': [],
            'mag_err': [],
            'success': []
        }

        def log_error( sci_image, templ_image, x ):
            SNLogger.error( f"make_phot_info_dict failure for {sci_image.image.pointing} \
                                                              {sci_image.image.sca} - \
                                                              {templ_image.image.pointing} \
                                                              {templ_image.image.sca}: {x}" )
            self.failures['make_lightcurve'].append({'science': f'{sci_image.image.band} \
                                                                  {sci_image.image.pointing} \
                                                                  {sci_image.image.sca}',
                                                     'template': f'{templ_image.image.band} \
                                                                   {templ_image.image.pointing} \
                                                                   {templ_image.image.sca}'
                                                    })

        if self.nprocs > 1:
            with Pool( self.nprocs ) as pool:
                for sci_image in self.science_images:
                    for templ_image in self.template_images:
                        logerr_partial = partial(log_error, sci_image, templ_image)
                        pool.apply_async( self.make_phot_info_dict, (sci_image, templ_image), {},
                                          self.add_to_results_dict,
                                          error_callback=logerr_partial )
                pool.close()
                pool.join()
        else:
            for i, sci_image in enumerate( self.science_images ):
                SNLogger.debug( f"Doing science image {i} of {len(self.science_images)}" )
                for templ_image in self.template_images:
                    self.add_to_results_dict( self.make_phot_info_dict( sci_image, templ_image ) )

        imgprov = Provenance.get_by_id( self.science_images[0].image.provenance_id, dbclient=self.dbclient )
        objprov = Provenance.get_by_id( self.diaobj.provenance_id, dbclient=self.dbclient )
        
        phrosty_version = phrosty.__version__
        major = int(phrosty_version.split('.')[0])
        minor = int(phrosty_version.split('.')[1])

        ltcvprov = Provenance( process='phrosty',
                                major=major, 
                                minor=minor,
                                params=Config.get(),
                                keepkeys=[ 'photometry.phrosty' ],
                                omitkeys=None,
                                upstreams=[imgprov, objprov],
                                )

        self.metadata['provenance_id'] = ltcvprov.id
        lc_obj = Lightcurve(data=self.results_dict, meta=self.metadata)
        lc_obj.diaobj = self.diaobj
        lc_obj.provenance_object = ltcvprov

        if self.dbsave:
            SNLogger.debug( "Saving results to database..." )
            ltcvprov.save_to_db( tag=self.ltcv_prov_tag )
            lc_obj.write()
            lc_obj.save_to_db( dbclient=self.dbclient )

        else:
            SNLogger.debug( "Saving results using paths..." )
            if self.diaobj.id is not None:
                save_basename = str(self.diaobj.id)
            else:
                save_basename = str(self.oid)

            filepath = pathlib.Path(f'data/{self.oid}/{save_basename}.pq')
            results_savepath = f'{self.ltcv_dir}/{filepath}'
            lc_obj.write( base_dir=self.ltcv_dir, filepath=filepath, overwrite=True)

            SNLogger.info(f'Results saved to {results_savepath}.')

            return results_savepath

    def write_fits_file( self, data, header, savepath ):
        """Helper function for writing fits files.

        Parameters
        ----------
        data : np.array
            Image array.
        header : _type_
            FITS header.
        savepath : str
            Savepath for FITS file.
        """

        try:
            if header is not None:
                hdr_dict = dict(header.items())
            else:
                hdr_dict = None
            fitsio.write( savepath, data, header=hdr_dict, clobber=True )
        except Exception as e:
            SNLogger.exception( f"Exception writing FITS image {savepath}: {e}" )
            raise

    def clear_contents( self, directory ):
        """Delete contents of a directory. Used to clear temporary
        files.

        Parameters
        ----------
        directory : pathlib.Path
            Path to directory to empty.
        """
        for f in directory.iterdir():
            try:
                if f.is_dir():
                    shutil.rmtree( f )
                else:
                    f.unlink()

            except Exception as e:
                print( f'Oops! Deleting {f} from {directory} did not work.\nReason: {e}' )

    def __call__( self, through_step=None ):
        """Run the pipeline.

        Parameters
        ----------
          through_step: str, default None
             Which step to run thorough?  Runs them all if not given.

             Steps in order are:
             * sky_subtract
             * get_psfs
             * align_and_preconvolve
             * subtract
             * find_decorrelation
             * apply_decorrelation
             * make_stamps
             * make_lightcurve

        Returns
        -------
          ltcvpath : pathlib.Path or None
            The path to the output lightcurve file if make_lightcurve
            was run, otherwise None.

        """
        if self.mem_trace:
            tracemalloc.start()
            tracemalloc.reset_peak()

        if through_step is None:
            through_step = 'make_lightcurve'

        steps = [ 'sky_subtract', 'get_psfs', 'align_and_preconvolve', 'subtract', 'find_decorrelation',
                  'apply_decorrelation', 'make_stamps', 'make_lightcurve' ]
        stepdex = steps.index( through_step )
        if stepdex < 0:
            raise ValueError( f"Unknown step {through_step}" )
        steps = steps[:stepdex+1]

        if 'sky_subtract' in steps:
            # After this step is done, all images (both science and template)
            #   will have the following fields set:
            #     .skysub_img : sub-subtracted image
            #     .detmask_img : deteciton mask image
            #     .skyrms : float, median of sky image calculated by SExtractor
            SNLogger.info( "Running sky subtraction" )
            with nvtx.annotate( "skysub", color=0xff8888 ):
                self.sky_sub_all_images()

        if self.mem_trace:
            SNLogger.info( f"After sky_subtract, memory usage = {tracemalloc.get_traced_memory()[1]/(1024**2):.2f} MB" )

        if 'get_psfs' in steps:
            # After this step, all images (both science and template) will
            #   have the .psf_data field set with the image-resolution PSF stamp data
            #   (from a PSF.get_stamp() call.)
            SNLogger.info( "Getting PSFs" )
            with nvtx.annotate( "getpsfs", color=0xff8888 ):
                self.get_psfs()

        if self.mem_trace:
            SNLogger.info( f"After get_psfs, memory usage = {tracemalloc.get_traced_memory()[1]/(1024**2):.2f} MB" )

        # Create a process pool to write fits files
        with Pool( self.nwrite ) as fits_writer_pool:

            def log_fits_write_error( savepath, x ):
                SNLogger.error( f"Exception writing FITS file {savepath}: {x}" )
                # raise?

            # Do the hardcore processing

            for templ_image in self.template_images:
                for sci_image in self.science_images:
                    SNLogger.info( f"Processing {sci_image.image.name} minus {templ_image.image.name}" )
                    sfftifier = None
                    fail_info = {'science': f'{sci_image.image.band} \
                                              {sci_image.image.pointing} \
                                              {sci_image.image.sca}',
                                 'template': f'{templ_image.image.band} \
                                               {templ_image.image.pointing} \
                                               {templ_image.image.sca}'
                                }
                    i_failed = False

                    if 'align_and_preconvolve' in steps:
                        # After this step, sfftifier will be a SpaceSFFT_CupyFlow
                        #  object with the following fields filled.  All cupy arrays
                        #  are float64 (even the DMASK!), and all of them are Transposed
                        #  from what's stored in the PipelineImage object
                        #    hdr_target : FITS header of science image
                        #    hdr_object : FITS header of template image
                        #    target_skyrms : float, median of sky for science image
                        #    object_skyrms : float, median of sky for template image
                        #    PixA_target_GPU  : science image data on GPU
                        #    PixA_Ctarget_GPU : cross-convolved [with what?] science image on GPU
                        #    PixA_object_GPU  : template image data on GPU
                        #    PixA_resamp_object_GPU : warped template image data on GPU
                        #    PixA_Cresampl_object_GPU : cross-convolved template image on GPU
                        #    BlankMask_GPU : boolean array, True where PixA_resampl_object_GPU is 0.
                        #    PixA_targetVar_GPU : noise CHECK THIS image for science image, on GPU
                        #    PixA_objectVar_GPU : noise CHECK THIS image for template image, on GPU
                        #    PiXA_resampl_objectVar_GPU : warped template variance data on GPU
                        #    PixA_target_DMASK_GPU : detmask for science image, on GPU
                        #    PixA_object_DMASK_GPU : detmask for template image, on GPU
                        #    PixA_resampl_object_DMASK_GPU : warped dmask for template image, on GPU
                        #    PSF_target_GPU : PSF stamp for science image
                        #    PSF_object_GPU : PSF stamp for template image
                        #    PSF_resampl_object_GPU : PSF stamp for template image rotated & resampled, on GPU
                        #    sci_is_target : True
                        #    GKerHW : 9 (int half-width of matching kernel, full width is 2*GKerHW+1)
                        #    KerPolyOrder : config value of photometry.phrosty.kerpolyorder
                        #    BGPolyOrder : 0
                        #    ConstPhotRatio : True
                        #    CUDA_DEVICE_4SUBTRACT : '0'
                        #    GAIN : 1.0
                        #    RANDOM_SEED : 10086
                        SNLogger.info( "...align_and_preconvolve" )
                        with nvtx.annotate( "align_and_pre_convolve", color=0x8888ff ):
                            try:
                                sfftifier = self.align_and_pre_convolve( templ_image, sci_image )
                            except:
                                i_failed = True
                                self.failures['align_and_preconvolve'].append(fail_info)

                    if 'subtract' in steps:
                        # After this step is done, two more fields in sfftifier are set:
                        #    Solution_GPU  : --something--??
                        #    PixA_DIFF_GPU : difference image ??
                        # Get Lei to write docs on PureCupy_Customized_Packet.PCCP so
                        #   we can figure out what these are
                        SNLogger.info( "...subtract" )
                        with nvtx.annotate( "subtraction", color=0x44ccff ):
                            try:
                                sfftifier.sfft_subtraction()
                            except:
                                i_failed = True
                                self.failures['subtract'].append(fail_info)

                    if 'find_decorrelation' in steps:
                        # This step does ...
                        # After it's done, the following fields of sfftifier are set:
                        #   Solution : CPU copy of Solution_GPU
                        #   FKDECO_GPU : result of PureCupy_Decorrelation_Calculator.PCDC
                        # Get Lei to write documentation on P:ureCupy_DeCorrelation_Calculator
                        #   so we can figure out what this is, but I THINK it's a
                        #   kernel that is used to convolve with things to "decorrelate".
                        #   (From what?)
                        # In addition the two local varaibles diff_var and diff_var_path are set.
                        #   diff_var : variance in difference image (I THINK), on GPU
                        #   diff_var_path : where we want to write diff_var in self.dia_out_dir
                        SNLogger.info( "...find_decorrelation" )
                        with nvtx.annotate( "find_decor", color=0xcc44ff ):
                            try:
                                sfftifier.find_decorrelation()
                            except:
                                i_failed = True
                                self.failures['find_decorrelation'].append(fail_info)

                        SNLogger.info( "...generate variance image" )
                        with nvtx.annotate( "variance", color=0x44ccff ):
                            try:
                                diff_var = sfftifier.create_variance_image()
                                mess = f"{sci_image.image.name}-{templ_image.image.name}"
                                diff_var_path = self.dia_out_dir / f"diff_var_{mess}"
                            except:
                                i_failed = True
                                self.failures['variance'].append(fail_info)

                    if 'apply_decorrelation' in steps and not i_failed:
                        try:
                            mess = f"{sci_image.image.name}-{templ_image.image.name}"
                            decorr_psf_path = self.dia_out_dir / f"decorr_psf_{mess}"
                            decorr_zptimg_path = self.dia_out_dir / f"decorr_zptimg_{mess}"
                            decorr_diff_path = self.dia_out_dir / f"decorr_diff_{mess}"

                            images =    [ sfftifier.PixA_DIFF_GPU,    diff_var,
                                        sfftifier.PixA_Ctarget_GPU, sfftifier.PSF_target_GPU ]
                            savepaths = [ decorr_diff_path,           diff_var_path,
                                        decorr_zptimg_path,         decorr_psf_path ]
                            headers =   [ sfftifier.hdr_target,       sfftifier.hdr_target,
                                        sfftifier.hdr_target,       None ]

                            for img, savepath, hdr in zip( images, savepaths, headers ):
                                with nvtx.annotate( "apply_decor", color=0xccccff ):
                                    SNLogger.info( f"...apply_decor to {savepath}" )
                                    decorimg = sfftifier.apply_decorrelation( img )
                                with nvtx.annotate( "submit writefits", color=0xff8888 ):
                                    SNLogger.info( f"...writefits {savepath}" )
                                    fits_writer_pool.apply_async( self.write_fits_file,
                                                                ( cp.asnumpy( decorimg ).T, hdr, savepath ), {},
                                                                error_callback=partial(log_fits_write_error, savepath) )
                            sci_image.decorr_psf_path[ templ_image.image.name ] = decorr_psf_path
                            sci_image.decorr_zptimg_path[ templ_image.image.name ] = decorr_zptimg_path
                            sci_image.decorr_diff_path[ templ_image.image.name ] = decorr_diff_path
                            sci_image.diff_var_path[ templ_image.image.name ] = diff_var_path
                        except:
                            i_failed = True
                            self.failures['apply_decorrelation'].append(fail_info)

                    if self.keep_intermediate and not i_failed:
                        # Each key is the file prefix addition.
                        # Each list has [descriptive filetype, image file name, data, header].

                        # TODO: Include multiprocessing.
                        # In the future, we may want to write these things right after they happen
                        # instead of saving it all for the end of the SFFT stuff.

                        write_filepaths = {'aligned': [['img',
                                                        f'{templ_image.image.name}_-_{sci_image.image.name}.fits',
                                                        cp.asnumpy(sfftifier.PixA_resamp_object_GPU.T),
                                                        sfftifier.hdr_target],
                                                        ['var',
                                                         f'{templ_image.image.name}_-_{sci_image.image.name}.fits',
                                                         cp.asnumpy(sfftifier.PixA_resamp_objectVar_GPU.T),
                                                         sfftifier.hdr_target],
                                                       ['psf',
                                                        f'{templ_image.image.name}_-_{sci_image.image.name}.fits',
                                                        cp.asnumpy(sfftifier.PSF_resamp_object_GPU.T),
                                                        sfftifier.hdr_target],
                                                       ['detmask',
                                                        f'{sci_image.image.name}_-_{templ_image.image.name}.fits',
                                                        cp.asnumpy(sfftifier.PixA_resamp_object_DMASK_GPU.T),
                                                        sfftifier.hdr_target]
                                                       ],
                                           'convolved': [['img',
                                                         f'{sci_image.image.name}_-_{templ_image.image.name}.fits',
                                                         cp.asnumpy(sfftifier.PixA_Ctarget_GPU.T),
                                                         sfftifier.hdr_target],
                                                         ['img',
                                                         f'{templ_image.image.name}_-_{sci_image.image.name}.fits',
                                                         cp.asnumpy(sfftifier.PixA_Cresamp_object_GPU.T),
                                                         sfftifier.hdr_target]
                                                        ],
                                           'diff':     [['img',
                                                        f'{sci_image.image.name}_-_{templ_image.image.name}.fits',
                                                        cp.asnumpy(sfftifier.PixA_DIFF_GPU.T),
                                                        sfftifier.hdr_target]
                                                       ],
                                           'decorr':   [['kernel',
                                                        f'{sci_image.image.name}_-_{templ_image.image.name}.fits',
                                                        cp.asnumpy(sfftifier.FKDECO_GPU.T),
                                                        sfftifier.hdr_target]
                                                       ]
                                          }
                        # Write the aligned images
                        for key in write_filepaths.keys():
                            for (imgtype, name, data, header) in write_filepaths[key]:
                                savepath = self.scratch_dir / f'{key}_{imgtype}_{name}'
                                self.write_fits_file( data, header, savepath=savepath )

                    SNLogger.info( f"DONE processing {sci_image.image.name} minus {templ_image.image.name}" )
                    if self.mem_trace:
                        SNLogger.info( f"After preprocessing, subtracting, and postprocessing \
                                        a science image, memory usage = \
                                         {tracemalloc.get_traced_memory()[1]/(1024**2):.2f} MB" )

                    sci_image.free()

                SNLogger.info( f"DONE with all science images for template {templ_image.image.name}" )
                templ_image.free()

            SNLogger.info( "Waiting for FITS writer processes to finish" )
            with nvtx.annotate( "fits_write_wait", color=0xff8888 ):
                fits_writer_pool.close()
                fits_writer_pool.join()
            SNLogger.info( "...FITS writer processes done." )

        if 'make_stamps' in steps and not i_failed:
            SNLogger.info( "Starting to make stamps..." )
            with nvtx.annotate( "make stamps", color=0xff8888 ):

                def log_stamp_err( sci_image, templ_image, x ):
                    SNLogger.error( f"do_stamps failure for {sci_image.image.pointing} \
                                                            {sci_image.image.sca} - \
                                                            {templ_image.image.pointing} \
                                                            {templ_image.image.sca}: {x} " )
                    self.failures['make_stamps'].append({'science': f'{sci_image.image.band} \
                                                                      {sci_image.image.pointing} \
                                                                      {sci_image.image.sca}',
                                                         'template': f'{templ_image.image.band} \
                                                                       {templ_image.image.pointing} \
                                                                       {templ_image.image.sca}'
                                                        })

                partialstamp = partial(stampmaker, self.diaobj.ra, self.diaobj.dec, np.array([100, 100]))
                # template path, savedir, savename
                templstamp_args = ( (ti.image, self.dia_out_dir, f'stamp_{str(ti.image.name)}')
                                    for ti in self.template_images )
                if self.nwrite > 1:
                    with Pool( self.nwrite ) as templ_stamp_pool:
                        templ_stamp_pool.starmap_async( partialstamp, templstamp_args,
                                                        error_callback=log_stamp_err )
                        templ_stamp_pool.close()
                        templ_stamp_pool.join()

                    with Pool( self.nwrite ) as sci_stamp_pool:
                        for sci_image in self.science_images:
                            for templ_image in self.template_images:
                                pair = (sci_image, templ_image)
                                stamperr_partial = partial(log_stamp_err, sci_image, templ_image)
                                sci_stamp_pool.apply_async( self.do_stamps, pair, {},
                                                            callback = partial(self.save_stamp_paths,
                                                                               sci_image, templ_image),
                                                            error_callback=stamperr_partial )
                        sci_stamp_pool.close()
                        sci_stamp_pool.join()

                else:
                    for tsargs in templstamp_args:
                        partialstamp(*tsargs)

                    for sci_image in self.science_images:
                        for templ_image in self.template_images:
                            stamp_paths = self.do_stamps( sci_image, templ_image)
                            self.save_stamp_paths( sci_image, templ_image, stamp_paths )

            SNLogger.info('...finished making stamps.')

        if self.mem_trace:
            SNLogger.info( f"After make_stamps, memory usage = {tracemalloc.get_traced_memory()[1]/(1024**2):.2f} MB" )

        lightcurve_path = None
        if 'make_lightcurve' in steps and not i_failed:
            with nvtx.annotate( "make_lightcurve", color=0xff8888 ):
                lightcurve_path = self.make_lightcurve()

        if self.mem_trace:
            SNLogger.info( f"After make_lightcurve, memory usage = \
                            {tracemalloc.get_traced_memory()[1]/(1024**2):.2f} MB" )

        if self.remove_temp_dir:
            self.clear_contents( self.temp_dir )

        if np.sum( [len(x) for x in self.failures.values()] ) > 0:
            SNLogger.info( f"There were some failures here! They were:\n{self.failures}" )

        else:
            SNLogger.info( f"No failures here! Fail list:\n{self.failures}" )

        if lightcurve_path is None:
            SNLogger.info( f"Light curves saved to database!" )

        return lightcurve_path


# ======================================================================


def main():
    # Run one arg pass just to get the config file, so we can augment
    #   the full arg parser later with config options
    configparser = argparse.ArgumentParser( add_help=False )
    configparser.add_argument( '-c', '--config-file', default=None,
                               help=( "Location of the .yaml config file; defaults to the value of the "
                                      "SNPIT_CONFIG environment varaible." ) )
    args, leftovers = configparser.parse_known_args()

    try:
        cfg = Config.get( args.config_file, setdefault=True )
    except RuntimeError as e:
        if str(e) == 'No default config defined yet; run Config.init(configfile)':
            sys.stderr.write( "Error, no configuration file defined.\n"
                              "Either run phrosty with -c <configfile>\n"
                              "or set the SNPIT_CONFIG environment variable.\n" )
            sys.exit(1)
        else:
            raise

    parser = argparse.ArgumentParser()
    # Put in the config_file argument, even though it will never be found, so it shows up in help
    parser.add_argument( '-c', '--config-file', help="Location of the .yaml config file" )

    # Running options
    parser.add_argument( '-p', '--nprocs', type=int, default=1,
                         help="Number of process for multiprocessing steps (e.g. skysub)" )
    parser.add_argument( '-w', '--nwrite', type=int, default=5,
                         help="Number of parallel FITS writing processes" )
    parser.add_argument( '-v', '--verbose', action='store_true', default=False,
                         help="Show debug log info" )
    parser.add_argument( '--through-step', default='make_lightcurve',
                         help="Stop after this step; one of (see above)" )
    parser.add_argument( '--dbsave', action='store_true',
                         help="Toggle saving to the database." )
    parser.add_argument( '--memtrace', action='store_true',
                         help="Toggle memory tracing with tracemalloc.")

    # Object collections
    parser.add_argument( '-oc', '--object-collection', default='snpitdb',
                         help='Collection of the object.  Currently, "snpitdb", "ou2024", and "manual" supported.' )
    parser.add_argument( '-os', '--object-subset', default=None,
                         help="Collection subset.  Not used by all collections." )

    # SN and observation information
    parser.add_argument( '--oid', type=int, required=True,
                         help="Object ID.  Meaning is collection-dependent." )
    parser.add_argument( '-r', '--ra', type=float, default=None,
                         help="Object RA.  By default, uses the one found for the object." )
    parser.add_argument( '-d', '--dec', type=float, default=None,
                         help="Object Dec.  By default, uses the one found for the object." )
    parser.add_argument( '-b', '--band', type=str, required=True,
                         help="Band: R062, Z087, Y106, J129, H158, F184, or K213" )

    # Required args for using SN PIT database
    # DiaObj
    parser.add_argument( '-did', '--diaobject-id', type=str, default=None,
                         help="ID for DiaObject. Required to use SN PIT database. \
                               Invalid if --image-collection is not snpitdb." )
    parser.add_argument( '-dpt', '--diaobject-provenance-tag', type=str, default=None,
                         help="Provenance tag for DiaObject. Required to use SN PIT database. \
                               Invalid if --image-collection is not snpitdb." )
    parser.add_argument( '-dp', '--diaobject-process', type=str, required=False, default=None,
                         help="Process for DiaObject. Required to use SN PIT database. \
                               Invalid if --image-collection is not snpitdb.")
    parser.add_argument( '-dppt', '--diaobject-position-provenance-tag', type=str, default=None,
                         help="Provenance tag for the position of the DiaObject." )
    parser.add_argument( '-dpp', '--diaobject-position-process', type=str, default=None,
                         help="The process where the DiaObject position originated." )

    # Lightcurve
    parser.add_argument( '-lpi', '--ltcv-provenance-id', type=str, default=None,
                         help="Provenance ID for lightcurve. Required to use SN PIT database. \
                               Invalid if --image-collection is not snpitdb." )
    parser.add_argument( '-lpt', '--ltcv-provenance-tag', type=str, default=None,
                         help="Provenance tag for lightcurve. Required to use SN PIT database. \
                               Invalid if --image-collection is not snpitdb." )
    parser.add_argument( '-lp', '--ltcv-process', type=str, default='phrosty',
                         help="Process for light curve. Required to use SN PIT database. \
                               Invalid if --image-collection is not snpitdb." )

    # Image collection
    parser.add_argument( '-ic', '--image-collection', default='snpitdb',
                         help="Collection of the images we're using. For SN PIT database, use snpitdb. \
                               Currently supported: ou2024, manual_fits, snpitdb (default)." )
    parser.add_argument( '-is', '--image-subset', default=None,
                         help="Image collection subset. To use SN PIT database, must be None." )

    # Image
    parser.add_argument( '-ipt', '--image-provenance-tag', default=None,
                         help='Provenance tag for images. Required to use SN PIT database. \
                               Invalid if --image-collection is not snpitdb.' )
    parser.add_argument( '-ip', '--image-process', default=None,
                         help='Image process. Required to use SN PIT database. \
                               Invalid if --image-collection is not snpitdb.' )

    # Path-based options
    parser.add_argument( '--base-path', type=str, default=None,
                         help='Base path for reading images. Required for "manual_fits" image collection.' )
    parser.add_argument( '-t', '--template-images', type=str, default=None,
                         help="Path to file with, per line, ( path_to_image, pointing, sca, mjd, band )" )
    parser.add_argument( '-s', '--science-images', type=str, default=None,
                         help="Path to file with, per line, ( path_to_image, pointing, sca, mjd, band )" )

    cfg.augment_argparse( parser )
    args = parser.parse_args( leftovers )
    cfg.parse_args( args )

    if args.base_path is None and args.image_collection == 'manual_fits':
        SNLogger.error( 'Must provide --base-path if --image-collection is manual_fits.' )
        raise ValueError( f'args.base_path is {args.base_path}.' )

    if args.image_collection == 'snpitdb' and args.image_provenance_tag is None:
        SNLogger.error( 'Must provide --image-provenance-tag if --image-collection is snpitdb.' )
        raise ValueError( f'args.image_provenance_tag is {args.image_provenance_tag}.' )

    dbclient = SNPITDBClient()
    # Get the DiaObject, update the RA and Dec
    if args.diaobject_id is None:
        diaobjs = DiaObject.find_objects( collection=args.object_collection, 
                                        #   subset=args.object_subset,
                                          provenance_tag=args.diaobject_provenance_tag,
                                          process=args.diaobject_process,
                                          name=args.oid, ra=args.ra, dec=args.dec )

    if len( diaobjs ) == 0:
        raise ValueError( f"Could not find DiaObject with id={args.id}, ra={args.ra}, dec={args.dec}." )
    if len( diaobjs ) > 1:
        raise ValueError( f"Found multiple DiaObject with id={args.id}, ra={args.ra}, dec={args.dec}." )
    diaobj = diaobjs[0]
    if args.ra is not None:
        if np.fabs( args.ra - diaobj.ra ) > 1. / 3600. / np.cos( diaobj.dec * np.pi / 180. ):
            SNLogger.warning( f"Given RA {args.ra} is far from DiaObject nominal RA {diaobj.ra}" )
        diaobj.ra = args.ra
    if args.dec is not None:
        if np.fabs( args.dec - diaobj.dec ) > 1. / 3600.:
            SNLogger.warning( f"Given Dec {args.dec} is far from DiaObject nominal Dec {diaobj.dec}" )
        diaobj.dec = args.dec

    # Get the image collection
    imgcol = ImageCollection.get_collection( collection=args.image_collection,
                                             subset=args.image_subset,
                                             provenance_tag=args.image_provenance_tag,
                                             process=args.image_process,
                                             base_path=args.base_path,
                                             dbclient=dbclient
                                           )

    if args.image_collection == 'snpitdb':
        found_images = imgcol.find_images( ra=diaobj.ra,
                                           dec=diaobj.dec,
                                           band=args.band,
                                           dbclient=dbclient
                                         )
    else:
        found_images = imgcol.find_images( ra=diaobj.ra,
                                           dec=diaobj.dec,
                                           band=args.band
                                         )

    fetched_prov = Provenance.get_provs_for_tag( tag=args.diaobject_provenance_tag,
                                                 process=args.diaobject_process
                                               )

    # Create and launch the pipeline
    pipeline = Pipeline( diaobj, imgcol, args.band,
                         science_csv=args.science_images,
                         template_csv=args.template_images,
                         oid=args.oid,
                         ltcv_prov_tag=args.ltcv_provenance_tag,
                         dbsave=args.dbsave,
                         dbclient=dbclient,
                         nprocs=args.nprocs,
                         nwrite=args.nwrite,
                         verbose=args.verbose,
                         memtrace=args.memtrace )

    pipeline( args.through_step )


# ======================================================================
if __name__ == "__main__":
    main()
