__all__ = [ 'PipelineImage', 'Pipeline' ]

# Imports STANDARD
import sys
import argparse
import cupy as cp
from functools import partial
import logging
from multiprocessing import Pool
import numpy as np
import nvtx
import pathlib
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
from phrosty.imagesubtraction import sky_subtract, stampmaker
from sfft.SpaceSFFTCupyFlow import SpaceSFFT_CupyFlow
from snappl.diaobject import DiaObject
from snappl.imagecollection import ImageCollection
from snappl.image import FITSImageOnDisk
from snappl.psf import PSF
from snpit_utils.config import Config
from snpit_utils.logger import SNLogger


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
            self.save_dir = pathlib.Path( self.config.value( 'photometry.phrosty.paths.scratch_dir' ) )
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
        try:
            return sky_subtract( self.image )
        except Exception as ex:
            SNLogger.exception( ex )
            raise

    def save_sky_subtract_info( self, info ):
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
                  nprocs=1, nwrite=5,
                  verbose=False ):

        """Create the a pipeline object.

        Parameters
        ----------
           diaobj : DiaObject
             The object we're building a lightcurve for

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

           nprocs: int, default 1
             Number of cpus for the CPU multiprocessing segments of the pipeline.
             (GPU segments will run a single process.)

           nwrite: int, default 5
             Number of asynchronous FITS writer processes.

        """

        SNLogger.setLevel( logging.DEBUG if verbose else logging.INFO )
        self.config = Config.get()
        self.imgcol = imgcol
        self.diaobj = diaobj
        self.band = band

        self.dia_out_dir = pathlib.Path( self.config.value( 'photometry.phrosty.paths.dia_out_dir' ) )
        self.scratch_dir = pathlib.Path( self.config.value( 'photometry.phrosty.paths.scratch_dir' ) )
        self.temp_dir_parent = pathlib.Path( self.config.value( 'photometry.phrosty.paths.temp_dir' ) )
        self.temp_dir = self.temp_dir_parent / str(uuid.uuid1())
        self.temp_dir.mkdir()
        self.ltcv_dir = pathlib.Path( self.config.value( 'photometry.phrosty.paths.ltcv_dir' ) )

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

        self.nprocs = nprocs
        self.nwrite = nwrite

        self.keep_intermediate = self.config.value( 'photometry.phrosty.keep_intermediate' )
        self.remove_temp_dir = self.config.value( 'photometry.phrosty.remove_temp_dir' )
        self.mem_trace = self.config.value( 'photometry.phrosty.mem_trace' )


    def _read_csv( self, csvfile ):
        imlist = []
        with open( csvfile ) as ifp:
            hdrline = ifp.readline()
            if not re.search( r"^\s*path\s+pointing\s+sca\s+mjd\s+band\s*$", hdrline ):
                raise ValueError( f"First line of list file {csvfile} didn't match what was expected." )
            for line in ifp:
                path, pointing, sca, mjd, band = line.split()
                if band == self.band:
                    # This should yell at us if the pointing or sca doesn't match what is read from the path
                    imlist.append( self.imgcol.get_image( path=path, pointing=pointing, sca=sca, band=band ) )
        return imlist


    def sky_sub_all_images( self ):
        # Currently, this writes out a bunch of FITS files.  Further refactoring needed
        #   to support more general image types.
        all_imgs = self.science_images.copy()     # shallow copy
        all_imgs.extend( self.template_images )

        self._omg = False

        def log_error( img, x ):
            SNLogger.error( f"Sky subtraction subprocess failure: {x} for image {img.image.path}" )
            self._omg = True

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

        if self._omg:
            raise RuntimeError( "Sky subtraction errors." )


    def get_psfs( self ):
        all_imgs = self.science_images.copy()     # shallow copy
        all_imgs.extend( self.template_images )

        self._omg = False

        def log_error( x ):
            SNLogger.error( f"get_psf subprocess failure: {x}" )
            self._omg = True    # noqa: F841

        if self.nprocs > 1:
            with Pool( self.nprocs ) as pool:
                for img in all_imgs:
                    # callback_partial = partial( img.save_psf_path, all_imgs )
                    pool.apply_async( img.get_psf, (self.diaobj.ra, self.diaobj.dec), {},
                                      img.keep_psf_data, log_error )
                pool.close()
                pool.join()
        else:
            for img in all_imgs:
                img.keep_psf_data( img.get_psf(self.diaobj.ra, self.diaobj.dec) )

        if self._omg:
            raise RuntimeError( "get_psf errors." )


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
                        'aperture_sum': init['flux_init'][0],
                        'flux_fit': flux,
                        'flux_fit_err': flux_err,
                        'mag_fit': mag,
                        'mag_fit_err': mag_err
                        }

        return results_dict

    def make_phot_info_dict( self, sci_image, templ_image, ap_r=4 ):
        """"Do things.

        Parmaeters
        ----------
          sci_image: PipelineImage
            science image wrapper

          temp_image: PipelineImage
            template image wrapper

          ap_r: float, default 4
             Radius of aperture to use in something

        Returns
        -------
          something: dict
            It has things in it

        """

        # Do photometry on stamp because it will read faster.
        # (We hope.  But CFS latency will kill you at 1 byte.)
        SNLogger.debug( "...make_phot_info_dict reading stamp and psf" )
        diff_img = FITSImageOnDisk( sci_image.diff_stamp_path[ templ_image.image.name ],
                                    noisepath=sci_image.diff_var_stamp_path[ templ_image.image.name ] )
        psf_img = FITSImageOnDisk( sci_image.decorr_psf_path[ templ_image.image.name ] )

        results_dict = {}
        results_dict['sci_name'] = sci_image.image.name
        results_dict['templ_name'] = templ_image.image.name
        results_dict['success'] = False
        results_dict['ra'] = self.diaobj.ra
        results_dict['dec'] = self.diaobj.dec
        results_dict['mjd'] = sci_image.image.mjd
        results_dict['filter'] = self.band
        results_dict['pointing'] = sci_image.image.pointing
        results_dict['sca'] = sci_image.image.sca
        results_dict['template_pointing'] = templ_image.image.pointing
        results_dict['template_sca'] = templ_image.image.sca

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
            results_dict['zpt'] = np.nan
            results_dict['ap_zpt'] = np.nan
            results_dict['aperture_sum'] = np.nan
            results_dict['flux_fit'] = np.nan
            results_dict['flux_fit_err'] = np.nan
            results_dict['mag_fit'] = np.nan
            results_dict['mag_fit_err'] = np.nan

        SNLogger.debug( "...make_phot_info_dict getting psf" )
        coord = SkyCoord(ra=self.diaobj.ra * u.deg, dec=self.diaobj.dec * u.deg)
        pxcoords = skycoord_to_pixel( coord, diff_img.get_wcs().get_astropy_wcs() )
        psf = PSF.get_psf_object( 'OversampledImagePSF',
                                  x=diff_img.data.shape[1]/2., y=diff_img.data.shape[1]/2.,
                                  oversample_factor=1.,
                                  data=psf_img.data )
        SNLogger.debug( "...make_phot_info_dict doing photometry" )
        results_dict.update( self.phot_at_coords( diff_img, psf, pxcoords=pxcoords, ap_r=ap_r) )

        # Add additional info to the results dictionary so it can be merged into a nice file later.
        SNLogger.debug( "...make_phot_info_dict getting zeropoint" )
        results_dict['zpt'] = sci_image.image.zeropoint
        results_dict['success'] = True

        SNLogger.debug( "...make_phot_info_dict done." )
        return results_dict

    def add_to_results_dict( self, one_pair ):
        for key, arr in self.results_dict.items():
            arr.append( one_pair[ key ] )
        SNLogger.debug( "Done adding to results dict" )

    def save_stamp_paths( self, sci_image, templ_image, paths ):
        sci_image.zpt_stamp_path[ templ_image.image.name ] = paths[0]
        sci_image.diff_stamp_path[ templ_image.image.name ] = paths[1]
        sci_image.diff_var_stamp_path[ templ_image.image.name ] = paths[2]

    def do_stamps( self, sci_image, templ_image ):

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

    def make_lightcurve( self ):
        SNLogger.info( "Making lightcurve." )

        self.results_dict = {
            'ra': [],
            'dec': [],
            'mjd': [],
            'filter': [],
            'pointing': [],
            'sca': [],
            'template_pointing': [],
            'template_sca': [],
            'zpt': [],
            'aperture_sum': [],
            'flux_fit': [],
            'flux_fit_err': [],
            'mag_fit': [],
            'mag_fit_err': [],
        }

        self._omg = False

        def log_error( x ):
            SNLogger.error( f"make_phot_info_dict subprocess failure: {x}" )
            self._omg = True

        if self.nprocs > 1:
            with Pool( self.nprocs ) as pool:
                for sci_image in self.science_images:
                    for templ_image in self.template_images:
                        pool.apply_async( self.make_phot_info_dict, (sci_image, templ_image), {},
                                          self.add_to_results_dict,
                                          error_callback=log_error )
                pool.close()
                pool.join()
        else:
            for i, sci_image in enumerate( self.science_images ):
                SNLogger.debug( f"Doing science image {i} of {len(self.science_images)}" )
                for templ_image in self.template_images:
                    self.add_to_results_dict( self.make_phot_info_dict( sci_image, templ_image ) )

        if self._omg:
            raise RuntimeError( "make_phot_info_dict failed" )

        SNLogger.debug( "Saving results..." )
        results_tab = Table(self.results_dict)
        results_tab.sort('mjd')
        results_savedir = self.ltcv_dir / 'data' / str(self.diaobj.id)
        results_savedir.mkdir( exist_ok=True, parents=True )
        results_savepath = results_savedir / f'{self.diaobj.id}_{self.band}_all.csv'
        results_tab.write(results_savepath, format='csv', overwrite=True)
        SNLogger.info(f'Results saved to {results_savepath}')

        return results_savepath

    def write_fits_file( self, data, header, savepath ):
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
            SNLogger.info( "Running sky subtraction" )
            with nvtx.annotate( "skysub", color=0xff8888 ):
                self.sky_sub_all_images()

        if self.mem_trace:
            SNLogger.info( f"After sky_subtract, memory usage = {tracemalloc.get_traced_memory()[1]/(1024**2):.2f} MB" )

        if 'get_psfs' in steps:
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

                    if 'align_and_preconvolve' in steps:
                        SNLogger.info( "...align_and_preconvolve" )
                        with nvtx.annotate( "align_and_pre_convolve", color=0x8888ff ):
                            sfftifier = self.align_and_pre_convolve( templ_image, sci_image )

                    if 'subtract' in steps:
                        SNLogger.info( "...subtract" )
                        with nvtx.annotate( "subtraction", color=0x44ccff ):
                            sfftifier.sfft_subtraction()

                    if 'find_decorrelation' in steps:
                        SNLogger.info( "...find_decorrelation" )
                        with nvtx.annotate( "find_decor", color=0xcc44ff ):
                            sfftifier.find_decorrelation()

                        SNLogger.info( "...generate variance image" )
                        with nvtx.annotate( "variance", color=0x44ccff ):
                            diff_var = sfftifier.create_variance_image()
                            mess = f"{sci_image.image.name}-{templ_image.image.name}"
                            diff_var_path = self.dia_out_dir / f"diff_var_{mess}"

                    if 'apply_decorrelation' in steps:
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

                    if self.keep_intermediate:
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

        if 'make_stamps' in steps:
            SNLogger.info( "Starting to make stamps..." )
            with nvtx.annotate( "make stamps", color=0xff8888 ):
                self._omg = False

                def log_stamp_err( x ):
                    SNLogger.error( f"do_stamps subprocess failure: {x} " )
                    self._omg = True

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
                                sci_stamp_pool.apply_async( self.do_stamps, pair, {},
                                                            callback = partial(self.save_stamp_paths,
                                                                               sci_image, templ_image),
                                                            error_callback=log_stamp_err )
                        sci_stamp_pool.close()
                        sci_stamp_pool.join()

                else:
                    for tsargs in templstamp_args:
                        partialstamp(*tsargs)

                    for sci_image in self.science_images:
                        for templ_image in self.template_images:
                            stamp_paths = self.do_stamps( sci_image, templ_image)
                            self.save_stamp_paths( sci_image, templ_image, stamp_paths )

                if self._omg:
                    raise RuntimeError( "Error making stamps." )

            SNLogger.info('...finished making stamps.')

        if self.mem_trace:
            SNLogger.info( f"After make_stamps, memory usage = {tracemalloc.get_traced_memory()[1]/(1024**2):.2f} MB" )

        lightcurve_path = None
        if 'make_lightcurve' in steps:
            with nvtx.annotate( "make_lightcurve", color=0xff8888 ):
                lightcurve_path = self.make_lightcurve()

        if self.mem_trace:
            SNLogger.info( f"After make_lightcurve, memory usage = \
                            {tracemalloc.get_traced_memory()[1]/(1024**2):.2f} MB" )

        if self.remove_temp_dir:
            self.clear_contents( self.temp_dir )

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
    parser.add_argument( '--object-collection', '--oc', required=True,
                         help='Collection of the object.  Currently only "ou2024" and "manual" supported.' )
    parser.add_argument( '--object-subset', '--os', default=None,
                         help="Collection subset.  Not used by all collections." )
    parser.add_argument( '--oid', type=int, required=True,
                         help="Object ID.  Meaning is collection-dependent." )
    parser.add_argument( '-r', '--ra', type=float, default=None,
                         help="Object RA.  By default, uses the one found for the object." )
    parser.add_argument( '-d', '--dec', type=float, default=None,
                         help="Object Dec.  By default, uses the one found for the object." )
    parser.add_argument( '-b', '--band', type=str, required=True,
                         help="Band: R062, Z087, Y106, J129, H158, F184, or K213" )
    parser.add_argument( '--image-collection', '--ic', required=True, help="Collection of the images we're using" )
    parser.add_argument( '--image-subset', '--is', default=None, help="Image collection subset" )
    parser.add_argument( '--base-path', type=str, default=None,
                         help='Base path for images.  Required for "manual_fits" image collection' )
    parser.add_argument( '-t', '--template-images', type=str, required=True,
                         help="Path to file with, per line, ( path_to_image, pointing, sca, mjd, band )" )
    parser.add_argument( '-s', '--science-images', type=str, required=True,
                         help="Path to file with, per line, ( path_to_image, pointing, sca, mjd, band )" )
    parser.add_argument( '-p', '--nprocs', type=int, default=1,
                         help="Number of process for multiprocessing steps (e.g. skysub)" )
    parser.add_argument( '-w', '--nwrite', type=int, default=5, help="Number of parallel FITS writing processes" )
    parser.add_argument( '-v', '--verbose', action='store_true', default=False, help="Show debug log info" )
    parser.add_argument( '--through-step', default='make_lightcurve',
                         help="Stop after this step; one of (see above)" )

    cfg.augment_argparse( parser )
    args = parser.parse_args( leftovers )
    cfg.parse_args( args )

    # Get the DiaObject, update the RA and Dec

    diaobjs = DiaObject.find_objects( collection=args.object_collection, subset=args.object_subset,
                                      id=args.oid, ra=args.ra, dec=args.dec )
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

    imgcol = ImageCollection.get_collection( collection=args.image_collection, subset=args.image_subset,
                                             base_path=args.base_path )

    # Create and launch the pipeline

    pipeline = Pipeline( diaobj, imgcol, args.band,
                         science_csv=args.science_images, template_csv=args.template_images,
                         nprocs=args.nprocs, nwrite=args.nwrite, verbose=args.verbose )
    pipeline( args.through_step )


# ======================================================================
if __name__ == "__main__":
    main()
