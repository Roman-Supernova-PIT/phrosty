import nvtx

import re
import pathlib
import argparse
import logging
from multiprocessing import Pool
from functools import partial

import numpy as np
import cupy as cp

import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.table
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u

from sfft.SpaceSFFTCupyFlow import SpaceSFFT_CupyFlow

from snpit_utils.logger import SNLogger
from snpit_utils.config import Config
from snappl.psf import PSF
from snappl.image import OpenUniverse2024FITSImage
from phrosty.utils import read_truth_txt, get_exptime
from phrosty.imagesubtraction import sky_subtract, stampmaker
from phrosty.photometry import ap_phot, psfmodel, psf_phot

from galsim import roman


class PipelineImage:
    """Holds a snappl.image.Image, with some other stuff the pipeline needs."""

    def __init__( self, imagepath, pointing, sca ):
        """Create a PipelineImage

        Parameters:
        -----------
           imagepath : str or Path
              A path to the image.  This will be passed on to an Image
              subclass constructor; which subclass depends on the config
              option photometry.phrosty.image_type

           pointing : str or int
              An identifier of the pointing of this image.  Used e.g. to pull PSFs.

           pipeline : phrosty.pipeline.Pipeline
              The pipeline that owns this image.

        """
        # self.psf is a object of a subclass of snappl.psf.PSF
        self.config = Config.get()
        self.temp_dir = pathlib.Path( self.config.value( 'photometry.phrosty.paths.temp_dir' ) )

        if self.config.value( 'photometry.phrosty.image_type' ) == 'ou2024fits':
            self.image = OpenUniverse2024FITSImage( imagepath, None, sca )
        else:
            raise RuntimeError( "At the moment, phrosty only works with ou2024fits images. "
                                "We hope this will change soon." )

        self.pointing = pointing

        self.decorr_psf_path = {}
        self.decorr_zptimg_path = {}
        self.decorr_diff_path = {}
        self.zpt_stamp_path = {}
        self.diff_stamp_path = {}

        self.skysub_path = None
        self.detmas_path = None
        self.skyrms = None
        self.psfobj = None
        self.psf_data = None

    def run_sky_subtract( self, mp=True ):
        # Eventually, we may not want to save the sky subtracted image, but keep
        #   it in memory.  (Reduce I/O.)  Will require SFFT changes.
        #   (This may not be practical, as it will increase memory usage a *lot*.
        #   We may still need to write files.)
        try:
            imname = self.image.name
            # HACK ALERT : we're stripping the .gz off of the end of filenames
            #   if they have them, and making sure filenames end in .fits,
            #   because that's what SFFT needs.  This can go away if we
            #   refactor to pass data.
            if imname[-3:] == '.gz':
                imname = imname[:-3]
            if imname[-5:] != '.fits':
                imname = f'{imname}.fits'
            if mp:
                SNLogger.multiprocessing_replace()
            SNLogger.debug( f"run_sky_subtract on {imname}" )
            self.skysub_path = self.temp_dir / f"skysub_{imname}"
            self.detmask_path = self.temp_dir / f"detmask_{imname}"
            self.skyrms = sky_subtract( self.image.path, self.skysub_path, self.detmask_path, temp_dir=self.temp_dir,
                                        force=self.config.value( 'photometry.phrosty.force_sky_subtract' ) )
            SNLogger.debug( f"...done running sky subtraction on {self.image.name}" )
            return ( self.skysub_path, self.detmask_path, self.skyrms )
        except Exception as ex:
            SNLogger.exception( ex )
            raise

    def save_sky_subtract_info( self, info ):
        SNLogger.debug( f"Saving sky_subtract info for path {info[0]}" )
        self.skysub_path = info[0]
        self.detmask_path = info[1]
        self.skyrms = info[2]

    def get_psf( self, ra, dec, dump_file=False ):
        """Get the at the right spot on the image.

        Parameters
        ----------
          ra, dec : float
             The coordinates in decimal degrees where we want the PSFD.

          dump_file : bool
             If True, write out the psf as a FITS file in dia_out_dir
             (for diagnostic purposes; these files are not read again
             interally by the pipeline).

        """

        # TODO: right now snappl.psf.PSF.get_psf_object just
        #   passes the keyword arguments on to whatever makes
        #   the psf... and it's different for each type of
        #   PSF.  We need to fix that... somehow....

        if self.psfobj is None:
            psftype = self.config.value( 'photometry.phrosty.psf.type' )
            self.psfobj = PSF.get_psf_object( psftype, pointing=self.pointing, sca=self.image.sca )

        wcs = self.image.get_wcs()
        x, y = wcs.world_to_pixel( ra, dec )
        stamp = self.psfobj.get_stamp( x, y )
        if dump_file:
            outfile = pathlib.Path( self.config.value( "photometry.phrosty.paths.dia_out_dir" ) )
            outfile = outfile / f"psf_{self.image.name}.fits"
            fits.writeto( outfile, stamp, overwrite=True )
        return stamp

    def keep_psf_data( self, psf_data ):
        self.psf_data = psf_data


class Pipeline:
    def __init__( self, object_id, ra, dec, band, science_images, template_images, nprocs=1, nwrite=5,
                  nuke_temp_dir=False, verbose=False ):

        """Create the a pipeline object.

        Parameters
        ----------
           object_id: int

           ra, dec: float
             Position of transient in decimal degrees

           band: str
             One of R062, Z087, Y106, J129, H158, F184, K213

           science_images: list of tuple
               ( path_to_image, pointing, sca )

           template_images: list of tuple
               ( path_to_image, pointing, sca )

           nprocs: int, default 1
             Number of cpus for the CPU multiprocessing segments of the pipeline.
             (GPU segments will run a single process.)

           nwrite: int, default 5
             Number of asynchronous FITS writer processes.


        """

        SNLogger.setLevel( logging.DEBUG if verbose else logging.INFO )
        self.config = Config.get()
        self.image_base_dir = pathlib.Path( self.config.value( 'photometry.phrosty.paths.image_base_dir' ) )
        self.dia_out_dir = pathlib.Path( self.config.value( 'photometry.phrosty.paths.dia_out_dir' ) )
        self.ltcv_dir = pathlib.Path( self.config.value( 'photometry.phrosty.paths.ltcv_dir' ) )

        self.object_id = object_id
        self.ra = ra
        self.dec = dec
        self.band = band
        self.science_images = ( [ PipelineImage( self.image_base_dir / ppsmb[0], ppsmb[1], ppsmb[2] )
                                  for ppsmb in science_images if ppsmb[4] == self.band ] )
        self.template_images = ( [ PipelineImage( self.image_base_dir / ppsmb[0], ppsmb[1], ppsmb[2] )
                                   for ppsmb in template_images if ppsmb[4] == self.band ] )
        self.nprocs = nprocs
        self.nwrite = nwrite
        self.nuke_temp_dir = nuke_temp_dir
        if self.nuke_temp_dir:
            SNLogger.warning( "nuke_temp_dir not implemented" )


    def sky_sub_all_images( self ):
        # Currently, this writes out a bunch of FITS files.  Further refactoring needed
        #   to support more general image types.
        all_imgs = self.science_images.copy()     # shallow copy
        all_imgs.extend( self.template_images )

        def log_error( img, x ):
            SNLogger.error( f"Sky subtraction subprocess failure: {x} for image {img.image.path}" )

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
        all_imgs = self.science_images.copy()     # shallow copy
        all_imgs.extend( self.template_images )

        if self.nprocs > 1:
            with Pool( self.nprocs ) as pool:
                for img in all_imgs:
                    # callback_partial = partial( img.save_psf_path, all_imgs )
                    pool.apply_async( img.get_psf, (self.ra, self.dec), {},
                                      img.keep_psf_data,
                                      lambda x: SNLogger.error( f"get_psf subprocess failure: {x}" ) )
                pool.close()
                pool.join()
        else:
            for img in all_imgs:
                img.keep_psf_data( img.get_psf(self.ra, self.dec) )


    def align_and_pre_convolve(self, templ_image, sci_image ):
        """Align and pre convolve a single template/science pair.

        Parameters
        ----------
          sci_image: phrosty.Image
            The science (new) image.

          templ_image: phrosty.Image
            The template (ref) image that will be subtracted from sci_image.

        """

        with fits.open( sci_image.skysub_path ) as hdul:
            hdr_sci = hdul[0].header
            data_sci = cp.array( np.ascontiguousarray(hdul[0].data.T), dtype=cp.float64 )
        with fits.open( templ_image.skysub_path ) as hdul:
            hdr_templ = hdul[0].header
            data_templ = cp.array( np.ascontiguousarray(hdul[0].data.T), dtype=cp.float64 )

        sci_psf = cp.ascontiguousarray( cp.array( sci_image.psf_data.T, dtype=cp.float64 ) )
        templ_psf = cp.ascontiguousarray( cp.array( templ_image.psf_data.T, dtype=cp.float64 ) )

        with fits.open( sci_image.detmask_path ) as hdul:
            sci_detmask = cp.array( np.ascontiguousarray( hdul[0].data.T ) )

        with fits.open( templ_image.detmask_path ) as hdul:
            templ_detmask = cp.array( np.ascontiguousarray( hdul[0].data.T ) )

        sfftifier = SpaceSFFT_CupyFlow(
            hdr_sci, hdr_templ,
            sci_image.skyrms, templ_image.skyrms,
            data_sci, data_templ,
            sci_detmask, templ_detmask,
            sci_psf, templ_psf
        )
        sfftifier.resampling_image_mask_psf()
        sfftifier.cross_convolution()

        return sfftifier

    def phot_at_coords( self, img, psf, pxcoords=(50, 50), ap_r=4 ):
        """Do photometry at forced set of pixel coordinates."""

        forcecoords = Table([[float(pxcoords[0])], [float(pxcoords[1])]], names=["x", "y"])
        init = ap_phot(img, forcecoords, ap_r=ap_r)
        init['flux_init'] = init['aperture_sum']
        final = psf_phot(img, psf, init, forced_phot=True)

        flux = final['flux_fit'][0]
        flux_err = final['flux_err'][0]
        mag = -2.5 * np.log10(final["flux_fit"][0])
        mag_err = (2.5 / np.log(10)) * np.abs(final["flux_err"][0] / final["flux_fit"][0])

        results_dict = {
                        'flux_fit': flux,
                        'flux_fit_err': flux_err,
                        'mag_fit': mag,
                        'mag_fit_err': mag_err
                        }

        return results_dict

    def get_stars(self, truthpath, nx=4088, ny=4088, transform=False, wcs=None):
        """Get the stars in the science images.

        Optional to transform to another WCS.
        """
        truth_tab = read_truth_txt(path=truthpath)
        truth_tab['mag'].name = 'mag_truth'
        truth_tab['flux'].name = 'flux_truth'

        if transform:
            assert wcs is not None, 'You need to provide a WCS to transform to!'
            truth_tab['x'].name, truth_tab['y'].name = 'x_orig', 'y_orig'
            worldcoords = SkyCoord(ra=truth_tab['ra'] * u.deg, dec=truth_tab['dec'] * u.deg)
            x, y = skycoord_to_pixel(worldcoords, wcs)
            truth_tab['x'] = x
            truth_tab['y'] = y

        if not transform:
            truth_tab['x'] -= 1
            truth_tab['y'] -= 1

        idx = np.where(truth_tab['obj_type'] == 'star')[0]
        stars = truth_tab[idx]
        stars = stars[np.logical_and(stars["x"] < nx, stars["x"] > 0)]
        stars = stars[np.logical_and(stars["y"] < ny, stars["y"] > 0)]

        return stars

    def get_galsim_values(self):
        exptime = get_exptime(self.band)
        area_eff = roman.collecting_area
        gs_zpt = roman.getBandpasses()[self.band].zeropoint

        return {'exptime': exptime, 'area_eff': area_eff, 'gs_zpt': gs_zpt}

    def get_zpt(self, zptimg, psf, band, stars, ap_r=4, ap_phot_only=False,
                zpt_plot=None, oid=None, sci_pointing=None, sci_sca=None):

        # TODO : Need to move this code all over into snappl Image.  It sounds like
        #   for Roman images we may have to do our own zeropoints (which is what
        #   is happening here), but for actual Roman we're going to use
        #   calibration information we get from elsewhere, so we don't want to bake
        #   doing the calibration into the pipeline as we do here.
        #
        # Also Issue #70

        """Get the zeropoint based on the stars."""

        # First, need to do photometry on the stars.
        init_params = ap_phot(zptimg, stars, ap_r=ap_r)
        init_params['flux_init'] = init_params['aperture_sum']
        final_params = psf_phot(zptimg, psf, init_params, forced_phot=True)

        # Do not need to cross match. Can just merge tables because they
        # will be in the same order.
        photres = astropy.table.join(stars, init_params, keys=['object_id', 'ra', 'dec', 'realized_flux',
                                                               'flux_truth', 'mag_truth', 'obj_type'])
        if not ap_phot_only:
            photres = astropy.table.join(photres, final_params, keys=['id'])

        # Get the zero point.
        galsim_vals = self.get_galsim_values()
        star_fit_mags = -2.5 * np.log10(photres['flux_fit'])
        star_truth_mags = ( -2.5 * np.log10(photres['flux_truth']) + galsim_vals['gs_zpt']
                            + 2.5 * np.log10(galsim_vals['exptime'] * galsim_vals['area_eff']) )

        # Eventually, this should be a S/N cut, not a mag cut.
        zpt_mask = np.logical_and(star_truth_mags > 20, star_truth_mags < 23)
        zpt = np.nanmedian(star_truth_mags[zpt_mask] - star_fit_mags[zpt_mask])

        if zpt_plot is not None:
            assert oid is not None, 'If zpt_plot=True, oid must be provided.'
            assert sci_pointing is not None, 'If zpt_plot=True, sci_pointing must be provided.'
            assert sci_sca is not None, 'If zpt_plot=True, sci_sca must be provided.'

            savedir = self.ltcv_dir / f'figs/{oid}/zpt_plots'
            savedir.mkdir(parents=True, exist_ok=True)
            savepath = savedir / f'zpt_stars_{band}_{sci_pointing}_{sci_sca}.png'

            plt.figure(figsize=(8, 8))
            yaxis = star_fit_mags + zpt - star_truth_mags

            plt.plot(star_truth_mags, yaxis, marker='o', linestyle='')
            plt.axhline(0, linestyle='--', color='k')
            plt.xlabel('Truth mag')
            plt.ylabel('Fit mag - zpt + truth mag')
            plt.title(f'{band} {sci_pointing} {sci_sca}')
            plt.savefig(savepath, dpi=300, bbox_inches='tight')
            plt.close()

            SNLogger.info(f'zpt debug plot saved to {savepath}')

            # savepath = os.path.join(savedir, f'hist_truth-fit_{band}_{sci_pointing}_{sci_sca}.png')
            # plt.hist(star_truth_mags[zpt_mask] - star_fit_mags[zpt_mask])
            # plt.title(f'{band} {sci_pointing} {sci_sca}')
            # plt.xlabel('star_truth_mags[zpt_mask] - star_fit_mags[zpt_mask]')
            # plt.savefig(savepath, dpi=300, bbox_inches='tight')
            # plt.close()

        return zpt

    def make_phot_info_dict( self, sci_image, templ_image, ap_r=4 ):
        # Do photometry on stamp because it will read faster
        diff_img_stamp_path = sci_image.diff_stamp_path[ templ_image.image.name ]

        results_dict = {}
        results_dict['sci_name'] = sci_image.image.name
        results_dict['templ_name'] = templ_image.image.name
        results_dict['success'] = False
        results_dict['ra'] = self.ra
        results_dict['dec'] = self.dec
        results_dict['mjd'] = sci_image.image.mjd
        results_dict['filter'] = self.band
        results_dict['pointing'] = sci_image.pointing
        results_dict['sca'] = sci_image.image.sca
        results_dict['template_pointing'] = templ_image.pointing
        results_dict['template_sca'] = templ_image.image.sca

        if diff_img_stamp_path.is_file():
            # Load in the difference image stamp.
            with fits.open(diff_img_stamp_path) as diff_hdu:
                diffimg = diff_hdu[0].data
                wcs = WCS(diff_hdu[0].header)

            # Load in the decorrelated PSF.
            psfpath = sci_image.decorr_psf_path[ templ_image.image.name ]
            with fits.open( psfpath ) as hdu:
                psf = psfmodel( hdu[0].data )
            coord = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg)
            pxcoords = skycoord_to_pixel(coord, wcs)
            results_dict.update( self.phot_at_coords(diffimg, psf, pxcoords=pxcoords, ap_r=ap_r) )

            # Get the zero point from the decorrelated, convolved science image.
            # First, get the table of known stars.

            # TODO -- take this galsim-specific code out, move it to a separate module.  Define a general
            #  zeropointing interface, of which the galsim-speicifc one will be one instance
            truthpath = str( self.image_base_dir /
                             f'RomanTDS/truth/{self.band}/{sci_image.pointing}/'
                             f'Roman_TDS_index_{self.band}_{sci_image.pointing}_{sci_image.image.sca}.txt' )
            stars = self.get_stars(truthpath)
            # Now, calculate the zero point based on those stars.
            zptimg_path = sci_image.decorr_zptimg_path[ templ_image.image.name ]
            with fits.open(zptimg_path) as hdu:
                zptimg = hdu[0].data
            zpt = self.get_zpt(zptimg, psf, self.band, stars, oid=self.object_id,
                               sci_pointing=sci_image.pointing, sci_sca=sci_image.image.sca)

            # Add additional info to the results dictionary so it can be merged into a nice file later.
            results_dict['zpt'] = zpt
            results_dict['success'] = True

        else:
            SNLogger.warning( f"Post-processed image files for "
                              f"{self.band}_{sci_image.pointing}_{sci_image.image.sca}-"
                              f"{self.band}_{templ_image.pointing}_{templ_image.image.sca} "
                              f"do not exist.  Skipping." )
            results_dict['zpt'] = np.nan
            results_dict['flux_fit'] = np.nan
            results_dict['flux_fit_err'] = np.nan
            results_dict['mag_fit'] = np.nan
            results_dict['mag_fit_err'] = np.nan

        return results_dict

    def add_to_results_dict( self, one_pair ):
        for key, arr in self.results_dict.items():
            arr.append( one_pair[ key ] )

    def save_stamp_paths( self, sci_image, templ_image, paths ):
        sci_image.zpt_stamp_path[ templ_image.image.name ] = paths[0]
        sci_image.diff_stamp_path[ templ_image.image.name ] = paths[1]

    def do_stamps( self, sci_image, templ_image ):

        zptname = sci_image.decorr_zptimg_path[ templ_image.image.name ]
        zpt_stampname = stampmaker( self.ra, self.dec, np.array([100, 100]),
                                    zptname,
                                    savedir=self.dia_out_dir,
                                    savename=f"stamp_{zptname.name}" )

        diffname = sci_image.decorr_diff_path[ templ_image.image.name ]
        diff_stampname = stampmaker( self.ra, self.dec, np.array([100, 100]),
                                 diffname,
                                 savedir=self.dia_out_dir,
                                 savename=f"stamp_{diffname.name}" )

        SNLogger.info(f"Decorrelated stamp path: {pathlib.Path( diff_stampname )}")
        SNLogger.info(f"Zpt image stamp path: {pathlib.Path( zpt_stampname )}")

        return pathlib.Path( zpt_stampname ), pathlib.Path( diff_stampname )

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
            'flux_fit': [],
            'flux_fit_err': [],
            'mag_fit': [],
            'mag_fit_err': [],
        }

        if self.nprocs > 1:
            with Pool( self.nprocs ) as pool:
                for sci_image in self.science_images:
                    for templ_image in self.template_images:
                        pool.apply_async( self.make_phot_info_dict, (sci_image, templ_image), {},
                                          self.add_to_results_dict,
                                          error_callback=lambda x: SNLogger.error( f"make_phot_info_dict "
                                                                                   f"subprocess failure: {x}" )
                                         )
                pool.close()
                pool.join()
        else:
            for i, sci_image in enumerate( self.science_images ):
                SNLogger.debug( f"Doing science image {i} of {len(self.science_images)}" )
                for templ_image in self.template_images:
                    self.add_to_results_dict( self.make_phot_info_dict( sci_image, templ_image ) )

        results_tab = Table(self.results_dict)
        results_tab.sort('mjd')
        results_savedir = self.ltcv_dir / 'data' / str(self.object_id)
        results_savedir.mkdir( exist_ok=True, parents=True )
        results_savepath = results_savedir / f'{self.object_id}_{self.band}_all.csv'
        results_tab.write(results_savepath, format='csv', overwrite=True)
        SNLogger.info(f'Results saved to {results_savepath}')

    def write_fits_file( self, data, header, savepath ):
        fits.writeto( savepath, data, header=header, overwrite=True )

    def __call__( self, through_step=None ):
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

        if 'get_psfs' in steps:
            SNLogger.info( "Getting PSFs" )
            with nvtx.annotate( "getpsfs", color=0xff8888 ):
                self.get_psfs()

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

                    if 'apply_decorrelation' in steps:
                        mess = ( f"{self.band}_{sci_image.pointing}_{sci_image.image.sca}_-"
                                 f"_{self.band}_{templ_image.pointing}_{templ_image.image.sca}.fits" )
                        decorr_psf_path = self.dia_out_dir / f"decorr_psf_{mess}"
                        decorr_zptimg_path = self.dia_out_dir / f"decorr_zptimg_{mess}"
                        decorr_diff_path = self.dia_out_dir / f"decorr_diff_{mess}"
                        for img, savepath, hdr in zip(
                                [ sfftifier.PixA_DIFF_GPU, sfftifier.PixA_Ctarget_GPU, sfftifier.PSF_target_GPU ],
                                [ decorr_diff_path,        decorr_zptimg_path,         decorr_psf_path ],
                                [ sfftifier.hdr_target,    sfftifier.hdr_target,       None ]
                        ):
                            with nvtx.annotate( "apply_decor", color=0xccccff ):
                                SNLogger.info( f"...apply_decor to {savepath}" )
                                decorimg = sfftifier.apply_decorrelation( img )
                            with nvtx.annotate( "submit writefits", color=0xff8888 ):
                                SNLogger.info( f"...writefits {savepath}" )
                                fits_writer_pool.apply_async( self.write_fits_file,
                                                              ( cp.asnumpy( decorimg ).T, hdr, savepath ), {},
                                                              error_callback=partial(log_fits_write_error, savepath) )
                        sci_image.decorr_psf_path[ templ_image.image.name ]= decorr_psf_path
                        sci_image.decorr_zptimg_path[ templ_image.image.name ]= decorr_zptimg_path
                        sci_image.decorr_diff_path[ templ_image.image.name ]= decorr_diff_path

                        SNLogger.info( f"DONE processing {sci_image.image.name} minus {templ_image.image.name}" )

            SNLogger.info( "Waiting for FITS writer processes to finish" )
            with nvtx.annotate( "fits_write_wait", color=0xff8888 ):
                fits_writer_pool.close()
                fits_writer_pool.join()
            SNLogger.info( "...FITS writer processes done." )

        if 'make_stamps' in steps:
            SNLogger.info( "Starting to make stamps..." )
            with nvtx.annotate( "make stamps", color=0xff8888 ):
                if self.nwrite > 1:
                    partialstamp = partial(stampmaker, self.ra, self.dec, np.array([100, 100]))
                    # template path, savedir, savename
                    templstamp_args = ( (ti, self.dia_out_dir, f'stamp_{ti}') for ti in self.template_images)

                    with Pool( self.nwrite ) as templ_stamp_pool:
                        templ_stamp_pool.starmap_async( partialstamp, templstamp_args )
                        templ_stamp_pool.close()
                        templ_stamp_pool.join()

                    with Pool( self.nwrite ) as sci_stamp_pool:
                        for sci_image in self.science_images:
                            for templ_image in self.template_images:
                                pair = (sci_image, templ_image)
                                sci_stamp_pool.apply_async( self.do_stamps, pair, {},
                                                            callback = partial(self.save_stamp_paths,
                                                                               sci_image, templ_image),
                                                            error_callback=partial(SNLogger.error,
                                                                                   "do_stamps subprocess failure: {x}")
                                                           )

                        sci_stamp_pool.close()
                        sci_stamp_pool.join()

                else:
                    for templ_image in self.template_images:
                        stamp_name = stampmaker( self.ra, self.dec, np.array([100, 100]), templ_image.image.path,
                                                savedir=self.dia_out_dir, savename=f"stamp_{templ_image.image.name}" )

                    for sci_image in self.science_images:
                        for templ_image in self.template_images:
                            zptname = sci_image.decorr_zptimg_path[ templ_image.image.name ]
                            diffname = sci_image.decorr_diff_path[ templ_image.image.name ]
                            stamp_name = stampmaker( self.ra, self.dec, np.array([100, 100]), zptname,
                                                    savedir=self.dia_out_dir, savename=f"stamp_{zptname.name}" )
                            sci_image.zpt_stamp_path[ templ_image.image.name ] = pathlib.Path( stamp_name )
                            stamp_name = stampmaker( self.ra, self.dec, np.array([100, 100]), diffname,
                                                    savedir=self.dia_out_dir, savename=f"stamp_{diffname.name}" )
                            sci_image.diff_stamp_path[ templ_image.image.name ] = pathlib.Path( stamp_name )

            SNLogger.info('...finished making stamps.')

        if 'make_lightcurve' in steps:
            SNLogger.info( "Making lightcurve" )
            with nvtx.annotate( "make_lightcurve", color=0xff8888 ):
                self.make_lightcurve()


                # ======================================================================

def main():
    # Run one arg pass just to get the config file, so we can augment
    #   the full arg parser later with config options
    configparser = argparse.ArgumentParser( add_help=False )
    configparser.add_argument( '-c', '--config-file', required=True, help="Location of the .yaml config file" )
    args, leftovers = configparser.parse_known_args()

    cfg = Config.get( args.config_file, setdefault=True )

    parser = argparse.ArgumentParser()
    # Put in the config_file argument, even though it will never be found, so it shows up in help
    parser.add_argument( '-c', '--config-file', help="Location of the .yaml config file" )
    parser.add_argument( '--oid', type=int, required=True, help="Object ID" )
    parser.add_argument( '-r', '--ra', type=float, required=True, help="Object RA" )
    parser.add_argument( '-d', '--dec', type=float, required=True, help="Object Dec" )
    parser.add_argument( '-b', '--band', type=str, required=True,
                         help="Band: R062, Z087, Y106, J129, H158, F184, or K213" )
    parser.add_argument( '-t', '--template-images', type=str, required=True,
                         help="Path to file with, per line, ( path_to_image, pointing, sca )" )
    parser.add_argument( '-s', '--science-images', type=str, required=True,
                         help="Path to file with, per line, ( path_to_image, pointing, sca )" )
    parser.add_argument( '-p', '--nprocs', type=int, default=1,
                         help="Number of process for multiprocessing steps (e.g. skysub)" )
    parser.add_argument( '-w', '--nwrite', type=int, default=5, help="Number of parallel FITS writing processes" )
    parser.add_argument( '-v', '--verbose', action='store_true', default=False, help="Show debug log info" )
    parser.add_argument( '--through-step', default='make_lightcurve',
                         help="Stop after this step; one of (see above)" )

    cfg.augment_argparse( parser )
    args = parser.parse_args( leftovers )
    cfg.parse_args( args )

    science_images = []
    template_images = []
    for infile, imlist in zip( [ args.science_images, args.template_images ],
                               [ science_images, template_images ] ):
        with open( infile ) as ifp:
            hdrline = ifp.readline()
            if not re.search( r"^\s*path\s+pointing\s+sca\s+mjd\s+band\s*$", hdrline ):
                raise ValueError( f"First line of list file {infile} didn't match what was expected." )
            for line in ifp:
                img, point, sca, mjd, band = line.split()
                imlist.append( ( pathlib.Path(img), int(point), int(sca), float(mjd), band ) )

    pipeline = Pipeline( args.oid, args.ra, args.dec, args.band, science_images, template_images,
                         nprocs=args.nprocs, nwrite=args.nwrite, nuke_temp_dir=False, verbose=args.verbose )
    pipeline( args.through_step )


# ======================================================================
if __name__ == "__main__":
    main()
