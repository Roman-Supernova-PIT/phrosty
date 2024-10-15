import nvtx

import os
import pathlib
import argparse
import logging
from multiprocessing import Pool
from functools import partial

import numpy as np
import cupy as cp

import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units

from sfft.SpaceSFFTCupyFlow import SpaceSFFT_CupyFlow

from phrosty.utils import set_logger, read_truth_txt, get_exptime
from phrosty.imagesubtraction import sky_subtract, get_imsim_psf, stampmaker
from phrosty.photometry import ap_phot, psfmodel, psf_phot

from galsim import roman

class Image:
    def __init__( self, path, pointing, sca, mjd, pipeline ):
        self.pipeline = pipeline
        self.image_path = pathlib.Path( path )
        self.image_name = self.image_path.name
        if self.image_name[-3:] == '.gz':
            self.image_name = self.image_name[:-3]
        if self.image_name[-5:] != '.fits':
            raise ValueError( f"Image name {self.image_name} doesn't end in .fits, I don't know how to cope." )
        self.basename = self.image_name[:-5]
        self.pointing = pointing
        self.sca = sca
        self.mjd = mjd
        self.psf_path = None
        self.detect_mask_path = None
        self.skyrms = None
        self.skysub_path = None

        self.decorr_psf_path = {}
        self.decorr_zptimg_path = {}
        self.decorr_diff_path = {}
        self.zpt_stamp_path = {}
        self.diff_stamp_path = {}

    def run_sky_subtract( self ):
        self.skysub_path = self.pipeline.temp_dir / f"skysub_{self.image_name}"
        self.detmask_path = self.pipeline.temp_dir / f"detmask_{self.image_name}"
        self.skyrms = sky_subtract( self.image_path, self.skysub_path, self.detmask_path,
                                    temp_dir=self.pipeline.temp_dir, force=self.pipeline.force_sky_subtract )
        return ( self.skysub_path, self.detmask_path, self.skyrms )

    def save_sky_subtract_info( self, info ):
        self.skysub_path = info[0]
        self.detmask_path = info[1]
        self.skyrms = info[2]


    def run_get_imsim_psf( self ):
        psf_path = self.pipeline.temp_dir / f"psf_{self.image_name}"
        get_imsim_psf( self.image_path, self.pipeline.ra, self.pipeline.dec, self.pipeline.band,
                       self.pointing, self.sca,
                       psf_path=psf_path, config_yaml_file=self.pipeline.galsim_config_file )
        return psf_path

    def save_psf_path( self, psf_path ):
        self.psf_path = psf_path


class Pipeline:
    def __init__( self, object_id, ra, dec, band, science_images, template_images, ncpus=1,
                  temp_dir='/phrosty_temp', out_dir='/dia_out_dir', ltcv_dir='/lc_out_dir', galsim_config_file=None,
                  force_sky_subtract=False, nuke_temp_dir=False, verbose=False ):
        """
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

           ncpus: int, default 1
             Number of cpus for the CPU multiprocessing segments of the pipeline.
             (GPU segments will run a single process.)


        """

        if galsim_config_file is None:
            raise RuntimeError( "Gotta give me galsim_config_file" )
        self.galsim_config_file = galsim_config_file

        self.object_id = object_id
        self.ra = ra
        self.dec = dec
        self.band = band
        self.science_images = [ Image( ppsm[0], ppsm[1], ppsm[2], ppsm[3], self ) for ppsm in science_images if self.band in ppsm[0].name ]
        self.template_images = [ Image( ppsm[0], ppsm[1], ppsm[2], ppsm[3], self ) for ppsm in template_images if self.band in ppsm[0].name ]
        self.ncpus = ncpus
        self.temp_dir = pathlib.Path(temp_dir )
        self.out_dir = pathlib.Path( out_dir )
        self.ltcv_dir = pathlib.Path( ltcv_dir ) if ltcv_dir is not None else self.out_dir
        self.nuke_temp_dir = nuke_temp_dir
        self.force_sky_subtract = force_sky_subtract

        self.logger = set_logger( 'phrosty', 'phrosty' )
        self.logger.setLevel( logging.DEBUG if verbose else logging.INFO )

        if self.nuke_temp_dir:
            self.logger.warning( "nuke_temp_dir not implemented" )

    def sky_sub_all_images( self ):
        all_imgs = self.science_images.copy()     # shallow copy
        all_imgs.extend( self.template_images )

        if self.ncpus > 1:
            with Pool( self.ncpus ) as pool:
                for img in all_imgs:
                    pool.apply_async( img.run_sky_subtract, (), {},
                                      callback=img.save_sky_subtract_info,
                                      error_callback=lambda x: logger.error( f"Sky subtraction subprocess failure: {x}" ) )
                pool.close()
                pool.join()
        else:
            for img in all_imgs:
                img.save_sky_subtract_info( img.run_sky_subtract() )



    def get_psfs( self ):
        all_imgs = self.science_images.copy()     # shallow copy
        all_imgs.extend( self.template_images )

        if self.ncpus > 1:
            with Pool( self.ncpus ) as pool:
                for img in all_imgs:
                    callback_partial = partial( img.save_psf_path, all_imgs )
                    pool.apply_async( img.run_get_imsim_psf, (), {},
                                      img.save_psf_path,
                                      lambda x: logger.error( f"get_imsim_psf subprocess failure: {x}" ) )
                pool.close()
                pool.join()
        else:
            for img in all_imgs:
                img.save_psf_path( img.run_get_imsim_psf() )




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

        with fits.open( sci_image.psf_path ) as hdul:
            sci_psf = cp.array( np.ascontiguousarray( hdul[0].data.T ), dtype=cp.float64 )

        with fits.open( templ_image.psf_path ) as hdul:
            templ_psf = cp.array( np.ascontiguousarray( hdul[0].data.T ), dtype=cp.float64 )

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
        """
        Do photometry at forced set of pixel coordinates.
        """

        forcecoords = Table([[float(pxcoords[0])], [float(pxcoords[1])]], names=["x", "y"])
        init = ap_phot(img, forcecoords, ap_r=ap_r)
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
        """
        Get the stars in the science images.

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
        """
        Get the zeropoint based on the stars.
        """

        # First, need to do photometry on the stars.
        init_params = ap_phot(zptimg, stars, ap_r=ap_r)
        final_params = psf_phot(zptimg, psf, init_params, forced_phot=True)

        # Do not need to cross match. Can just merge tables because they
        # will be in the same order.
        photres = astropy.table.join(stars,init_params,keys=['object_id','ra','dec','realized_flux',
                                                             'flux_truth','mag_truth','obj_type'])
        if not ap_phot_only:
            photres = astropy.table.join(photres,final_params,keys=['id'])

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

            plt.figure(figsize=(8,8))
            yaxis = star_fit_mags + zpt - star_truth_mags

            plt.plot(star_truth_mags,yaxis,marker='o',linestyle='')
            plt.axhline(0,linestyle='--',color='k')
            plt.xlabel('Truth mag')
            plt.ylabel('Fit mag - zpt + truth mag')
            plt.title(f'{band} {sci_pointing} {sci_sca}')
            plt.savefig(savepath,dpi=300,bbox_inches='tight')
            plt.close()

            self.logger.info(f'zpt debug plot saved to {savepath}')

            # savepath = os.path.join(savedir,f'hist_truth-fit_{band}_{sci_pointing}_{sci_sca}.png')
            # plt.hist(star_truth_mags[zpt_mask] - star_fit_mags[zpt_mask])
            # plt.title(f'{band} {sci_pointing} {sci_sca}')
            # plt.xlabel('star_truth_mags[zpt_mask] - star_fit_mags[zpt_mask]')
            # plt.savefig(savepath,dpi=300,bbox_inches='tight')
            # plt.close()

        return zpt

    def make_phot_info_dict( self, sci_image, templ_image, ap_r=4 ):
        # Do photometry on stamp because it will read faster
        diff_img_stamp_path = sci_image.diff_stamp_path[ templ_image.image_name ]

        results_dict = {}
        results_dict['sci_name'] = sci_image.image_name
        results_dict['templ_name'] = templ_image.image_name
        results_dict['success'] = False
        results_dict['ra'] = self.ra
        results_dict['dec'] = self.dec
        results_dict['mjd'] = sci_image.mjd
        results_dict['filter'] = self.band
        results_dict['pointing'] = sci_image.pointing
        results_dict['sca'] = sci_image.sca
        results_dict['template_pointing'] = templ_image.pointing
        results_dict['template_sca'] = templ_image.sca

        if diff_img_stamp_path.is_file():
            # Load in the difference image stamp.
            with fits.open(diff_img_stamp_path) as diff_hdu:
                diffimg = diff_hdu[0].data
                wcs = WCS(diff_hdu[0].header)

            # Load in the decorrelated PSF.
            psfpath = sci_image.decorr_psf_path[ templ_image.image_name ]
            with fits.open( psfpath ) as hdu:
                psf = psfmodel( hdu[0].data )
            coord = SkyCoord(ra=self.ra * astropy.units.deg, dec=self.dec * astropy.units.deg)
            pxcoords = skycoord_to_pixel(coord,wcs)
            results_dict.update( self.phot_at_coords(diffimg, psf, pxcoords=pxcoords, ap_r=ap_r) )

            # Get the zero point from the decorrelated, convolved science image.
            # First, get the table of known stars.

            # TODO -- take this galsim-specific code out, move it to a separate module.  Define a general
            #  zeropointing interface, of which the galsim-speicifc one will be one instance
            truthpath = os.path.join(os.getenv('SIMS_DIR'),
                                     ( f'RomanTDS/truth/{self.band}/{sci_image.pointing}/'
                                       f'Roman_TDS_index_{self.band}_{sci_image.pointing}_{sci_image.sca}.txt') )
            stars = self.get_stars(truthpath)
            # Now, calculate the zero point based on those stars.
            zptimg_path = sci_image.decorr_zptimg_path[ templ_image.image_name ]
            with fits.open(zptimg_path) as hdu:
                zptimg = hdu[0].data
            zpt = self.get_zpt(zptimg, psf, self.band, stars, oid=self.object_id,
                               sci_pointing=sci_image.pointing, sci_sca=sci_image.sca)

            # Add additional info to the results dictionary so it can be merged into a nice file later.
            results_dict['zpt'] = zpt
            results_dict['success'] = True

        else:
            self.logger.warning( f"Post-processed image files for {self.band}_{sci_image.pointing}_{sci_image.sca}-"
                                 f"{self.band}_{templ_image.pointing}_{templ_image.sca} do not exist.  Skipping." )
            results_dict['zpt'] = np.nan
            results_dict['flux_fit'] = np.nan
            results_dict['flux_fit_err'] = np.nan
            results_dict['mag_fit'] = np.nan
            results_dict['mag_fit_err'] = np.nan

        return results_dict

    def add_to_results_dict( self, one_pair ):
        for key, arr in self.results_dict.items():
            arr.append( one_pair[ key ] )


    def make_lightcurve( self ):
        self.logger.info( "Making lightcurve." )

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

        if self.ncpus > 1:
            with Pool( self.ncpus ) as pool:
                for sci_image in self.science_images:
                    for templ_image in self.template_images:
                        pool.apply_async( self.make_phot_info_dict, (sci_image, templ_image), {},
                                          self.add_to_results_dict,
                                          error_callback=lambda x: logger.error( "Make_phot_info_dict subprocess failure: {x}" )
                                         )
                pool.close()
                pool.join()
        else:
            for i, sci_image in enumerate( self.science_images ):
                self.logger.debug( f"Doing science image {i} of {len(self.science_images)}" )
                for templ_image in self.template_images:
                    self.add_to_results_dict( self.make_phot_info_dict( sci_image, templ_image ) )

        results_tab = Table(self.results_dict)
        results_tab.sort('mjd')
        results_savedir = self.ltcv_dir / 'data' / str(self.object_id)
        results_savedir.mkdir( exist_ok=True, parents=True )
        results_savepath = results_savedir / f'{self.object_id}_{self.band}_all.csv'
        results_tab.write(results_savepath, format='csv', overwrite=True)
        self.logger.info(f'Results saved to {results_savepath}')



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
            with nvtx.annotate( "skysub", color=0xff8888 ):
                self.sky_sub_all_images()

        if 'get_psfs' in steps:
            with nvtx.annotate( "getpsfs", color=0xff8888 ):
                self.get_psfs()

        for templ_image in self.template_images:
            for sci_image in self.science_images:
                sfftifier = None

                if 'align_and_preconvolve' in steps:
                    with nvtx.annotate( "align_and_pre_convolve", color=0x8888ff ):
                        sfftifier = self.align_and_pre_convolve( templ_image, sci_image )

                if 'subtract' in steps:
                    with nvtx.annotate( "subtraction", color=0x44ccff ):
                        sfftifier.sfft_subtraction()

                if 'find_decorrelation' in steps:
                    with nvtx.annotate( "find_decor", color=0xcc44ff ):
                        sfftifier.find_decorrelation()

                if 'apply_decorrelation' in steps:
                    mess = f"{self.band}_{sci_image.pointing}_{sci_image.sca}_-_{self.band}_{templ_image.pointing}_{templ_image.sca}.fits"
                    decorr_psf_path = self.out_dir / f"decorr_psf_{mess}"
                    decorr_zptimg_path = self.out_dir / f"decorr_zptimg_{mess}"
                    decorr_diff_path = self.out_dir / f"decorr_diff_{mess}"
                    for img, savepath, hdr in zip( [ sfftifier.PixA_DIFF_GPU, sfftifier.PixA_Ctarget_GPU, sfftifier.PSF_target_GPU ],
                                                   [ decorr_diff_path,        decorr_zptimg_path,         decorr_psf_path ],
                                                   [ sfftifier.hdr_target,    sfftifier.hdr_target,       None ] ):
                        with nvtx.annotate( "apply_decor", color=0xccccff ):
                            decorimg = sfftifier.apply_decorrelation( img )
                        with nvtx.annotate( "writefits", color=0xff8888 ):
                            fits.writeto( savepath, cp.asnumpy( decorimg ).T, header=hdr, overwrite=True )
                    sci_image.decorr_psf_path[ templ_image.image_name ]= decorr_psf_path
                    sci_image.decorr_zptimg_path[ templ_image.image_name ]= decorr_zptimg_path
                    sci_image.decorr_diff_path[ templ_image.image_name ]= decorr_diff_path


        if 'make_stamps' in steps:
            with nvtx.annotate( "make stamps", color=0xff8888 ):
                for templ_image in self.template_images:
                    stamp_name = stampmaker( self.ra, self.dec, templ_image.image_path, shape=np.array([100,100]),
                                             savedir=self.out_dir, savename=f"stamp_{templ_image.image_name}" )

                for sci_image in self.science_images:
                    for templ_image in self.template_images:
                        zptname = sci_image.decorr_zptimg_path[ templ_image.image_name ]
                        diffname = sci_image.decorr_diff_path[ templ_image.image_name ]
                        stamp_name = stampmaker( self.ra, self.dec, zptname, shape=np.array([100,100]),
                                                 savedir=self.out_dir, savename=f"stamp_{zptname.name}" )
                        sci_image.zpt_stamp_path[ templ_image.image_name ] = pathlib.Path( stamp_name )
                        stamp_name = stampmaker( self.ra, self.dec, diffname, shape=np.array([100,100]),
                                                 savedir=self.out_dir, savename=f"stamp_{diffname.name}" )
                        sci_image.diff_stamp_path[ templ_image.image_name ] = pathlib.Path( stamp_name )

        if 'make_lightcurve' in steps:
            with nvtx.annotate( "make_lightcurve", color=0xff8888 ):
                self.make_lightcurve()

# ======================================================================

def main():
    parser = argparse.ArgumentParser( 'phrosty pipeline' )
    parser.add_argument( '--oid', type=int, required=True, help="Object ID" )
    parser.add_argument( '-r', '--ra', type=float, required=True, help="Object RA" )
    parser.add_argument( '-d', '--dec', type=float, required=True, help="Object Dec" )
    parser.add_argument( '-b', '--band', type=str, required=True, help="Band: R062, Z087, Y106, J129, H158, F184, or K213" )
    parser.add_argument( '-t', '--template-images', type=str, required=True,
                         help="Path to file with, per line, ( path_to_image, pointing, sca )" )
    parser.add_argument( '-s', '--science-images', type=str, required=True,
                         help="Path to file with, per line, ( path_to_image, pointing, sca )" )
    parser.add_argument( '-n', '--ncpus', type=int, default=1, help="Number of CPUs for CPU multiprocessing steps" )
    parser.add_argument( '-v', '--verbose', action='store_true', default=False, help="Show debug log info" )
    parser.add_argument( '--out-dir', default="/dia_out_dir", help="Output dir, default /dia_out_dir" )
    parser.add_argument( '--ltcv-dir', default="/lc_out_dir", help="Output dir for lightcurves, default /lc_out_dir" )
    parser.add_argument( '--temp-dir', default="/phrosty_temp", help="Temporary working dir, default /phrosty_temp" )
    parser.add_argument( '--through-step', default='make_lightcurve',
                         help="Stop after this step; one of (see above)" )
    parser.add_argument( '--force-sky-subtract', action='store_true', default=False,
                         help='Redo sky subtraction even if the right file is sitting in the temp dir.' )

    args = parser.parse_args()

    science_images = []
    template_images = []
    for infile, imlist in zip( [ args.science_images, args.template_images ],
                               [ science_images, template_images ] ):
        with open( infile ) as ifp:
            for line in ifp:
                # WARNING: hardcoded header format
                if 'path pointing sca mjd' not in line:
                    img, point, sca, mjd = line.split()
                    imlist.append( ( pathlib.Path(img), int(point), int(sca), float(mjd) ) )

    # WARNING, hardcoded, fix this later
    galsim_config = pathlib.Path( os.getenv("SN_INFO_DIR" ) ) / "tds.yaml"

    pipeline = Pipeline( args.oid, args.ra, args.dec, args.band, science_images, template_images,
                         ncpus=args.ncpus, temp_dir=args.temp_dir, out_dir=args.out_dir, ltcv_dir=args.ltcv_dir,
                         galsim_config_file=galsim_config, force_sky_subtract=args.force_sky_subtract,
                         nuke_temp_dir=False, verbose=args.verbose )
    pipeline( args.through_step )

# ======================================================================
if __name__ == "__main__":
    main()
