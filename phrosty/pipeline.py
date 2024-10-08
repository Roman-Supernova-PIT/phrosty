import pathlib
import argparse
import logging

from sfft.SpaceSFFTCupyFlow import SpaceSFFT_CupyFlow

from phrosty.utils import set_logger
from phrosty.imagesubtraction import sky_subtract

class Image:
    def __init__( self, path, pointing, sca, pipeline ):
        self.pipeline = pipeline
        self.image_path = pathlib.Path( path )
        self.basename = self.image_path.name
        if self.basename[-3:] == '.gz':
            self.basename = self.basename[:-3]
        self.pointing = pointing
        self.sca = sca
        self.psf_path = None
        self.detect_mask_path = None
        self.skyrms = None
        self.skysub_path = None
        
    def run_sky_subtract( self ):
        self.skysub_path = self.pipeline.temp_dir / f"skysub_{self.basename}"
        self.detmask_path = self.pipeline.temp_dir / f"detmask_{self.basename}"
        self.skyrms = sky_subtract( self.image_path, self.skysub_path, self.detmask_path,
                                    temp_dir=self.pipeline.temp_dir, force=True )
        return ( self.skysub_path, self.detmask_path, self.skyrms )

    def save_sky_subtract_info( self, info ):
        self.skysub_path = info[0]
        self.detmask_path = info[1]
        self.skyrms_path = info[2]

class Pipeline:
    def __init__( self, object_id, ra, dec, band, science_images, template_images, ncpus=1,
                  temp_dir='/phrosty_temp', out_dir='/dia_out_dir', nuke_temp_dir=False, verbose=False ):
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
        self.object_id = object_id
        self.ra = ra
        self.dec = dec
        self.band = band
        self.science_images = [ Image( pps[0], pps[1], pps[2], self ) for pps in science_images ]
        self.template_images = [ Image( pps[0], pps[1], pps[2], self ) for pps in template_images ]
        self.ncpus = ncpus
        self.temp_dir = pathlib.Path(temp_dir )
        self.out_dir = pathlib.Path( out_dir )
        self.nuke_temp_dir = nuke_temp_dir

        self.logger = set_logger( 'phrosty', 'phrosty' )
        self.logger.setLevel( logging.DEBUG if verbose else logging.INFO )

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

        
    def align_and_pre_convolve(self, templ_image, sci_image ):
        """Align and pre convolve a single template/science pair.

        Parameters
        ----------
          sci_image: phrosty.Image
            The science (new) image.

          templ_image: phrosty.Image
            The template (ref) image that will be subtracted from sci_image.

        """

        # TODO : read these from Image
        # TODO : fix this whole function
        sci_image.skysub_path = self.temp_dir / f"skysub_{sci_image.image_path.name}"
        templ_image.skysub_path = self.temp_dir / f"skysub_{templ_image.image_path.name}"

        self.logger.debug(f'Path to sky-subtracted science image: \n {sci_image.skysub_path}')
        self.logger.debug(f'Path to sky-subtracted template image: \n {templ_image.skysub_path}')

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

    def __call__( self, through_step=None ):
        if through_step is None:
            through_step = 'make_lightcurve'

        steps = [ 'sky_subtract', 'get_psf', 'align_and_preconvolve', 'subtract', 'find_decorrelation',
                  'apply_decorrelation', 'make_stamps', 'make_lightcurve' ]
        stepdex = steps.index( through_step )
        if stepdex < 0:
            raise ValueError( f"Unknown step {through_step}" )
        steps = steps[:stepdex+1]

        if 'sky_subtract' in steps:
            self.sky_sub_all_images()

        if 'get_psf' in steps:
            raise NotImplementedError( "Gotta do this" )

        for templ_image in self.template_images:
            for sci_image in self.science_images:
                sfftifier = None
                if 'align_and_preconvolve' in steps:
                    sfftifier = self.align_and_pre_convolve( templ_image, sci_image )
            
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
    parser.add_argument( '--temp-dir', default="/phrosty_temp", help="Temporary working dir, default /phrosty_temp" )
    parser.add_argument( '--through-step', default='make_lightcurve',
                         help="Stop after this step; one of (see above)" )
    args = parser.parse_args()

    science_images = []
    template_images = []
    for infile, imlist in zip( [ args.science_images, args.template_images ],
                               [ science_images, template_images ] ):
        with open( infile ) as ifp:
            for line in ifp:
                img, point, sca = line.split()
                imlist.append( ( pathlib.Path(img), int(point), int(sca) ) )
    
    pipeline = Pipeline( args.oid, args.ra, args.dec, args.band, science_images, template_images,
                         ncpus=args.ncpus, temp_dir=args.temp_dir, out_dir=args.out_dir, nuke_temp_dir=False,
                         verbose=args.verbose )
    pipeline( args.through_step )

# ======================================================================
if __name__ == "__main__":
    main()
