# ropeat
Basic package for Roman photometry, associated with Aldoroty et al. 2024 in prep. 

Current state of functionality: Relies on truth files. i.e., star/galaxy separation not implemented, source finder has initial parameters, magnitudes are zeropointed to the truth. 

To-do: 
- Implement star/galaxy separation
- Ensure accurate source detection
- Zero point for real
- Replace current DAOStarFinder ePSF stamp exclusion method with one based on nearness to other coordinates. 

Usage example:
```
import ropeat.photometry as rp
from glob import glob 

# Note: truth file from https://roman.ipac.caltech.edu/data/sims/sn_image_sims/galaxia_akari.fits.gz
# But it'll work as long as you have a fits file with columns 'ra', 'dec', 'Y106', 'J129', 'H158', and 'F184' in radians and magnitudes. 
# But don't attach an astropy unit to them. It converts the radians to degrees.  
truth_filepath = 'imgs/truth/galaxia_akari.fits'
science_filepath = glob('imgs/rotate_update_Y106_*_1.fits')

# Define truth object and pull out the coordinates as astropy SkyCoord objects.
# View the table you imported.
tr = rp.truth(truth_filepath)
tc = tr.truthcoords()
tt = tr.table()

# Define a science image object.
for file in science_filepath:
    sci = rp.scienceimg(file,tc,band='Y106',pointing='141',chip='1')
    sci.plot_truth() # Plot the locations of the truth coordinates on your image.
    ap = sci.ap_phot(3.0) # Do aperture photometry on the object.
    psf = sci.psf_phot('epsf',ap,plot_epsf=True) # Construct a
                                                 # PSF from field   
                                                 # stars and do PSF photometry. 
                                                 # Also, plot the PSF. This 
                                                 # returns a table as well. 
    savepath = f'tables/test/test_{sci.band}_{sci.pointing}_{sci.chip}_tab.csv'
    sci.save_results(savepath,truth_table=tt,truth_zpt=True,overwrite=True)

# Now, crossmatch the science image results to the truth catalog.
# This returns a table where all measurements for one star are in a list
# all in the same cell.
results_filepath = glob('tables/test/*)
cm = rp.crossmatch_truth(truth_filepath,results_filepath)
```