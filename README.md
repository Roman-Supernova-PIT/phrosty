# ropeat
Basic package for Roman photometry, associated with Aldoroty et al. 2024 in prep. 

Current state of functionality: Relies on truth files. i.e., star/galaxy separation not implemented, source finder has initial parameters, magnitudes are zeropointed to the truth. 

To-do: 
- Implement star/galaxy separation
- Ensure accurate source detection
- Zero point for real

Usage example:
```
import ropeat.photometry as rp
# Note: truth file from https://roman.ipac.caltech.edu/data/sims/sn_image_sims/galaxia_akari.fits.gz
truth_filepath = 'imgs/truth/galaxia_akari.fits'
science_filepath = 'imgs/rotate_update_Y106_141_1.fits'

# Define truth object and pull out the coordinates as astropy SkyCoord objects.
# View the table you imported.
tr = rp.truth(truth_filepath)
tc = tr.truthcoords()
print(tr.table())

# Define a science image object.
sci = rp.scienceimg(science_filepath,tc,band='Y106',pointing='141',chip='1')
sci.plot_truth()
ap = sci.ap_phot(3.0) # Do aperture photometry on the object.
sci.psf_phot('epsf',ap,plot_epsf=True) # Construct a PSF from field stars and
                                       # do PSF photometry. Also, plot the 
                                       # PSF. This returns a table as well. 
```