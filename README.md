# phrosty

## NOTE: If you use this, do a git pull every day. I change this almost daily. 

"phrosty": PHotometry for ROman Simulations protoTYpe.

Basic package for working with the Roman-DESC simulations, associated with Aldoroty et al. 2024 in prep. 


Install by cloning this directory, navigating to that directory in your local terminal, and then using
```
pip install -e .
```
## Using the image subtraction module:
Make sure you update the output_files_rootdir variable in imagesubtraction.py, 
as well as rootdir (points to original RomanDESC images) in utils.py. Otherwise, 
you will get path errors.

DEPENDENCIES: 
You need to clone and install this directory:
https://github.com/thomasvrussell/sfft

You also need GalSim, but you can't pip install it at this time because the features
you need are not yet pushed to pip. So, clone and install. 
https://github.com/GalSim-developers/GalSim

Lastly, you need to clone and install the v2.0 branch of roman_imsim.
https://github.com/matroxel/roman_imsim

```
from phrosty.imagesubtraction import *
from phrosty.plotting import showimage

# Which images are we using? 
band = 'H158'
refvisit = 19138
refsca = 14
scivisit = 36445
scisca = 1
ra = 8.037774
dec = -42.752337

# Display a raw image. 
showimage(band='H158',pointing=scivisit,sca=scisca)

# Do sky subtraction on both the science and reference images, then display one.
sci_skysub_path = sky_subtract(band='H158',pointing=scivisit,sca=scisca)
ref_skysub_path = sky_subtract(band='H158',pointing=refvisit,sca=refsca)
showimage(path=sci_skysub_path,data_ext=0)

# Align the science image to the reference image and display. 
sci_imalign_path = imalign(ref_skysub_path,sci_skysub_path)
showimage(path=sci_imalign_path,data_ext=0)

# Retrieve PSFs at the chosen RA, dec location for each image, and rotate to the alignment. 
sci_psf_path = get_psf(ra,dec,sci_imalign_path,sci_skysub_path,'H158',
                        scivisit,scisca,'H158',refvisit,refsca)
showimage(path=sci_psf_path,data_ext=0)
ref_psf_path = get_psf(ra,dec,ref_skysub_path,ref_skysub_path,'H158',
                        refvisit,refsca,'H158',refvisit,refsca)
showimage(path=ref_psf_path,data_ext=0)

# Cross-convolve the PSFs and images.
convolvedpaths = crossconvolve(sci_imalign_path, sci_psf_path, ref_skysub_path, ref_psf_path)
for path in convolvedpaths:
    showimage(path,data_ext=0)

# Cut out 1000x1000 stamps centered at the chosen RA, dec. 
spaths = []
for path in convolvedpaths:
    spath = stampmaker(ra, dec, path)
    spaths.append(spath)
    showimage(spath,data_ext=0)

# Do SFFT subtraction. 
scipath, refpath = spaths
diff, soln = difference(scipath, refpath, sci_psf_path, ref_psf_path)
showimage(diff, data_ext=0)

# Make your decorrelation kernel.
dcker_path = decorr_kernel(scipath, refpath, sci_psf_path, ref_psf_path, diff, soln)

# Decorrelate the difference image. 
decorr_diff_path = decorr_img(diff, dcker_path, imgtype='difference')
showimage(decorr_path, data_ext=0)

# Decorrelate the cross-convolved science image so you can get the zeropoint from this correctly.
decorr_zpt_path = decorr_img(spaths[0], dcker_path, imgtype='science')

# Calculate your final PSF to be used on your difference image and decorrelated cross-convolved science image. 
psfpath = calc_psf(spaths[0], spaths[1], sci_psf_path, ref_psf_path, dcker_path)
```

## Old readme contents:
Currently, you will also need GalSim (clone from git repo), roman_imsim (also clone from git repo), and SFFT (clone from my forked repo **FOR NOW. A pull request on the original repo will be merged soon.)

You will likely need to change the "rootdir" variable in utils.py, and output_files_rootdir in imagesubtraction.py. 

This package is modular. Modules are compatible with each other, and they interface
easily, but you do not need to use all the modules. It contains:
```
sourcedetection.py
objectidentification.py
photometry.py
utils.py
plotaesthetics.py
plotting.py
imagesubtraction.py
```

You'll need:
- Your science image
- The truth file corresponding to that science image
- WCS information 

Photometry usage example:
(THIS IS OLD! PROBABLY DOESN'T WORK!)
```
# Standard imports
import numpy as np
import pandas as pd

# Astro imports
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

# This package
from phrosty.sourcedetection import detect_sources, plot_sources, catalog_matching
from phrosty.photometry import ap_phot

imgpath = '/path/to/science/image.fits.gz'
truthpath = '/path/to/truth/file.txt'

hdu = fits.open(imgpath)
wcs = WCS(hdu[1].header)
img = hdu[1].data

objects = detect_sources(img)
# plot_sources(img, objects)

# Note: column names in your truth object table must be distinct from the names of the output columns,
# and the suffix must be '_truth'. e.g. don't put 'x' in your truth object table, put 'x_truth'.
truth_colnames = ['object_id', 'ra_truth', 'dec_truth', 'x_truth', 'y_truth', 'realized_flux', 'flux_truth', 'mag_truth']
truthread = pd.read_csv(truthpath, skipinitialspace=True, names=truth_colnames, sep=' ', comment='#')
cat = Table(Table.from_pandas(truthread), masked=True)

matched_catalog = catalog_matching(objects, cat, wcs)
detected_coords = matched_catalog[matched_catalog['detection'] == 1]

apphot = ap_phot(img, detected_coords)
```