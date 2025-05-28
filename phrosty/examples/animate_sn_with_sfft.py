# ruff: noqa
# Re-enable ruff when this file has been made to work again.

raise RuntimeError( "This whole file is currently broken, it refers to deprecated functions.  Needs updating." )

"""Saves a gif of the bright SN Ia with object ID 20202893. H158."""

from phrosty.imagesubtraction import sky_subtract, get_psf, imalign, crossconvolve, stampmaker
from phrosty.plotting import animate_stamps
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# NOTE TO SELF: UPDATE THIS SO THE INSTANCES OF THE OBJECT ARE LOCATED WITH THE CODE AS WELL.
# DO NOT WANT ANY PRE-MADE FILES HERE BECAUSE THAT IS USELESS AND UNHELPFUL.
tab = Table.read('/hpc/group/cosmology/lna18/roman_sim_imgs/Roman_Rubin_Sims_2024/'
                 '20202893/20202893_instances_nomjd.csv')
tab.sort('pointing')
tab = tab[tab['filter'] == 'H158']

# I already know that the reference image is the first row in the table (i.e., earliest image).
refvisit = 614
refsca = 4
ra = 8.037774
dec = -42.752337

# Process reference images.
ref_skysub_path = sky_subtract(band='H158', pointing=refvisit, sca=refsca)
ref_psf_path = get_psf(ra,dec,ref_skysub_path,ref_skysub_path,'H158',refvisit,refsca,'H158',refvisit,refsca)
print('Reference image sky subtracted and PSF retrieved.')

# Now, do the science images.
animate_panels = []
animate_raw_panels = []
animate_pointings = []
for row in tab[1:]:
    band = row['filter']
    pointing = row['pointing']
    sca = row['sca']

    print(band, pointing, sca)

    sci_skysub_path = sky_subtract(band=band,pointing=pointing,sca=sca)
    sci_imalign_path = imalign(ref_skysub_path,sci_skysub_path)
    sci_psf_path = get_psf(ra,dec,sci_imalign_path,sci_skysub_path,'H158',pointing,sca,'H158',refvisit,refsca)

    convolvedpaths = crossconvolve(sci_imalign_path, sci_psf_path, ref_skysub_path, ref_psf_path)
    stamp_paths = []
    for path in convolvedpaths:
        spath = stampmaker(ra,dec,path)
        stamp_paths.append(spath)

    # Also make stamp of sky subtracted original science image.
    rawstamp_savepath = f'/work/lna18/imsub_out/rawstamps/20202893_{band}_{pointing}_{sca}.fits'
    rawstamp = stampmaker(ra,dec,sci_imalign_path,savepath=rawstamp_savepath)

    scipath, refpath = stamp_paths
    diff, soln = sfft(scipath, refpath, sci_psf_path, ref_psf_path)
    decorr_path = decorr(scipath, refpath, sci_psf_path, ref_psf_path, diff, soln)
    finalimg = fits.open(decorr_path)[0].data
    finalstamp = finalimg[450:550,450:550]
    rawimg = fits.open(rawstamp_savepath)[0].data
    rawstamp = rawimg[450:550,450:550]

    animate_panels.append(finalstamp)
    animate_raw_panels.append(rawstamp)
    animate_pointings.append(pointing)

    zscale = ZScaleInterval()
    z1, z2 = zscale.get_limits(rawstamp)
    fig, ax = plt.subplots(figsize=(5,5))
    im = ax.imshow(rawstamp, cmap='Greys', vmin=z1, vmax=z2)
    ax.text(0.05,0.95,pointing,color='white',transform=ax.transAxes,va='top',ha='left')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')
    plt.savefig(f'figs/sfft/RAW_20202893_{band}_{pointing}_{sca}.png', dpi=300, bbox_inches='tight')
    plt.close()

# Make animation.
savepath = '/hpc/group/cosmology/lna18/roman_sim_imgs/Roman_Rubin_Sims_2024/20202893_SFFT.gif'
metadata = dict(artist='Lauren Aldoroty')
animate_stamps(animate_panels, savepath, metadata, labels=animate_pointings, staticlabel='H158')

savepath = '/hpc/group/cosmology/lna18/roman_sim_imgs/Roman_Rubin_Sims_2024/20202893_RAW.gif'
metadata = dict(artist='Lauren Aldoroty')
animate_stamps(animate_raw_panels, savepath, metadata, labels=animate_pointings, staticlabel='H158')
