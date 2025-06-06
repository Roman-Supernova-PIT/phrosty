# IMPORTS Standard:
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import animation
from matplotlib.patches import Rectangle
import numpy as np

# IMPORTS Astro:
from astropy.io import fits
from astropy.visualization import ZScaleInterval

# IMPORTS Internal:
from phrosty.utils import _build_filepath
import phrosty.plotaesthetics  # noqa: F401   (We need this for the commands run when it's imported.)

roman_bands = ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'W146', 'K213']


def showimage(path=None, band=None, pointing=None, sca=None,
              xycoords=None, box=False, boxsize=1000, data_ext=1,
              cmap='Greys', return_fig=False, **kwargs):
    """Quickly display an original RomanDESC image from either a filepath or a filter, pointing, and SCA ID."""

    path =_build_filepath(path=path, band=band, pointing=pointing, sca=sca, filetype='image')
    hdu = fits.open(path)
    img = hdu[data_ext].data
    z1, z2 = ZScaleInterval(n_samples=1000,
                            contrast=0.25,
                            max_reject=0.5,
                            min_npixels=5,
                            krej=2.5,
                            max_iterations=5,
                            ).get_limits(img)

    fig, ax = plt.subplots()
    imgax = ax.imshow(img, vmin=z1, vmax=z2, cmap=cmap, **kwargs)

    if xycoords is not None:
        ax.plot(xycoords[0], xycoords[1], markeredgecolor='r', marker='o', linestyle='',
                markersize=10, fillstyle='none')

        if box:
            rect = Rectangle((int(xycoords[0]-boxsize/2), int(xycoords[1]-boxsize/2)), width=boxsize,
                             height=boxsize, fill=None, edgecolor='r')
            ax.add_artist(rect)

    fig.colorbar(imgax, ax=ax)

    if return_fig:
        return fig
    else:
        plt.show()


class MidpointNormalize(mpl.colors.Normalize):
    """Normalise the colorbar."""
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def roman_sca_plot(data_array, sca_order, ptype='image', cmap='Greys', residual_plot=True, clabel=None, title=None,
                   vlims = None, return_fig=False, savefig=False, show_sca_id=False, savepath='roman_scas.png'):

    detector = plt.figure(figsize=(10, 6), dpi=300)
    nrows, ncols = 55, 91
    grid = detector.add_gridspec(nrows=nrows, ncols=ncols, figure=detector,
                                 width_ratios=[1]*ncols, height_ratios=[1]*nrows,
                                 hspace=0, wspace=0.1)
    row_begins = np.array([10, 3, 0, 0, 3, 10])
    row_ends = np.array([x+14 for x in row_begins])
    col_begins = np.arange(0, ncols, 14)
    col_ends = np.array([x+14 for x in col_begins])
    add_distance = [15, 16, 16]

    axs = []
    for row in add_distance:
        for i in range(len(row_begins)):
            ax = detector.add_subplot(grid[row_begins[i]:row_ends[i], col_begins[i]+1:col_ends[i]])
            ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)
            axs.append(ax)

        row_begins += row
        row_ends += row

    # Argument data_array should be an array of len(N SCAs) containing arrays:
    # fake_data = np.array([np.random.rand(14, 14)]*len(axs))
    if ptype=='image':
        if vlims is None:
            vmin, vmax = ZScaleInterval().get_limits(data_array.ravel())
        else:
            vmin, vmax = vlims

    sortidx = sca_order.argsort()
    sca_order = sca_order[sortidx]
    data_array = data_array[sortidx]
    imsim_sca_order = np.array([9, 6, 3, 12, 15, 18, 8, 5, 2, 11, 14, 17, 7, 4, 1, 10, 13, 16])-1

    for i, sca in enumerate(imsim_sca_order):
        if ptype=='image':
            if residual_plot:
                ends = np.nanmax(np.array([abs(vmin), abs(vmax)]))
                im = axs[i].imshow(data_array[sca], cmap=cmap,
                                       norm=MidpointNormalize(midpoint=0, vmin=-ends, vmax=ends))
            else:
                im = axs[i].imshow(data_array[sca], cmap=cmap, vmin=vmin, vmax=vmax)
        elif ptype=='scatter':
            axs[i].plot(data_array[sca][0], data_array[sca][1], marker='.', linestyle='', color='k', markersize=2)
            if residual_plot:
                axs[i].axhline(0, color='k', linestyle='--')

        if show_sca_id:
            axs[i].annotate(sca+1, xy=(0, 1), fontsize=12)

    if ptype=='image':
        cbar_ax = detector.add_subplot(grid[:, -4:-1])
        cbar = plt.colorbar(im, cax=cbar_ax)
        if clabel is not None:
            cbar.set_label(clabel, labelpad=20, fontsize=18, rotation=270)
    if title is not None:
        plt.suptitle(title, y=0.93, fontsize=18)
    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')

    if return_fig:
        return detector
    else:
        plt.show()


def classification_contours(data, model,
                            cmap='Set3',
                            plotxlims=None, plotylims=None,
                            contourxlims=(0, 1e5+5000), contourylims=(0, 9e6),
                            e_slices = [0.0, 0.2, 0.4, 0.6],
                            title=None,
                            savefig=False, savepath='figs/classification_model.png'):

    fig, ax = plt.subplots(2, 2, figsize=(12, 12), dpi=300, sharex=True, sharey=True)
    plt.subplots_adjust(hspace=0, wspace=0)

    x_min, x_max = contourxlims
    y_min, y_max = contourylims

    h = 1000
    xx, yy, zz = np.meshgrid(np.linspace(x_min, x_max, h),
                        np.linspace(y_min, y_max, h),
                            e_slices)

    predicted = model.predict(np.c_[xx.ravel(), yy.ravel(), zz.ravel()])
    predicted_grid = predicted.reshape(xx.shape)

    for i, a in enumerate(ax.flatten()):
        a.pcolormesh(xx[:, :, i], yy[:, :, i], predicted_grid[:, :, i], cmap=cmap, alpha=0.5)
        a.contour(xx[:, :, i], yy[:, :, i], predicted_grid[:, :, i])

        for t in ['train', 'test']:
            e_mask = np.logical_and(data[f'X_{t}']['ellipticity'] > e_slices[i] - 0.1,
                                    data[f'X_{t}']['ellipticity'] < e_slices[i] + 0.1)
            galaxies = data[f'X_{t}'][data[f'y_{t}'] == 0][e_mask]
            stars = data[f'X_{t}'][data[f'y_{t}'] == 1][e_mask]

            if t == 'train':
                markerfill = 'none'
            elif t == 'test':
                markerfill = 'full'

            a.plot(stars['peak'], stars['flux'], color='gold', fillstyle=markerfill, markeredgecolor='k',
                   linestyle='', marker='o', label=f'Stars, {t}', zorder=100, alpha=0.7)
            a.plot(galaxies['peak'], galaxies['flux'], color='green', fillstyle=markerfill, markeredgecolor='k',
                   linestyle='', marker='^', label=f'Galaxies, {t}', zorder=101, alpha=0.7)
            a.ticklabel_format(axis='both', style='sci', scilimits= (5, 5))
            a.tick_params(labelsize=24)
        a.annotate(f'$e = {{{e_slices[i]}}}$', xy=(.93, .9), xycoords='axes fraction', horizontalalignment='right')

    if plotxlims is not None:
        ax[0, 0].set_xlim(min(plotxlims), max(plotxlims))
    if plotylims is not None:
        ax[0, 0].set_ylim(min(plotylims), max(plotylims))
    fig.supxlabel('Maximum pixel value')
    fig.supylabel('Flux', x=.04)
    ax[0, 0].legend(loc='upper left')

    if title is not None:
        fig.suptitle(title)

    if savefig:
        plt.savefig(savepath, bbox_inches='tight', dpi=300)

    plt.show()


def animate_stamps(stamps, savepath, metadata=dict(), no_whitespace=True, labels=[], labelxy=(0.05, 0.95),
                   staticlabel=None, **kwargs):
    """_summary_

    :param stamps: Must be in chronological order.
    :type stamps: List of stamps from get_stamps or get_object_instances.
    :param savepath: _description_
    :type savepath: _type_
    :param metadata: _description_
    :type metadata: _type_
    """

    if no_whitespace:
        # NOTE: Your first axis (first indexing value) should return one stamp. e.g. stamps[0] is the first stamp.
        with_whitespace = np.invert(np.any((np.isnan(np.array(stamps))), axis=(1, 2)))
        idx_whitespace = np.where(with_whitespace)[0]
        stamps = np.array(stamps)[idx_whitespace]
        if len(labels) != 0:
            labels = np.array(labels)[idx_whitespace]

    fig, ax = plt.subplots(figsize=(5, 5))
    plt.tight_layout()
    ax.axis('off')
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    plt.xticks([])
    plt.yticks([])

    im = ax.imshow(stamps[0], animated=True)
    if staticlabel is not None:
        stxt = ax.text(0.95, 0.05, staticlabel, animated=True, color='white', transform=ax.transAxes,
                       va='bottom', ha='right', **kwargs)
    if len(labels) != 0:
        txt = ax.text(labelxy[0], labelxy[1], labels[0], animated=True, color='white', transform=ax.transAxes,
                      va='top', ha='left', **kwargs)

    def animate(i):
        im.set_array(stamps[i])
        if len(labels) != 0:
            txt.set_text(labels[i])

            return [im] + [txt] + [stxt]
        else:
            return [im]

    writer = animation.PillowWriter(metadata=metadata)
    anim = animation.FuncAnimation(fig, animate, interval=600, frames=len(stamps))
    anim.save(savepath, writer=writer)


# def photometric_repeatability(crossmatch_catalogs, savepath, stdev_endpoints=(20, 23),
#                               bins=np.arange(10, 35, 0.5), title=None, figsize=figsize):
#     """
#     Argument crossmatch_catalogs should be a dictionary such that
#     dict = {filter name: corresponding crossmatched catalog from ropeat.photometry.crossmatch_truth,
#                     imported as an astropy table}

#     THIS FUNCTION IS UNTESTED.

#     """
#     if not isinstance(crossmatch_catalogs, dict):
#         raise ValueError:
#             print('Argument crossmatch_catalogs must be a dictionary.')
#     elif any(crossmatch_catalogs.keys()) not in roman_bands:
#         raise ValueError:
#             print(f'All dictionary keys must be in {roman_bands}.')
#     for b in crossmatch_catalogs.keys():
#         if not isinstance(crossmatch_catalogs[b], astropy.table.table.Table):
#             raise ValueError:
#                 print('All dictionary values must be type astropy.table.table.Table.')

#     fig, ax = plt.subplots(len(crossmatch_catalogs), 1, sharex=True, sharey=True, dpi=300, figsize=figsize)
#     fig.subplots_adjust(hspace=0)

#     for b in crossmatch_catalogs.keys():
#         tab = crossmatch_catalogs[b]
#         allmags = []
#         allresidmean = []
#         for row in tab:
#         if isinstance(row[f'{b}_psf_mag'], str):
#             mags_str = row[f'{b}_psf_mag'].split(',')
#             make_float = lambda x : float(x)
#             mags = np.array(list(map(make_float, mags_str)))
#             nanmask = np.isnan(mags)
#             mags = mags[~nanmask]
#             allmags.append(mags)
#             mean = np.mean(mags)
#             resids = (mags - mean)
#             resids_div_mean = resids/mean
#             allresidmean.append(resids_div_mean)
#             ax[i].plot([mean]*len(mags), resids_div_mean, color='k',
#                     marker='.', linestyle='', markersize=2, alpha=0.5)
#             if mean < stdev_endpoints[1] and mean > stdev_endpoints[0]:
#                 check_stdev.append(resids_div_mean)

#         ax[i].axhline(0, color='k', linestyle='--', alpha=0.5)
#         flattened_mags = np.array([item for sublist in allmags for item in sublist])
#         flattened_resids_div_mean = np.array([item for sublist in allresidmean for item in sublist])
#         flattened_check_stdev = np.array([item for sublist in check_stdev for item in sublist])
#         scatter = 1.48*np.median(np.absolute(flattened_check_stdev))
#         print(f'The scatter for {b} is {scatter}.')
#         bins = bins
#         digitized = np.digitize(flattened_mags, bins)
#         binmeans = [np.mean(flattened_resids_div_mean[digitized==j]) for j in range(len(bins))]
#         binstdevs = [np.std(flattened_resids_div_mean[digitized==j]) for j in range(len(bins))]

#         bbox = ax[i].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#         width, height = bbox.width*fig.dpi, bbox.height*fig.dpi

#         ax[i].annotate(r'$\sigma = $' + f'{scatter:.2E}',
#                     xy=(100, height-100), xycoords='axes pixels', ha='left', fontsize=22)
#         ax[i].annotate(b,
#                     xy=(width-100, height-100), xycoords='axes pixels', ha='right', fontsize=22)
#         ax[i].errorbar(bins, binmeans, yerr=binstdevs, color='red',
#                   marker='o', markersize=5, capsize=5, linestyle='')
#         ax[i].set_ylabel(rf'$\frac{{m_{{{b}, ij}} - \langle m_{{{b}, i}} \rangle}}{{\langle m_{{{b}, i}} \rangle}}$')
#         ax[i].set_xlabel(r'$\langle m_{f, i} \rangle$')

#     if title is not None:
#         ax[0].set_title(title)
#     plt.savefig(savepath, bbox_inches='tight', dpi=300)
