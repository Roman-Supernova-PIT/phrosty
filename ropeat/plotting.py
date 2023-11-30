import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
import matplotlib as mpl
from matplotlib.colors import Normalize
import numpy as np
from astropy.table import Table
from . import plotaesthetics

roman_bands = ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'W146', 'K213']

# def photometric_repeatability(crossmatch_catalogs, savepath, stdev_endpoints=(20,23), bins=np.arange(10,35,0.5), title=None, figsize=figsize):
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
    
#     fig,ax = plt.subplots(len(crossmatch_catalogs),1, sharex=True, sharey=True, dpi=300, figsize=figsize)
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
#                     marker='.', linestyle='',markersize=2, alpha=0.5)
#             if mean < stdev_endpoints[1] and mean > stdev_endpoints[0]:
#                 check_stdev.append(resids_div_mean)
                
#         ax[i].axhline(0,color='k',linestyle='--',alpha=0.5)
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
#                     xy=(100,height-100), xycoords='axes pixels', ha='left', fontsize=22)
#         ax[i].annotate(b, 
#                     xy=(width-100,height-100), xycoords='axes pixels', ha='right', fontsize=22)
#         ax[i].errorbar(bins,binmeans,yerr=binstdevs, color='red',
#                   marker='o',markersize=5,capsize=5,linestyle='')
#         ax[i].set_ylabel(rf'$\frac{{m_{{{b}, ij}} - \langle m_{{{b}, i}} \rangle}}{{\langle m_{{{b}, i}} \rangle}}$')
#         ax[i].set_xlabel(r'$\langle m_{f, i} \rangle$')
        
#     if title is not None:
#         ax[0].set_title(title)
#     plt.savefig(savepath, bbox_inches='tight', dpi=300)
    
    
class MidpointNormalize(mpl.colors.Normalize):
    """Normalise the colorbar."""
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def roman_sca_plot(data_array,sca_order,ptype='image',residual_plot=True,clabel=None,title=None,
                   savefig=False,show_sca_id=False,savepath='roman_scas.png'):
    
    detector = plt.figure(figsize=(10,6),dpi=300)
    nrows, ncols = 55,91
    grid = detector.add_gridspec(nrows=nrows,ncols=ncols,figure=detector, 
                                 width_ratios=[1]*ncols, height_ratios=[1]*nrows,
                                 hspace=0,wspace=0.1)
    row_begins = np.array([10,3,0,0,3,10])
    row_ends = np.array([x+14 for x in row_begins])
    col_begins = np.arange(0,ncols,14)
    col_ends = np.array([x+14 for x in col_begins])
    add_distance = [15,16,16]

    axs = []
    for row in add_distance:
        for i in range(len(row_begins)):
            ax = detector.add_subplot(grid[row_begins[i]:row_ends[i],col_begins[i]+1:col_ends[i]])
            ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)
            axs.append(ax)

        row_begins += row
        row_ends += row

    # Argument data_array should be an array of len(N SCAs) containing arrays:
    # fake_data = np.array([np.random.rand(14,14)]*len(axs))
    if ptype=='image':
        vmin = np.nanmin(data_array.ravel())
        vmax = np.nanmax(data_array.ravel())

    sortidx = sca_order.argsort()
    sca_order = sca_order[sortidx]
    data_array = data_array[sortidx]
    imsim_sca_order = np.array([9,6,3,12,15,18,8,5,2,11,14,17,7,4,1,10,13,16])-1

    for i, sca in enumerate(imsim_sca_order):
        if ptype=='image':
            if residual_plot:
                ends = np.nanmax(np.array([abs(vmin),abs(vmax)]))
                im = axs[i].imshow(data_array[sca], cmap='seismic',
                                       norm=MidpointNormalize(midpoint=0,vmin=-ends,vmax=ends))
            else:
                im = axs[i].imshow(data_array[sca], cmap='plasma', vmin=vmin,vmax=vmax)
        elif ptype=='scatter':
            axs[i].plot(data_array[sca][0],data_array[sca][1],marker='.',linestyle='',color='k',markersize=2)
            if residual_plot:
                axs[i].axhline(0,color='k',linestyle='--')
            
        if show_sca_id:
            axs[i].annotate(sca+1, xy=(0,2), fontsize=12)
    
    if ptype=='image':
        cbar_ax = detector.add_subplot(grid[:,-4:-1])
        cbar = plt.colorbar(im, cax=cbar_ax)
        if clabel is not None:
            cbar.set_label(clabel, labelpad=20, fontsize=18, rotation=270)
    if title is not None:
        plt.suptitle(title, y=0.93, fontsize=18)
    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')

    plt.show()