import matplotlib.pyplot as plt
import numpy as np
from astropy.Table import Table

roman_bands = ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'W146', 'K213']

def photometric_repeatability(crossmatch_catalogs, savepath, stdev_endpoints=(20,23), bins=np.arange(10,35,0.5), title=None, figsize=figsize):
    """
    Argument crossmatch_catalogs should be a dictionary such that
    dict = {filter name: corresponding crossmatched catalog from ropeat.photometry.crossmatch_truth,
                    imported as an astropy table}
    """
    if not isinstance(crossmatch_catalogs, dict):
        raise ValueError: 
            print('Argument crossmatch_catalogs must be a dictionary.')
    elif any(crossmatch_catalogs.keys()) not in roman_bands:
        raise ValueError:
            print(f'All dictionary keys must be in {roman_bands}.')
    for b in crossmatch_catalogs.keys():
        if not isinstance(crossmatch_catalogs[b], astropy.table.table.Table):
            raise ValueError:
                print('All dictionary values must be type astropy.table.table.Table.')
    
    fig,ax = plt.subplots(len(crossmatch_catalogs),1, sharex=True, sharey=True, dpi=300, figsize=figsize)
    fig.subplots_adjust(hspace=0)
    
    for b in crossmatch_catalogs.keys():
        tab = crossmatch_catalogs[b]
        allmags = []
        allresidmean = []
        for row in tab:
        if isinstance(row[f'{b}_psf_mag'], str):
            mags_str = row[f'{b}_psf_mag'].split(',')
            make_float = lambda x : float(x)
            mags = np.array(list(map(make_float, mags_str)))
            nanmask = np.isnan(mags)
            mags = mags[~nanmask]
            allmags.append(mags)
            mean = np.mean(mags)
            resids = (mags - mean)
            resids_div_mean = resids/mean
            allresidmean.append(resids_div_mean)
            ax[i].plot([mean]*len(mags), resids_div_mean, color='k',
                    marker='.', linestyle='',markersize=2, alpha=0.5)
            if mean < stdev_endpoints[1] and mean > stdev_endpoints[0]:
                check_stdev.append(resids_div_mean)
                
        ax[i].axhline(0,color='k',linestyle='--',alpha=0.5)
        flattened_mags = np.array([item for sublist in allmags for item in sublist])
        flattened_resids_div_mean = np.array([item for sublist in allresidmean for item in sublist])
        flattened_check_stdev = np.array([item for sublist in check_stdev for item in sublist])
        scatter = 1.48*np.median(np.absolute(flattened_check_stdev))
        print(f'The scatter for {b} is {scatter}.')
        bins = bins
        digitized = np.digitize(flattened_mags, bins)
        binmeans = [np.mean(flattened_resids_div_mean[digitized==j]) for j in range(len(bins))]
        binstdevs = [np.std(flattened_resids_div_mean[digitized==j]) for j in range(len(bins))]

        bbox = ax[i].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width*fig.dpi, bbox.height*fig.dpi
        
        ax[i].annotate(r'$\sigma = $' + f'{scatter:.2E}', 
                    xy=(100,height-100), xycoords='axes pixels', ha='left', fontsize=22)
        ax[i].annotate(b, 
                    xy=(width-100,height-100), xycoords='axes pixels', ha='right', fontsize=22)
        ax[i].errorbar(bins,binmeans,yerr=binstdevs, color='red',
                  marker='o',markersize=5,capsize=5,linestyle='')
        ax[i].set_ylabel(rf'$\frac{{m_{{{b}, ij}} - \langle m_{{{b}, i}} \rangle}}{{\langle m_{{{b}, i}} \rangle}}$')
        ax[i].set_xlabel(r'$\langle m_{f, i} \rangle$')
        
    if title is not None:
        ax[0].set_title(title)
    plt.savefig(savepath, bbox_inches='tight', dpi=300)