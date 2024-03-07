import os
import numpy as np
import matplotlib.pyplot as plt

def standard_plot(filename, result, figsize=(12,10), xlims=None, max_eV=5):
    fig, sps = plt.subplots(nrows=3, ncols=3, figsize=figsize)
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    sps = sps.ravel()
    plt.suptitle(os.path.basename(filename))

    n_spikes = np.round(result['number_of_peaks']).astype(int)
    n_spikes_maxxed = np.clip(n_spikes, 0, 10)
    bar_x, counts = np.unique(n_spikes_maxxed, return_counts=True)
    ratios = counts/counts.sum()

    sp = sps[0]
    sp.set_title('Spike count (%i total)' % n_spikes.size)
    sp.set_xlabel('Number of spikes')
    sp.set_ylabel('Percentage')

    sp.bar(bar_x, ratios*100)
    sp.set_xticks(bar_x)
    mean_spikes = np.mean(n_spikes[n_spikes != 0])
    median_spikes = np.median(n_spikes[n_spikes != 0])
    sp.axvline(mean_spikes, label='Mean (excl. 0): %.1f' % mean_spikes, color='tab:orange')
    sp.axvline(median_spikes, label='Median (excl. 0): %.1f' % median_spikes, color='tab:green')
    xticklabels = ['%i' % x for x in bar_x]
    if n_spikes.max() > 10:
        xticklabels[-1] = xticklabels[-1]+'+'
    sp.set_xticklabels(xticklabels)
    sp.legend()

    all_spike_widths = []
    for list_ in result['all_spike_widths'].values():
        all_spike_widths.extend(list(list_))
    all_spike_widths0 = np.array(all_spike_widths)
    all_spike_widths = all_spike_widths0.clip(0, max_eV)

    sp = sps[1]
    sp.set_title('Spike widths (%i total, < %.1f eV)' % (len(all_spike_widths), max_eV))
    sp.set_xlabel('FWHM spike widths (eV)')
    sp.set_ylabel('Percentage')

    values, bin_edges = np.histogram(all_spike_widths, bins=10, density=True)
    sp.step(bin_edges[:-1], values*100, where='pre')
    mean = np.mean(all_spike_widths0)
    median = np.median(all_spike_widths0)
    sp.axvline(mean, label='Mean %.1f' % mean, color='tab:orange')
    sp.axvline(median, label='Median %.1f' % median, color='tab:green')
    sp.legend()

    data_min = np.mean(np.min(result['raw_data_intensity'],axis=0))
    data_max = (result['raw_data_intensity'] - data_min).max()
    #fit_max = np.max(result['fit_functions'])
    filtered_max = np.max(result['filtered_spectra'])

    avg = np.mean(result['filtered_spectra'], axis=0)

    for n_spectrum in range(7):
        sp = sps[n_spectrum+2]
        sp.set_title('Spectrum %i with %i spikes' % (n_spectrum, n_spikes[n_spectrum]))
        sp.set_xlabel('E (eV)')
        sp.set_ylabel('Intensity (arb. units)')

        _yy = result['raw_data_intensity'][n_spectrum]
        _yy_filtered = result['filtered_spectra'][n_spectrum]
        _yy_fit = result['fit_functions'][n_spectrum]
        ene = result['raw_data_energy']

        sp.plot(ene, _yy/data_max, label='Data')
        filter_plot = _yy_filtered/filtered_max
        sp.plot(ene, filter_plot, label='Filtered')
        _max = _yy_fit.max()
        if _max != 0:
            fit_plot = _yy_fit / _yy_fit.max()*filter_plot.max()
        else:
            fit_plot = _yy_fit
        sp.plot(ene, fit_plot, label='Fit')
        sp.plot(ene, avg/filtered_max, label='Avg', color='black', ls='--')

    sps[2].get_shared_x_axes().join(*sps[2:])
    sps[2].legend()
    if xlims is not None:
        sps[2].set_xlim(*xlims)

    return fig, sps

