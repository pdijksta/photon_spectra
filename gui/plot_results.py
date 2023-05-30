import os
import numpy as np
import matplotlib.pyplot as plt

def standard_plot(filename, result, figsize=(12,10), xlims=None):
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
    sp.axvline(np.mean(n_spikes[n_spikes != 0]), label='Mean (excl. 0)', color='tab:orange')
    sp.axvline(np.median(n_spikes[n_spikes != 0]), label='Median (excl. 0)', color='tab:green')
    xticklabels = ['%i' % x for x in bar_x]
    if bar_x[-1] > 10:
        xticklabels[-1] = xticklabels[-1]+'+'
    sp.set_xticklabels(xticklabels)
    sp.legend()

    all_spike_widths = []
    for list_ in result['all_spike_widths'].values():
        all_spike_widths.extend(list(list_))
    all_spike_widths0 = np.array(all_spike_widths)
    all_spike_widths = all_spike_widths0.clip(0, 15)

    sp = sps[1]
    sp.set_title('Spike widths (%i total, < 15 eV)' % len(all_spike_widths))
    sp.set_xlabel('FWHM spike widths (eV)')
    sp.set_ylabel('Percentage')

    values, bin_edges = np.histogram(all_spike_widths, bins=10, density=True)
    sp.step(bin_edges[:-1], values*100, where='pre')
    sp.axvline(np.mean(all_spike_widths0), label='Mean', color='tab:orange')
    sp.axvline(np.median(all_spike_widths0), label='Median', color='tab:green')
    sp.legend()

    data_min = np.mean(np.min(result['raw_data_intensity'],axis=0))
    data_max = (result['raw_data_intensity'] - data_min).max()
    fit_max = np.max(result['fit_functions'])
    filtered_max = np.max(result['filtered_spectra'])

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
        sp.plot(ene, _yy_filtered/filtered_max, label='Filtered')
        sp.plot(ene, _yy_fit/fit_max, label='Fit')

    sps[2].get_shared_x_axes().join(*sps[2:])
    sps[2].legend()
    if xlims is not None:
        sps[2].set_xlim(*xlims)

    return fig, sps

