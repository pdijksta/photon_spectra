import os
from multiprocessing import Pool
import numpy as np
from scipy.signal import butter, filtfilt, find_peaks

from PassiveWFMeasurement import h5_storage
from PassiveWFMeasurement import myplotstyle as ms

import spectrum
from photon_spectra.src.single_spectra import single_spectra
import plot_results

ms.closeall()

tool_analysis = False
butter_deg = 5

textbbox = {'boxstyle': 'square', 'alpha': 0.75, 'facecolor': 'white', 'edgecolor': 'gray'}

params = {
        'spike-width': 1,
        'auto-correlation': 0,
        'parameter optimization': 0,
        'pulse duration correction': 0,


        'height': 4,
        'prominence': 8,
        'snr': 1,
        'sbr': 1,
        'intensity_thresh': 0.1,
        'frequency_cutoff': 0.04,
        'background_prominence': 0,
        'cutoff_prominence': 5,
        'cutoff_height': 20.0,
        }

fig = ms.figure('Spectrum analysis')
subplot = ms.subplot_factory(5, 5)
sp_ctr = 1

def remove_duplicates(intensity0):
    _, indices = np.unique(intensity0, axis=0, return_index=True)
    intensity = intensity0[np.sort(indices)]
    intensity = intensity[np.sum(intensity, axis=1) != 0]
    return intensity

def butterfilter(single_intensity, lowpasscutoff):
    b, a = butter(butter_deg, lowpasscutoff, 'low')
    return filtfilt(b, a, single_intensity)

def load_sven_tool(filename):
    data = h5_storage.loadH5Recursive(filename)
    photon_energy = data['scan_1']['data']['SATOP31-PMOS132-2D']['SPECTRUM_X']
    intensity0 = data['scan_1']['data']['SATOP31-PMOS132-2D']['SPECTRUM_Y']
    return photon_energy, intensity0


dfile = '/sf/data/measurements/2024/06/14/SpectralAnalysis_2024_06_14_12_16_41_230992.h5'
photon_energy, intensity0 = load_sven_tool(dfile)

dfile_blank = '/sf/data/measurements/2024/06/14/SpectralAnalysis_2024_06_14_11_50_55_270269.h5'
photon_energy, intensity_blank = load_sven_tool(dfile_blank)

mean_blank = np.mean(intensity_blank, axis=0)

intensity = remove_duplicates(intensity0) - mean_blank

filtered_intensity = np.array([butterfilter(x, params['frequency_cutoff']) for x in intensity])

class SingleSpectrum(single_spectra):
    fact = 1

    def __init__(self, photon_energy, raw_intensity, filtered_intensity, background):
        self.photE = photon_energy
        self.raw_intensity = raw_intensity
        self.filtered_intensity = self.shiftraw = self.intensity = filtered_intensity
        self.average_intensity = filtered_intensity.mean()
        self.background = background
        self.max_intensity = filtered_intensity.max()
        self.indivtime = np.zeros(4)

    def find_peaks(self, params, maxint):
        self.minheight = params['height']*np.max(filtered_intensity)/100
        self.minprominence = params['prominence']*maxint/100
        self.xpeaks = find_peaks(self.filtered_intensity, height=self.minheight, prominence=self.minprominence)[0]
        self.ypeaks = self.filtered_intensity[self.xpeaks]
        self.n = len(self.xpeaks)

all_single_spectra = []

#for index, single_intensity in enumerate(intensity):
def analyze_single(single_data):
    single_spectrum = SingleSpectrum(photon_energy, single_data[0], single_data[1], background=0)
    single_spectrum.find_peaks(params, filtered_intensity.max())
    single_spectrum.first_estimate(params['cutoff_prominence'], params['cutoff_height'], globalminima=False)
    return single_spectrum

try:
    pool = Pool(8)
    all_single_spectra = pool.map(analyze_single, list(zip(intensity[:20], filtered_intensity[:20])))
finally:
    pool.close()
    pool.join()

nspectra = len(all_single_spectra)

new_result_dict = {
        'input_parameters': params,
        'number_of_spectra': nspectra,
        'raw_data_intensity': intensity,
        'raw_data_energy': photon_energy,
        'filtered_spectra': filtered_intensity,
        'fit_functions': np.empty([nspectra, len(photon_energy)]),
        'number_of_peaks': np.empty(nspectra),
        'all_spike_widths': {},
        }

for ctr, single_spectrum in enumerate(all_single_spectra):
    new_result_dict['number_of_peaks'][ctr] = single_spectrum.n
    new_result_dict['all_spike_widths'][str(ctr)] = single_spectrum.all_fwhm
    new_result_dict['fit_functions'][ctr] = single_spectrum.total_fit

for index, single_spectrum in enumerate(all_single_spectra[:25]):
    sp = subplot(sp_ctr, title='Spectrum %i' % index, xlabel='Photon energy (eV)', ylabel='Intensity (arb. units)')
    sp_ctr += 1

    sp.plot(photon_energy, single_spectrum.raw_intensity)
    sp.plot(photon_energy, single_spectrum.filtered_intensity)
    sp.plot(photon_energy, single_spectrum.total_fit)
    for peak in single_spectrum.xpeaks:
        sp.axvline(photon_energy[peak], color='black', ls='--')
        textstr = 'Min height %.1f\nProminence %.1f' % (single_spectrum.minheight, single_spectrum.minprominence)
        sp.text(0.05, 0.05, textstr, transform=sp.transAxes, verticalalignment='bottom', bbox=textbbox)

try:
    plot_results.standard_plot('Test', new_result_dict, norm_plot=False, remove_min_plot=False)
except Exception as e:
    print(e)

if tool_analysis:
    _, result_dict = spectrum.analyze_spectrum(dfile, params)
    plot_results.standard_plot(os.path.basename(dfile), result_dict)

ms.show()

