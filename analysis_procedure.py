import inspect
import os
import numpy as np
from uncertainties import unumpy

#from qtpy.QtCore import QObject, Signal

from .src.spectrum_dataset import GUISpectrumDataset
_pymodule = os.path.basename(__file__)

#Special values for progress_update Signal
_PROGRESS_BAR_THREAD_START = 0
_PROGRESS_BAR_THREAD_ERROR = 3
_PROGRESS_BAR_THREAD_END = 100

def __LINE__():
    """Macro to return the current line number.

    The current line number within the file is used when
    reporting messages to the message logging window.

    Returns:
        int: Current line number.
    """
    return inspect.currentframe().f_back.f_lineno


class AnalysisProcedure():
    #trigger_abort = Signal()
    #trigger_analysis_message = Signal(str, str)

    def __init__(self, parent=None):
        #super(AnalysisProcedure, self).__init__(parent)
        self.parent = parent
        self.abort = False
        #self.trigger_abort.connect(self.abort_update)
        self.dataset = None
        self.parallel = True
        self.pool = None

    def __del__(self):
        if self.pool is not None:
            self.pool.terminate()
            self.pool.close()

    def abort_update(self):
        self.abort = True


    def PerformMeas(self, input_parameters):

        print(input_parameters)

        if self.parent is not None:
            self.parent.trigger_progressbar.emit(_PROGRESS_BAR_THREAD_START)

        self.dataset = GUISpectrumDataset(self, input_parameters, self.parallel)

        if not self.dataset.open_dataset():
            print('Cannot open dataset')
            return

        self.dataset.get_general_information()

        if input_parameters['spike-width']:
            if not self.dataset.analyse():
                return

        if input_parameters['auto-correlation']:
            self.dataset.second_order_correlation()

        if self.abort:
            self.abort = False
            return
        if self.parent is not None:
            self.parent.trigger_progressbar.emit(_PROGRESS_BAR_THREAD_END)

        print(self.dataset.tauto)
        results_dict = {'input_parameters': input_parameters}
        results_dict['number_of_spectra'] = self.dataset.spectrum_number
        results_dict['average_spectrum_width'] = self.dataset.bandwidth
        results_dict['average_photon_energy'] = [self.dataset.average_photon_energy, self.dataset.std_photon_energy]
        results_dict['raw_data_energy'] = self.dataset.photE
        results_dict['raw_data_intensity'] = self.dataset.intense
        results_dict['fit_functions'] = self.dataset.fit_functions
        results_dict['filtered_spectra'] = self.dataset.filtered_spectra
        results_dict['tmin'] = self.dataset.tmin

        if self.dataset.tmin_avg is not None:
            results_dict['average_spike_width'] = [self.dataset.overall_average_FWHM, self.dataset.FWHM_std]
            results_dict['minimal_pulse_duration'] = [self.dataset.tmin_avg, self.dataset.tmin_std]
            results_dict['maximal_pulse_duration'] = [self.dataset.tmax_avg, self.dataset.tmax_std]
            results_dict['spike_widths'] = [self.dataset.FWHM_per_spike, self.dataset.FWHM_std_per_spike]
            results_dict['average_spike_number'] = self.dataset.average_peak_number
            results_dict['peak_histogram'] = self.dataset.histogram
        else:
            results_dict['average_spike_width'] = [0, 0]
            results_dict['minimal_pulse_duration'] = [0, 0]
            results_dict['maximal_pulse_duration'] = [0, 0]
            results_dict['spike_widths'] = [np.zeros(5), np.zeros(5)]
            results_dict['average_spike_number'] = 0
            results_dict['peak_histogram'] = np.zeros(5)

        if self.dataset.noisy_spectra is None:
            results_dict['noisy_spectra'] = 0
        else:
            results_dict['noisy_spectra'] = self.dataset.noisy_spectra

        results_dict['pulse_energy'] = self.dataset.average_energy
        results_dict['relative_energy'] = self.dataset.all_energy
        results_dict['intensity_fluctuations'] = self.dataset.intensity_fluctuations
        results_dict['average_spectrum'] = [self.dataset.average_spectra, self.dataset.fit.gausscurve]
        results_dict['photon_energy'] = self.dataset.photE

        if self.dataset.tauto is not None:
            print(unumpy.nominal_values(self.dataset.tauto), unumpy.nominal_values(self.dataset.resolution))
            results_dict['tauto'] = unumpy.nominal_values(self.dataset.tauto)
            results_dict['resolution'] = unumpy.nominal_values(self.dataset.resolution)
        else:
            results_dict['tauto'] = 0
            results_dict['resolution'] = 0

        results_dict["Figure1"] = self.dataset.get_average_figure()
        results_dict["Figure2"] = self.dataset.get_histogram()
        results_dict['Figure3'] = self.dataset.get_average_figure()
        results_dict['Figure4'] = self.dataset.get_histogram()
        results_dict['Figure5'] = self.dataset.get_auto_correlation_plot()

        results_dict['examples'] = self.dataset.examples
        results_dict['number_of_peaks'] = self.dataset.number_of_peaks
        results_dict['all_spike_widths'] = self.dataset.all_spike_widths
        results_dict['all_spike_amplitudes'] = self.dataset.all_spike_amplitudes
        results_dict['all_spike_centers'] = self.dataset.all_spike_centers

        return results_dict

    def GetFigures(self):
        figures = self.dataset.get_figure_list()
        self.dataset = None
        return figures

