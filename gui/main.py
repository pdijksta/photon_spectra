import os
import sys
import socket
import time
import base64

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import PyQt5
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QMainWindow, QApplication

from h5_storage import saveH5Recursive
import logbook

if '../../' not in sys.path:
    sys.path.append('../../')

import spectrum

if __name__ == '__main__' and (not os.path.isfile('./gui.py') or os.path.getmtime('./gui.ui') > os.path.getmtime('./gui.py')):
    cmd = 'bash ./ui2py.sh'
    print(cmd)
    os.system(cmd)

from gui import Ui_MainWindow

default_input_parameters = {
    'height': 20,
    'prominence': 5,
    'spike-width': 1,
    'auto-correlation': 0,
    'parameter optimization': 0,
    'pulse duration correction': 0,
    'snr': 1,
    'sbr': 0.1,
    'intensity_thresh': 0.25,
    'roughness': 0.5,
    'frequency_cutoff': 0.05,
    'background_prominence': 6,
    'cutoff_prominence': 8,
    'cutoff_height': 25,
    }

class Main(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.SelectFile.clicked.connect(self.select_file(self.ui.Filename))
        self.ui.DoAnalysis.clicked.connect(self.do_analysis)
        self.ui.DoLogbook.clicked.connect(self.do_logbook)

        if 'xfelbkr' in socket.gethostname():
            self.ui.Filename.setText('/Users/xfeloper/user/pySpectrometer/SASE2/20230413-19_10_59_waterflow.npz')
            self.default_path = '/Users/xfeloper/user/pySpectrometer/'
        else:
            self.ui.Filename.setText('./test_data/20230413-19_10_59_waterflow.npz')
            self.default_path = './test_data/'

        self.result_dict = None

    def do_analysis(self):
        self.result_dict = self.fig_savename = self.save_filename = None

        self.filename = filename = self.ui.Filename.text().strip()
        if not os.path.isfile(filename):
            raise ValueError('%s does not exist' % filename)
        parameters = self.get_parameters()
        print('Start analysis')
        time0 = time.time()
        _, result_dict = spectrum.analyze_spectrum(filename, parameters)
        time1 = time.time()
        print('End analysis after %.0f s' % (time1-time0))
        save_filename = os.path.abspath('./analyzed_data/'+os.path.basename(filename).replace('.npz', '_analyzed.h5'))
        saveH5Recursive(save_filename, result_dict)
        print('Saved %s' % save_filename)

        self.fig = self.do_plot(result_dict, filename)
        self.save_fig()
        self.result_dict = result_dict
        self.save_filename = save_filename

    def save_fig(self):
        self.fig_savename = os.path.abspath('./analyzed_data/'+os.path.basename(self.filename).replace('.npz', '_analyzed.png'))
        self.fig.savefig(self.fig_savename)
        print('Saved %s' % self.fig_savename)

    def do_plot(self, result, filename):
        self.fig, sps = plt.subplots(nrows=3, ncols=3, figsize=(12,10))
        self.fig.subplots_adjust(hspace=0.5, wspace=0.5)
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
        xticklabels = ['%i' % x for x in bar_x]
        xticklabels[-1] = xticklabels[-1]+'+'
        sp.set_xticklabels(xticklabels)

        all_spike_widths = []
        for list_ in result['all_spike_widths'].values():
            all_spike_widths.extend(list(list_))
        all_spike_widths = np.array(all_spike_widths)
        all_spike_widths = all_spike_widths[all_spike_widths < 15]

        sp = sps[1]
        sp.set_title('Spike widths (%i total, < 15 fs)' % len(all_spike_widths))
        sp.set_xlabel('FWHM spike widths (fs)')
        sp.set_ylabel('Percentage')

        values, bin_edges = np.histogram(all_spike_widths, bins=10, density=True)
        bin_edges2 = (bin_edges[:-1] + bin_edges[1:])/2
        sp.bar(bin_edges2, values*100)

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

            sp.plot(ene, _yy/data_max)
            sp.plot(ene, _yy_filtered/filtered_max)
            sp.plot(ene, _yy_fit/fit_max)

        sps[2].get_shared_x_axes().join(*sps[2:])

        plt.show(block=False)
        return self.fig

    def do_logbook(self):
        if self.result_dict is None:
            print('No result to log.')
            return
        comment = 'File: %s\nAnalyzed File: %s\n\n' % (self.filename, self.save_filename)
        comment += parameters_to_text(self.result_dict['input_parameters'])
        self.save_fig()
        with open(self.fig_savename, 'rb') as f:
            image = base64.b64encode(f.read()).decode('ascii')
        succeeded = logbook.send_to_desy_elog('SASE spectrum fit', os.path.basename(self.filename), 'INFO', comment, 'xfellog', image)
        if succeeded:
            print('Logbook save successful.')
        else:
            print('Logbook save unsuccessful.')

    def select_file(self, widget):
        def f():
            QFileDialog = PyQt5.QtWidgets.QFileDialog
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            filename, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", self.default_path, "Npz files (*.npz);;All Files (*)", options=options)
            if filename:
                widget.setText(filename)
        return f

    def get_parameters(self):
        parameters = default_input_parameters.copy()
        parameters.update({
                'height': self.ui.ParHeight.value(),
                'prominence': self.ui.ParProminence.value(),
                'snr': self.ui.ParSnr.value(),
                'sbr': self.ui.ParSbr.value(),
                'intensity_thresh': self.ui.ParIntThresh.value(),
                'roughness': self.ui.ParRoughness.value(),
                'frequency_cutoff': self.ui.ParFreqCutoff.value(),
                'background_prominence': self.ui.ParBackProm.value(),
                'cutoff_prominence': self.ui.ParFitProminence.value(),
                'cutoff_height': self.ui.ParFitHeight.value(),
                })
        return parameters

def parameters_to_text(parameters):
    outp = [
            'Used parameters',
            '',
            '*Parameters that reject individual spectra',
            ]
    for key, disp in [
            ('snr', 'Max. noise-to-signal ratio'),
            ('sbr', 'Max. background-to-signal ratio'),
            ('intensity_thresh', 'Min. intensity rel. to highest intensity spectrum in dataset'),
            ]:
        outp.append(' *%s: %.5f' % (disp, parameters[key]))
    outp.append('*Parameters for the lowpass frequency filter')
    for key, disp in [
            ('roughness', 'Roughness'),
            ('frequency_cutoff', 'Frequency cutoff'),
            ]:
        outp.append(' *%s: %.5f' % (disp, parameters[key]))
    outp.append('*Parameters for the peak finding')
    for key, disp in [
            ('height', 'Peak min. relative height'),
            ('prominence', 'Peak min. relative prominence'),
            ]:
        outp.append(' *%s: %.5f' % (disp, parameters[key]))
    outp.append('*Parameters for the Gaussian function fitting')
    for key, disp in [
            ('cutoff_height', 'Fit peak min. relative height'),
            ('cutoff_prominence', 'Fit peak min. relative prominence'),
            ]:
        outp.append(' *%s: %.5f' % (disp, parameters[key]))
    return '\n'.join(outp)

if __name__ == "__main__":
    # for pdb to work
    PyQt5.QtCore.pyqtRemoveInputHook()

    #make pyqt threadsafe
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

    app = QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('./spikyboy.jpeg'))
    MainWindow = Main()
    MainWindow.show()
    sys.exit(app.exec_())

