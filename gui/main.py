import argparse
import os
import sys
import socket
import time
from datetime import datetime
import base64

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import PyQt5
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QMainWindow, QApplication

from h5_storage import saveH5Recursive

if '../../' not in sys.path:
    sys.path.append('../../')

import spectrum
import plot_results


parser = argparse.ArgumentParser()
parser.add_argument('--facility', default='XFEL', choices=('XFEL', 'SwissFEL'))
args = parser.parse_args()

elog = None
if args.facility == 'XFEL':
    import logbook
    default_input_parameters = spectrum.default_input_parameters_xfel
    default_filename = './test_data/20230413-19_10_59_waterflow.npz'
elif args.facility == 'SwissFEL':
    try:
        import elog
    except ImportError:
        print('elog python module unavailable')
    default_input_parameters = spectrum.default_input_parameters_swissfel
    default_filename = '/sf/data/measurements/2023/10/22/fel_spectra_20231022_193353.h5'

if __name__ == '__main__' and (not os.path.isfile('./gui.py') or os.path.getmtime('./gui.ui') > os.path.getmtime('./gui.py')):
    cmd = 'bash ./ui2py.sh'
    print(cmd)
    os.system(cmd)

from gui import Ui_MainWindow


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

        self.key_widget_dict = dict([
            ('height', self.ui.ParHeight),
            ('prominence', self.ui.ParProminence),
            ('snr', self.ui.ParSnr),
            ('sbr', self.ui.ParSbr),
            ('intensity_thresh', self.ui.ParIntThresh),
            ('frequency_cutoff', self.ui.ParFreqCutoff),
            ('background_prominence', self.ui.ParBackProm),
            ('cutoff_prominence', self.ui.ParFitProminence),
            ('cutoff_height', self.ui.ParFitHeight),
            ])

        for key, widget in self.key_widget_dict.items():
            widget.setValue(default_input_parameters[key])

        self.result_dict = None
        if elog is not None:
            self.logbook = elog.open('https://elog-gfa.psi.ch/SwissFEL+commissioning+data/', user='robot', password='robot')
        self.ui.Filename.setText(default_filename)

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

        if args.facility == 'XFEL':
            save_filename = os.path.abspath('./analyzed_data/'+os.path.basename(filename).replace('.npz', '_analyzed.h5'))
        elif args.facility == 'SwissFEL':
            save_filename = os.path.abspath('./analyzed_data/')+datetime.now().strftime('%Y_%m_%d-%H_%M_%S_spike_counting.h5')
        saveH5Recursive(save_filename, result_dict)
        print('Saved %s' % save_filename)

        self.fig = self.do_plot(result_dict, filename)
        self.save_fig(save_filename)
        self.result_dict = result_dict
        self.save_filename = save_filename

    def save_fig(self, save_filename):

        self.fig_savename = save_filename.replace('.h5', '.png')
        self.fig.savefig(self.fig_savename)
        print('Saved %s' % self.fig_savename)

    def do_plot(self, result, filename):
        self.fig, _ = plot_results.standard_plot(filename, result)

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
        if args.facility == 'XFEL':
            succeeded = logbook.send_to_desy_elog('SASE spectrum fit', os.path.basename(self.filename), 'INFO', comment, 'xfellog', image)
            if succeeded:
                print('Logbook save successful.')
            else:
                print('Logbook save unsuccessful.')
        elif args.facility == 'SwissFEL':
            if elog is None:
                print('Cannot save to elog. elog module unavailable.')
                return

            dict_att = {'Author': 'Application: OfflineSpikeCounting', 'Application': 'OfflineSpikeCounting', 'Category': 'Measurement', 'Title': 'Spike counting analysis'}
            self.logbook.post(comment, attributes=dict_att, attachments=[self.fig_savename])

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
        for key, widget in self.key_widget_dict.items():
            parameters[key] = widget.value()
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

