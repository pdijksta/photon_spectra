import os
import sys

import PyQt5
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QMainWindow, QApplication

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

    data_dir = './data/'

    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.SelectFile.clicked.connect(self.select_file(self.ui.Filename))
        self.ui.DoAnalysis.clicked.connect(self.do_analysis)

    def do_analysis(self):
        filename = self.ui.Filename.text().strip()
        if not os.path.isfile(filename):
            raise ValueError('%s does not exist' % filename)
        parameters = self.get_parameters()
        spectrum.analyze_spectrum(filename, parameters)

    def select_file(self, widget):
        def f():
            QFileDialog = PyQt5.QtWidgets.QFileDialog
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            filename, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", './', "Npz files (*.npz);;All Files (*)", options=options)
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

