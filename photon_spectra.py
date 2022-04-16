import copy
import inspect
import os
import platform
import sys
import time
from timeit import default_timer as dtimer
import h5py
from numpy import array
import numpy as np
import re
from collections import OrderedDict
import datetime
#IF USING MATPLOTLIB
import matplotlib
# matplotlib is not thread-safe. Must serialize access to matplotlib figs etc.
matplotlib.use('Agg')  # Using a non-interactice backend for threads
# as other GUI Backends must be run from the main thread.

from qtpy import QtCore, QtGui
from qtpy.QtGui import QColor, QFont, QIcon
from qtpy.QtCore import __version__ as QT_VERSION_STR
from qtpy.QtCore import (PYQT_VERSION_STR, QDate, QFile, QIODevice, QMutex, 
                         QObject, Qt, QSize, QThread, QTime, QTimer, Signal, 
                         Slot)
from qtpy.QtWidgets import (QApplication, QCheckBox, QDateTimeEdit, QFileDialog, 
                            QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                            QLineEdit, QMessageBox, QPlainTextEdit, 
                            QProgressBar, QPushButton, QRadioButton, 
                            QSizePolicy, QSpacerItem, QSplitter, QTabWidget, 
                            QToolButton, QVBoxLayout, QWidget)

from bdbase.base import SFMainWindow
from bdbase.enumkind import DAQState, ElogSwissFEL, MsgSeverity, PanelTitle
from bdbase.savehdf import QSaveHDF
from bdbase.sendelog import QSendToELOG
from bdbase.utilities import check_version, check_pyqt_version

from customelements.additional_widgets import (QResultsWidget, QScreenshotPlus,
                                               QHelpLabel)

from customelements.interactive_plot import InteractiveSpectraPlot

from analysis_procedure import AnalysisProcedure

#If using matplotlib
if check_version(PYQT_VERSION_STR, '5.0'):
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
elif check_version(PYQT_VERSION_STR, '4.5'):
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
else:
    raise ImportError("Python requires PyQt4 >= 4.5, found {0}"
                      .format(PYQT_VERSION_STR))

#import pyqtgraph
#print("\npyqtgraph version", pyqtgraph.__version__)

_pymodule = os.path.basename(__file__)
_appname, _appext = _pymodule.split(".")
_appversion = "1.0.0"

_elog_title = "Photon Spectra"

_PROGRESS_BAR_THREAD_START = 0
_PROGRESS_BAR_THREAD_ABORTING = 1
_PROGRESS_BAR_THREAD_ABORTED = 2
_PROGRESS_BAR_THREAD_ERROR = 3
_PROGRESS_BAR_THREAD_END = 100

_RADIO_ALL = "ALL"
_RADIO_INTERACTIVE = "INTERACTIVE"
_RADIO_THRESHOLD = "THRESHOLD"


_INCLUDE_HEADER = True
_INCLUDE_PLOTS = True

def __LINE__():
    """Macro to return the current line number.

    The current line number within the file is used when
    reporting messages to the message logging window.

    Returns: 
    int: Current line number.
    """
    return inspect.currentframe().f_back.f_lineno



class QSaveHDF5Analysed(QSaveHDF):
   
    def save(self):
        #message = self.comment.document().toPlainText()
        #if not message:
        #    self.messagelbl.setText('Please enter a Comment')
        #    return
        user_dict=self.get_data()
        #Do the saving here
        self.h5_filename = user_dict['Destination']
        
        if self.parent.results_dict is None:
            QMessageBox.information(self, 'No analysed data', '''There is no analysed data yet. Please make a measurement first.''')
            self.close() 
            return

        if not self.h5_filename.endswith('.h5'):
            self.h5_filename += '.h5'
            
        newfile = h5py.File(self.h5_filename, 'w')
        print(newfile)
        for key in user_dict:
            if key not in [ 'Destination', 'Author', 'Application']:
                lowkey = key.lower()
                newfile.create_dataset('metadata/'+lowkey.replace(' ', '_'), data=user_dict[key])

        for key in self.parent.results_dict:
            if key not in ['Figure1', 'Figure2', 'Figure3', 'Figure4',  'spectras', 'examples', 'number_of_peaks']:
                newfile.create_dataset('data/' + key, data=array(self.parent.results_dict[key]))

        newfile.create_dataset('examples', data=self.parent.results_dict['examples'])
        newfile.close()

        self.close()

class StartMain(SFMainWindow):

    trigger_plot_widget = Signal(object)
    trigger_progressbar = Signal(int)
    trigger_log_message = Signal(str, str, int, str, dict)
    
    trigger_analysis_message = Signal(str, str)
   
    

    class DAQThread(QObject):
        trigger_daq_thread_event = Signal(dict)
        
        def __init__(self, parent):
            #QThread.__init__(self)
            super().__init__(parent)
            self.parent = parent
            self.cafe = self.parent.cafe
            self.cyca = self.parent.cyca
           
            self.bsd = None
            self.daq_state = None
                     
            self.mutex = QMutex()
            self.bsreader = self.parent.bsreader
            self.bsreader.trigger_daq_timer.connect(self.parent.receive_daq_timer)  
            self.bsreader.trigger_daq.connect(self.receive_daq)
           
            self.icount = 0

            self.max_count_limit = 1000
            self.max_count = 10 #int(self.parent.no_daq_events_le.text())
            #print("no DAQ Events", int(self.parent.no_daq_events_le.text()))
            self.max_count = min(self.max_count, self.max_count_limit)
            self.x_data_list =[None] * self.max_count 
            self.y_data_list =[None] * self.max_count
            self.pid_data_list = [None] * self.max_count
            
        def __del__(self):
            self.wait()

      
        @Slot(object, int)    
        def receive_daq(self, bsd, daq_state):
            self.mutex.lock()
            self.bsd = bsd
            self.has_new_event = True
            self.daq_state = daq_state
            self.mutex.unlock()
            #_pvd = self.bsd.getPV()
            self.max_count = min(int(self.parent.no_daq_events_le.text()), self.max_count_limit)

            if daq_state == self.cyca.ICAFE_DAQ_STOPPED:
                #self.x_data_list =[None] * self.max_count 
                #self.y_data_list =[None] * self.max_count 
                self.icount=0                
                return

            #Reset to new max_count
            if self.icount == 0:
                self.x_data_list = [] # [None] * self.max_count 
                self.y_data_list = [] #[None] * self.max_count 
                self.pid_data_list = [] # [None] * self.max_count
            print ("daq_state=", daq_state)    
            print("bsd_state=", bsd.status)
            
            if daq_state == self.cyca.ICAFE_DAQ_RUN:
                if bsd.status == self.cyca.ICAFE_NORMAL:
                    print("icount", self.icount)
                    self.bsd.show(wf=4)
                    #print('type', type(bsd.getPV(0).value))
                    #numbers = bsd.getPV(1).value
                    #neg_count = sum(1 for number in numbers if number < 0)
                    #print ("negative Int", neg_count)
                    #if neg_count > 0:
                    #    self.bsd.show()
                    #x_spectrum = _pvd[0].value
                    #y_spectrum = _pvd[1].value
                    #print(x_spectrum)
                    #print(y_spectrum)
                    if self.icount < self.max_count: 
                        #print(self.icount)                        
                        energytmp = np.array(bsd.getPV(0).value)                                           
                        A = np.where(energytmp != 0.0)[0]
                        #print (list(energytmp[A]))
                        ##self.x_data_list[self.icount] = list(energytmp[A]) #_x_spectrum
                        ##self.y_data_list[self.icount] = list(np.array(bsd.getPV(1).value)[A]) #_y_spectrum
                        ##self.pid_data_list[self.icount] = bsd.getPV(0).pulseID
                        self.x_data_list.append( list(energytmp[A]) ) #_x_spectrum
                        self.y_data_list.append( list(np.array(bsd.getPV(1).value)[A]) ) #_y_spectrum
                        self.pid_data_list.append( bsd.getPV(0).pulseID )

                        self.icount += 1
                    else:
                        self.parent.daqStop()
                        print("stopped after n events ", self.icount)
                    
     
            if self.parent.daqState in (DAQState.BS_STOP, DAQState.CA_STOP, 
                                        DAQState.BS_PAUSE, DAQState.CA_PAUSE):             
                spectrum_dict = {}
                spectrum_dict['X'] =  self.x_data_list
                spectrum_dict['Y'] =  self.y_data_list
                spectrum_dict['pulse ids'] = self.pid_data_list
                spectrum_dict['icount'] = self.icount
                        
                print("Collected n events ", self.icount)
                self.trigger_daq_thread_event.emit(spectrum_dict)
                        #self.parent.post_daq_action(spectrum_dict)
        
    class AnalysisThread(QThread):
        trigger_results_event = Signal(dict)
        trigger_figures_event = Signal(list)
        trigger_error_event = Signal(str, str)
        def __init__(self, parent):
            QThread.__init__(self)
            self.parent = parent
            self.analysis_procedure = self.parent.analysis_procedure

            self.input_parameters = {} 
            
        def __del__(self):
            self.wait()
            
        def run(self):
            
            #From GUI Input
            for key in self.parent.input.keys():
                if key in ['parameter optimization', 'pulse duration correction', 'spike-width', 'auto-correlation']:
                    self.input_parameters[key] = float(self.parent.input[key].isChecked())
                else:
                    try:
                        inp = float(self.parent.input[key].text())
                        if inp < 0:
                            self.show_error_message()
                            return 
                        elif key in ['snr', 'sbr', 'intensity_thresh', 'frequency_cutoff'] and inp > 1:
                            self.show_error_message()
                            return

                        elif key in ['height', 'prominence', 'cutoff_prominence'] and inp > 100:
                            self.show_error_message()
                            return
                        self.input_parameters[key] = inp
                        
                    except ValueError:
                        self.show_error_message()
                        return 
            print(self.parent.h5_filename)
            self.input_parameters['h5_filename'] = self.parent.h5_filename

            results_dict = self.analysis_procedure.PerformMeas(self.input_parameters)

            #Emit results
            if results_dict is not None:
                self.trigger_results_event.emit(results_dict)
                _mess = "Analysis completed"
                self.parent.trigger_log_message.emit(
                    MsgSeverity.INFO.name, _pymodule, __LINE__(), _mess, {})
                
            figures_dict = self.analysis_procedure.GetFigures()
            self.trigger_figures_event.emit(figures_dict)

        def show_error_message(self):
            self.trigger_error_event.emit('Wrong input format', '''At least one of the input parameters has the wrong format. 

All the inputs need to be a non-negative floating point numbers. There are further restrictions for some of the parameters. For more information on this have a look at the help pages.''')
                
                
    def __init__(self, parent=None):
        super(StartMain, self).__init__(parent, pymodule=_pymodule, 
                                        appversion=_appversion)
        #super().__init__(parent, pymodule=_pymodule, 
        #                 appversion=_appversion)
        
        self.top_mutex = QMutex()
        self.setObjectName("MainWindow")
        self.appname = _appname         
        self.showMessage(MsgSeverity.INFO, _pymodule, __LINE__(), 
                         "Application started")

        self.analysis_thread = None

        self.daq_thread = None
        self.daq_completed = False # set to true once measurement is done
        self.h5_filename = None

        #Use is_dirty nomenclature for tracking other processes
        #that may prevent clean shutdown
        self.is_dirty = False  
        
        self.font_gui = QFont("sans serif")  
        self.font_gui.setPixelSize(16) 
        
        if _INCLUDE_HEADER:
        #Prepare SwissFEL header
            self.pv_pulse_id = "SIN-TIMAST-EVG0:TX-PULSEID"
            self.panel_title = PanelTitle.OPERATION
            self.sf_op_header, self.sf_op_mode = \
                self.setOperationMode(panel=self.panel_title, 
                                      userPV=self.pv_pulse_id)
        
                
        #Add Progress bar to widget
        self.progressbar_simulation = "simulation"
        self.progressbar_standard = "blue"  
        self.progressbar_color = self.progressbar_standard

        self.progressbar = QProgressBar(self)
        self.progressbar.setObjectName(self.progressbar_color)
        self.progressbar.setRange(_PROGRESS_BAR_THREAD_START, 
                                  _PROGRESS_BAR_THREAD_END)
        self.progressbar.setTextVisible(True)     
        self.progressbar.setAlignment(Qt.AlignCenter)
        self.progressbar.setVisible(False)
        self.trigger_progressbar.connect(self.progress_update)
        self.statusbar.addPermanentWidget(self.progressbar, 0)
        
        self.screenshotDest = '/afs/psi.ch/intranet/SF/Beamdynamics/Adrian/photon_spectra/screenshots/'
        self.rawdataDest = '/sf/data/measurements/'
        self.results_dict = None

        #Prepare central widget
        self.central_widget = QTabWidget()
        self.central_widget.setFont(self.font_gui)
        
        #Retrieve log_tab widget
        self.removeDockWidget(self.logDockWidget)
        log_tab_widget = self.logDockWidget.widget()

        #Create plot_tab_widget
        plot_tab_widget = self.get_plot_tab_widget()
        #
        self.infowidth = 500
        self.plotheight = 600
        self.plotwidth = 800

        #Average spectrum/general information subtab

        average_summary_dict = {}
        average_summary_dict['number of spectra:'] = '0'
        average_summary_dict['analysed spectra:'] = '0'
        average_summary_dict['average photon energy [eV]:'] = '0'
        average_summary_dict['intensity fluctuations [%]:'] = '0'
        average_summary_dict['bandwidth FWHM [eV]:'] = '0'

        self.average_results = QResultsWidget(summary_dict=average_summary_dict)
        
        self.average_results_box = self.average_results.group_box('General information', varvalmin=50)
        self.average_results_box.setFixedWidth(self.infowidth)
        self.average_results_box.setFixedHeight(self.plotheight)

        self.plot1_layout.addWidget(self.average_results_box)
        
        self.canvas1 = FigureCanvasQTAgg(self.plain_figure('photon energy eV', 'intensity [arb. units]'))
        self.canvas1.setFixedWidth(self.plotwidth)
        self.canvas1.setFixedHeight(self.plotheight)

        self.plot1_layout.addWidget(self.canvas1)
        self.plot1_layout.setAlignment(Qt.AlignCenter)
        self.plot1_layout.setSpacing(40)
        self.plot1_widget.setLayout(self.plot1_layout)
        self.plot1_widget.setFixedWidth(self.infowidth + self.plotwidth+240)
        self.plot1_widget.setFixedHeight(self.plotheight + 40)

        #Spike statistics tab
        spike_summary_dict = {}
        spike_summary_dict['average spike width FWHM [eV]:'] = '0 \u00B1 0 '
        spike_summary_dict['average spike number:'] = '0'
        spike_summary_dict['noisy spectra:'] = '0/0'

        spike_width_table_dict = {}
        spike_width_table_dict['spike number'] = 'average spike width [eV]'
        for i in range(10):
            spike_width_table_dict[str(i + 1)] = '0 \u00B1 0'

        self.detailed_spike_results = QResultsWidget(summary_dict=spike_summary_dict,
                                                     table_dict=spike_width_table_dict)
        
        self.detailed_spike_results_box = self.detailed_spike_results.group_box('Spike statistics', varvalmin=115)
        self.detailed_spike_results_box.setFixedWidth(self.infowidth)
        self.detailed_spike_results_box.setFixedHeight(self.plotheight)

        self.plot2_layout.addWidget(self.detailed_spike_results_box)
        self.canvas2 = FigureCanvasQTAgg(self.plain_figure('number of spikes', 'fraction of shots [%]'))
        self.canvas2.setFixedWidth(self.plotwidth)
        self.canvas2.setFixedHeight(self.plotheight)
        
        #self.cvlayout2 =QHBoxLayout()
        #self.cvwidget2 = QWidget()
        
        self.plot2_layout.addWidget(self.canvas2)
        self.plot2_layout.setAlignment(Qt.AlignCenter)
        self.plot2_layout.setSpacing(40)
        self.plot2_widget.setLayout(self.plot2_layout)
        
        self.plot2_widget.setFixedWidth(self.infowidth + self.plotwidth+240)
        self.plot2_widget.setFixedHeight(self.plotheight + 40)

        #Examples with spectrum finder
        self.finder_input = {}
        self.current_canvas_idx = 0
        self.canvas_nplts = 0
        self.spectrum_range = np.array([0])
        self.figures_dict = None
        
        self.plot3_layout.addWidget(self.spectrum_finder())

        self.canvas3 = FigureCanvasQTAgg(self.plain_figure('photon energy eV', 'intensity (arb. units)'))
        self.canvas3.setFixedWidth(self.plotwidth)
        self.canvas3.setFixedHeight(self.plotheight)

        self.plot3_layout.addWidget(self.canvas3)
        self.plot3_layout.setAlignment(Qt.AlignCenter)
        self.plot3_layout.setSpacing(40)
        self.plot3_widget.setFixedWidth(self.infowidth + self.plotwidth+240)
        self.plot3_widget.setFixedHeight(self.plotheight + 40)

        #auto-correlation
        auto_summary_dict = {}
        auto_summary_dict['pulse duration FWHM [as]:'] = '0'
        auto_summary_dict['resolution [eV]:'] = '0'

        self.auto_correlation_results = QResultsWidget(summary_dict=auto_summary_dict)
        
        self.auto_correlation_results_box = self.auto_correlation_results.group_box('Auto-correlation', varvalmin=50)
        self.auto_correlation_results_box.setFixedWidth(self.infowidth)
        self.auto_correlation_results_box.setFixedHeight(self.plotheight)
        
        self.plot4_layout.addWidget(self.auto_correlation_results_box)

        self.canvas4 = FigureCanvasQTAgg(self.plain_figure('relative photon energy width', 'second order correlation $G_2$'))
        self.canvas4.setFixedWidth(self.plotwidth)
        self.canvas4.setFixedHeight(self.plotheight)

        self.plot4_layout.addWidget(self.canvas4)
        self.plot4_layout.setAlignment(Qt.AlignCenter)
        self.plot4_layout.setSpacing(40)
        self.plot4_widget.setFixedWidth(self.infowidth + self.plotwidth+240)
        self.plot4_widget.setFixedHeight(self.plotheight + 40)

        #Create measurement_tab_widget
        measurement_tab_widget = self.get_measurement_tab_widget()
        comparison_tab_widget = InteractiveSpectraPlot(self.font_gui)
        self.addGraphToImage(gtitle="Comparison",
                                 graph=comparison_tab_widget.glayout)

        self.central_widget.addTab(measurement_tab_widget, "Measurement")
       
        self.central_widget.addTab(plot_tab_widget, "Plot")
        self.central_widget.addTab(comparison_tab_widget, 'Compare')
        self.central_widget.addTab(log_tab_widget, "Log")
      
        #Analysis tab
        self.analysis_procedure = AnalysisProcedure(self)
        self.trigger_analysis_message.connect(self.show_analysis_error_message)
        
        self.input = {}
        self.expert_setting_height = 210
        
        _start_analysis_button = self.get_start_analysis_button(
            title="Start Analysis",
            tooltip="Start analysis measurement procedure",
            action=self.start_analysis_thread)

        self.analysis_inner_outer_hbox = QHBoxLayout()
        _analysis_inner_layout = QVBoxLayout()


        #_analysis_inner_hbox = QHBoxLayout()
        #_analysis_inner_hbox_widget = QWidget()  
          
        _analysis_inner_layout.addWidget(self.user_input_widgets()) 
     
        #_analysis_inner_hbox_widget.setLayout(_analysis_inner_hbox)
        #_analysis_inner_hbox.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        #_analysis_inner_hbox.setSpacing(20)

        _analysis_inner_layout.addWidget(_start_analysis_button)
      
        #daq box
        self.no_daq_events_widget = None
        self.no_daq_events_layout = None
        self.no_daq_events_label = None
        self.no_daq_events_le = None
        
        daqwidget = self.input_daq()
        _analysis_inner_layout.addWidget(daqwidget)
        
        if self.settings.data["bsread"].get("PV"):
            self.pvbs=self.settings.data["bsread"]["PV"]
            #self.bsModulo=100 
            self.setBSChannels(self.pvbs, daq_start=False)
            
            self.create_daq_thread() 

        _analysis_inner_layout.setSpacing(20)
        _analysis_inner_layout.setAlignment(Qt.AlignCenter | Qt.AlignTop)
        _analysis_inner_layout.setContentsMargins(9, 9, 9, 9)
        self.analysis_inner_outer_hbox.addLayout(_analysis_inner_layout)
        
        #Results group box
        description = f'Analysed dataset: - \n \nMemo:'
        summary_dict = {}
        summary_dict['number of spectra:'] = '0'
        summary_dict['pulse duration FWHM [as]:'] = '0 \u00B1 0'
        #summary_dict['average spike width FWHM [eV]:'] = '0 \u00B1 0'
        summary_dict['average photon energy [eV]:'] = '0'
    
        summary_dict['analysed spectra:'] = '0'
        #summary_dict['maximal pulse duration [as]:'] = '0 \u00B1 0'
        summary_dict['average spike number:'] = '0'
        summary_dict['intensity fluctuations [%]:'] = '0'

        figures_dict = {}
        figures_dict['avg spectra'] = self.plain_figure('photon energy [eV]', 'intensity [arb.units]')
        figures_dict['histogram'] = self.plain_figure('number of spikes', 'fraction of shots [%]')

        self.important_results = QResultsWidget(description=description, summary_dict=summary_dict, figures_dict=figures_dict)
        self.important_results.setBoldFace('minimal pulse duration [as]:')
        self.important_results_box = self.important_results.group_box('Results', maxlabel=3, varvalmin=125, labelwidth=300)

        self.analysis_inner_outer_hbox.addWidget(self.important_results_box)
        
        self.operator_widget.setLayout(self.analysis_inner_outer_hbox)
        #self.operator_widget.setContentsMargins(9, 9, 9, 9)
       
        self.expeditwidth = 135
        self.explabwidth = 350
        _expert_inner_layout = QHBoxLayout()
        #_expert_inner_layout.addWidget(self.no_daq_events_widget)
        _expert_inner_layout.addWidget(self.expert_settings())
        _expert_inner_layout.setSpacing(20)
        _expert_inner_layout.setAlignment(Qt.AlignCenter) 
        self.expert_widget.setLayout(_expert_inner_layout)
        self.expert_widget.setContentsMargins(9, 9, 9, 9)

        # Mainwidow layout
        _layout = QVBoxLayout()
        if _INCLUDE_HEADER:
            _layout.addWidget(self.sf_op_header)
        _layout.addWidget(self.central_widget)

        _window = QWidget()
        _window.setLayout(_layout)
        self.setCentralWidget(_window)
          
        #self.trigger_plot_widget.connect(self.canvas_update)
        self.trigger_log_message.connect(self.log_message_update)
      
        self.statusbar.showMessage("Application ready for user operation")
       
     
    #def update_no_daq_events(self):
    #    self.no_daq_events =  self.no_daq_events_le.text()

################################################################################
#   Add Widgets to Analysis Tab
################################################################################

    def user_input_widgets(self):    
        group_box = QGroupBox("Input Parameters")
        group_box.setObjectName("OUTERLEFT")
        _widget_height = 33
        _parameternames = ["height threshold", 'prominence threshold']
        _initialvalues = ['20', '10']
        _dictnames = ['height', 'prominence']
        _commands = [self.height_help, self.prom_help]
        _vbox = QVBoxLayout()
        _vbox_widget=QWidget()
        
        for name, init, dictname, com in zip(_parameternames, _initialvalues, _dictnames, _commands):
            _widget = QWidget()
            
            _label = QHelpLabel(name)
            _label.setFixedWidth(250)
            _label.setFixedHeight(_widget_height)
            _label.connect(com)
            _label.setFont(self.font_gui)
            self.input[dictname] = QLineEdit("")
            self.input[dictname].setFixedHeight(_widget_height)
            self.input[dictname].setFixedWidth(135)
            self.input[dictname].setProperty("Write", True)
            self.input[dictname].setText(init)
            #_icon = QToolButton()
            #_icon.setIcon(QIcon(':/helpabout.png'))
            #_icon.setAutoRaise(True)
            #_icon.setStyleSheet('background-color: red')
            #_icon.setFixedHeight(_help_height)
            #_icon.setFixedWidth(_help_height)
            _hbox = QHBoxLayout()
            _hbox.addWidget(_label)
            #_hbox.addWidget(_icon)
            _hbox.addWidget(self.input[dictname])
            _hbox.setAlignment(Qt.AlignCenter)
            _widget.setLayout(_hbox)

            _vbox.addWidget(_widget)

        group_box.setContentsMargins(0, 0, 0, 0)
        group_box.setFixedWidth(450)
        group_box.setFixedHeight(200)
        group_box.setFont(self.font_gui)
        group_box.setLayout(_vbox)
        return group_box 
    
    def height_help(self):
        QMessageBox.about(self, 'Height treshold',
"""<html>
<b>Description</b>
<hr>
The height of a peak is defined by its maximum intensity. The height threshold gives the minimal peak height in percentages of the maximum intensity of the whole spectra.
<br><br>
<b>Hint</b>
<hr>
For spectra with low noise to signal ratio the height threshold might be set lower than the default value of 20 percent. But keep in mind that the spike widths are weighted by the average peak height, so these value won't change that much.
</html>""")
    
    def prom_help(self):
        QMessageBox.about(self, 'Prominence treshold',
"""<html>
<b>Decription</b>
<hr>
The prominence of a peak is the minimal distance one has to go down before reaching a higher peak. Intuitively this describes how much a peak stands out. The prominence threshold gives the minimal prominence in percentages of the maximum intensity of the whole spectra.<br><br>
<br><br>
<b>Hint</b>
<hr>
For short, sub-femtosecond pulses we recommend the default value of 10 percent. For longer pulses with many narrow spikes and overlapping spikes a lower value of the prominence gives a more accurate fit. 
</html>""")

    def get_abort_button(self, action=None):
        
        _abort_button = QPushButton("Abort Measurement")
        _abort_button.setObjectName("Abort")
        _abort_button.setFixedWidth(230)
        _abort_button.setFixedHeight(50)
        _abort_button.setToolTip("Aborts measurement procedure")
        _abort_button.clicked.connect(action)
        _abort_button.setVisible(False)    
        return _abort_button

   
    def get_start_button(self, title=None, tooltip=None, 
                                     action=None):
        _start_button = QPushButton(title)
        _start_button.setObjectName("Action")
        _start_button.setFixedWidth(210)
        _start_button.setFixedHeight(50)
        _start_button.setToolTip(tooltip)
        _start_button.clicked.connect(action)
        return _start_button

    def get_start_analysis_button(self, title=None, tooltip=None, 
                                     action=None):
        _hbox = QHBoxLayout()
        _hbox_widget = QWidget()
        self.start_analysis_button = self.get_start_button(
            title=title, tooltip=tooltip, action=action)
        self.abort_analysis_button = self.get_abort_button(
            action=self.receive_abort_analysis) 

     
        _hbox.addWidget(self.start_analysis_button)
        _hbox.addWidget(self.abort_analysis_button)
           
        _hbox_widget.setLayout(_hbox)
        _hbox.setAlignment(Qt.AlignLeft)
      
       
        return _hbox_widget
    
    def input_daq(self):
        _groupbox = QGroupBox('Data acquisition')
        _groupbox.setObjectName('OUTERLEFT')
        _gblayout = QVBoxLayout()
        
        self.no_daq_events_layout = QHBoxLayout()
        self.no_daq_events_label = QLabel("No. DAQ Events:")
        self.no_daq_events_label.setFont(self.font_gui)
        self.no_daq_events_le = QLineEdit()
        self.no_daq_events_le.setObjectName("Write")
        self.no_daq_events_le.setText('200')
        self.no_daq_events_le.setAlignment(Qt.AlignTop | Qt.AlignRight)
       
        self.no_daq_events_layout.addWidget(self.no_daq_events_label)
        self.no_daq_events_layout.addStretch()
        self.no_daq_events_layout.addWidget(self.no_daq_events_le)
        self.no_daq_events_layout.setContentsMargins(0,0,0,0)

        self.no_daq_events_le.setFixedWidth(220)
        self.no_daq_events_le.setFixedHeight(35)
        self.no_daq_events_label.setFixedWidth(170)
        self.no_daq_events_label.setFixedHeight(35)
        self.no_daq_events_layout.setAlignment(Qt.AlignCenter | Qt.AlignRight)
       
        _gblayout.addLayout(self.no_daq_events_layout)
        
        _tmplayout = QHBoxLayout()
        _tmplabel = QLabel("Filename:")
        _tmplabel.setFont(self.font_gui)
        self.daq_filename_le = QPlainTextEdit()
        self.daq_filename_le.setObjectName("Write")
        self.daq_filename_le.insertPlainText('photon_spectra_measurement')
        self.daq_filename_le.setFont(self.font_gui)
       
        _tmplayout.addWidget(_tmplabel)
        _tmplayout.addStretch()
        _tmplayout.addWidget(self.daq_filename_le)
        _tmplayout.setContentsMargins(0,0,0,0)
        
        self.daq_filename_le.setFixedWidth(220)
        self.daq_filename_le.setFixedHeight(70)
        _tmplabel.setFixedWidth(170)
        _tmplabel.setFixedHeight(35)
        
        _gblayout.addLayout(_tmplayout)
        
        _tmplayout = QHBoxLayout()
        _tmplabel = QLabel("Memo:")
        _tmplabel.setFont(self.font_gui)
        self.daq_memo_le = QPlainTextEdit()
        self.daq_memo_le.setObjectName("Write")
        self.daq_memo_le.insertPlainText('')
        self.daq_memo_le.setFont(self.font_gui)
       
        _tmplayout.addWidget(_tmplabel)
        _tmplayout.addStretch()
        _tmplayout.addWidget(self.daq_memo_le)
        _tmplayout.setContentsMargins(0,0,0,0)
        
        self.daq_memo_le.setFixedWidth(220)
        self.daq_memo_le.setFixedHeight(105)
        _tmplabel.setFixedWidth(170)
        _tmplabel.setFixedHeight(35)
        
        _gblayout.addLayout(_tmplayout)
        _gblayout.addStretch()
        _gblayout.setSpacing(15)
        _gblayout.setAlignment(Qt.AlignHCenter)
        _gblayout.setContentsMargins(20, 20, 20, 20)
        _groupbox.setLayout(_gblayout)
        _groupbox.setFont(self.font_gui)
        _groupbox.setFixedHeight(290)
        _groupbox.setFixedWidth(450)
        return _groupbox

    def expert_settings(self):
        group_box = QGroupBox("Advanced settings")
        group_box.setObjectName("OUTERLEFT")
        _expert_group_box_layout = QHBoxLayout()
        
        _expert_inner_layout1 = QVBoxLayout()
        _expert_inner_layout1.setSpacing(20)
        _ticks_layout = QHBoxLayout()
        
        _ticks_layout.addWidget(self.method_settings())
        _ticks_layout.addWidget(self.analysis_settings())
        _expert_inner_layout1.addLayout(_ticks_layout)
        _expert_inner_layout1.addWidget(self.noise_settings())
        
        _expert_group_box_layout.addLayout(_expert_inner_layout1)
        _expert_inner_layout2 = QVBoxLayout()
        _expert_inner_layout2.setSpacing(20)
        _expert_inner_layout2.addWidget(self.lowpass_filter_settings())
        _expert_inner_layout2.addWidget(self.peakfinder_settings())

        _expert_group_box_layout.addLayout(_expert_inner_layout2)
        _expert_group_box_layout.setAlignment(Qt.AlignCenter)

        group_box.setFixedWidth(1000)
        group_box.setFixedHeight(500)
        group_box.setContentsMargins(0, 0, 0, 0)
        group_box.setFont(self.font_gui)
        group_box.setLayout(_expert_group_box_layout)
        return group_box

    def lowpass_filter_settings(self):
        #group_box = QGroupBox("Lowpass filter settings")
        
        _widget_height = 30

        _parameternames = ['roughness', 'frequency cutoff']
        _initialvalues = ['200', '0.0001']
        _dictnames =['roughness', 'frequency_cutoff']

        _box = QVBoxLayout()
        _box_widget = QWidget()
        _box.setSpacing(10)
        _title = QHelpLabel('Lowpassfilter settings')
        _title.setStyleSheet("border-bottom-width: 1px; border-bottom-style: solid; border-radius: 0px;")
        _title.connect(self.help_lowpass)
        _title.setFont(self.font_gui)
        _box.addWidget(_title)
        for name, init, dictname in zip(_parameternames, _initialvalues, _dictnames):
            _widget = QWidget()

            _label = QLabel(name + ':')
            _label.setFixedWidth(self.explabwidth)
            _label.setFont(self.font_gui)
            self.input[dictname] = QLineEdit("")
            self.input[dictname].setFixedHeight(_widget_height)
            self.input[dictname].setProperty("Write", True)
            self.input[dictname].setText(init)
            self.input[dictname].setFixedWidth(self.expeditwidth)
            _hbox = QHBoxLayout()
            _hbox.setSpacing(10)
            _hbox.addWidget(_label)
            _hbox.addWidget(self.input[dictname])
            _widget.setLayout(_hbox)

            _box.addWidget(_widget)

        _box.setAlignment(Qt.AlignTop)
        _box_widget.setLayout(_box)
        
        #group_box.setFixedWidth(400)
        #group_box.setFixedHeight(140)
        #group_box.setContentsMargins(0, 0, 0, 0)
        #group_box.setFont(self.font_gui)
        #group_box.setLayout(_box)
        return _box_widget

    def help_lowpass(self):
        QMessageBox.about(self, 'Low-pass filter settings', 
"""<html>
<b>Description</b>
<hr>
The low-pass filter settings determine the parameters of the adaptive butterworth filter. You can control the following parameter.
<br><br> 
<b>roughness:</b>
<br>
The rougness determines how smooth the filtered spectra is. A high value of the roughness means that the filtered spectrum is smoother. This parameter is optimized in the parameter optimization.
<br>
<br><br>
<b>frequency cutoff:</b>
<br>
The frequency cutoff parameter determines the initial value of the frequency cutoff in the adaptive butterworth filter. It is given in terms of the Nyquist-frequency. The value of this parameter has to be between 0 and 1 and it is optimized during the parameter optimization.
<br><br>
<b>Hint</b>
<hr>
If lowering the prominence threshold couldn't help you to find more spikes the roughness parameter might be to high. Don't forget to turn off the parameter optimization before manually trying out lower roughness parameter. The lowest roughness value considered in the  parameter optimization is five.
</html>""")
    def peakfinder_settings(self):
        #group_box = QGroupBox("Peakfinder settings")
        #group_box.setObjectName("OUTERLEFT")
        _widget_height = 30

        _parameternames = ['background factor', 'cutoff prominence', 'cutoff height']
        _initialvalues = ['2.5', '5', '70']
        _dictnames =['background_prominence', 'cutoff_prominence', 'cutoff_height']

        _vbox = QVBoxLayout()
        _vbox_widget = QWidget()
        _vbox.setSpacing(10)
        _title = QHelpLabel('Peakfinder settings')
        _title.setStyleSheet("border-bottom-width: 1px; border-bottom-style: solid; border-radius: 0px;")
        _title.setFont(self.font_gui)
        _title.connect(self.help_peakfinder)
        _vbox.addWidget(_title)
        for name, init, dictname in zip(_parameternames, _initialvalues, _dictnames):
            _widget = QWidget()

            _label = QLabel(name + ':')
            _label.setFixedWidth(self.explabwidth)
            _label.setFont(self.font_gui)
            self.input[dictname] = QLineEdit("")
            self.input[dictname].setFixedHeight(_widget_height)
            self.input[dictname].setProperty("Write", True)
            self.input[dictname].setText(init)
            self.input[dictname].setFixedWidth(self.expeditwidth)
            _hbox = QHBoxLayout()
            _hbox.addWidget(_label)
            _hbox.addWidget(self.input[dictname])
            _widget.setLayout(_hbox)

            _vbox.addWidget(_widget)
        
        _vbox.setAlignment(Qt.AlignTop)
        _vbox_widget.setLayout(_vbox)

        #group_box.setContentsMargins(0, 0, 0, 0)
        #group_box.setFixedWidth(400)
        #group_box.setFixedHeight(self.expert_setting_height)
        #group_box.setFont(self.font_gui)
        #group_box.setLayout(_vbox)
        return _vbox_widget

    def help_peakfinder(self):
        QMessageBox.about(self, 'Peakfinder settings',
"""<html>
<b>Description</b>
<hr>
The peakfinder settings determine which peaks get fitted and at which points the spectra are sliced. Next to height and prominence parameter in the operator panel, you can control the following parameters.
<br><br> 
<b>background factor:</b>
<br>
The minimal peak height needs to be higher than the background factor times the intensity at the border of the spectrum window. 
<br><br>
<b>cutoff prominence:</b>
<br>
For fitting the spectrum is cutoff at a local minima between the edge peak and the border of the spectrum window. The cutoff prominence determines the minimal prominence of this minima in percentages of the maximal intensity.
<br><br>
<b>cutoff height:</b>
<br>
For fitting the spectrum is cutoff at a local minima between the edge peak and the border of the spectrum window. The cutoff height determines the maximal height of this minima in percentages of the edge peak's height.
 </html>""")

    def noise_settings(self):
        #group_box = QGroupBox("Noise detection settings")
        #group_box.setObjectName("OUTERLEFT")
        _widget_height = 30

        _parameternames = ['noise to signal  ratio', 'background to signal ratio', 'relative intensity threshold']
        _initialvalues = ['0.05', '0.15', '0.4']
        _dictnames =['snr', 'sbr', 'intensity_thresh']

        _vbox = QVBoxLayout()
        _vbox_widget = QWidget()
        _vbox.setSpacing(10)
        _title = QHelpLabel('Noise detection settings')
        _title.setStyleSheet("border-bottom-width: 1px; border-bottom-style: solid; border-radius: 0px;")
        _title.setFont(self.font_gui)
        _title.connect(self.help_noise)
        _vbox.addWidget(_title)
        for name, init, dictname in zip(_parameternames, _initialvalues, _dictnames):
            _widget = QWidget()

            _label = QLabel(name + ':')
            _label.setFixedWidth(self.explabwidth)
            _label.setFont(self.font_gui)
            self.input[dictname] = QLineEdit("")
            self.input[dictname].setFixedHeight(_widget_height)
            self.input[dictname].setProperty("Write", True)
            self.input[dictname].setText(init)
            self.input[dictname].setFixedWidth(self.expeditwidth)
            _hbox = QHBoxLayout()
            _hbox.addWidget(_label)
            _hbox.addWidget(self.input[dictname])
            _widget.setLayout(_hbox)

            _vbox.addWidget(_widget)
        
        _vbox.setAlignment(Qt.AlignTop)
        _vbox_widget.setLayout(_vbox)

        #group_box.setContentsMargins(0, 0, 0, 0)
        #group_box.setFixedWidth(400)
        #group_box.setFixedHeight(self.expert_setting_height)
        #group_box.setFont(self.font_gui)
        #group_box.setLayout(_vbox)
        return _vbox_widget
    def help_noise(self):
        QMessageBox.about(self, 'Noise detection settings', 
"""<html>
<b>Description</b>
<hr>
The noise detection settings determine which spectra are sorted out. You can control the following parameters:
<br><br>
<b>noise to signal ratio:</b>
<br>
Maximal noise to signal ratio of a spectrum considered for analysis. The noise is defined as the rms difference between the low-pass filtered spectrum and the orginal spectrum. The signal amplitude is the maximal intensity of the spectrum.
<br><br>
<b>background to signal ratio:</b>
<br>
Maximal background to signal ratio of a spectrum considered for analysis. The background is defined as the intensity at the edges of spectra window. The signal amplitude is the maximal intensity of the spectrum.
<br><br>
<b>relative intensity threshold:</b>
<br>
Minimal integrated intensity in percentages of the average integrated intensity.
<br><br>
<b>Hint</b>
<hr>
If you want to analyse all spectra then set the thresholds to
<ul>
<li>noise to signal ratio: 1 </li>
<li>background to signal ratio: 1 </li>
<li>relative intensity threshold: 0 </li>
</ul>
</html>""")
    def method_settings(self):
        #group_box = QGroupBox("Analysis settings")
        #group_box.setObjectName("OUTERLEFT")
        _widget_height = 30

        _parameternames = ['spike-width', 'auto-correlation']

        _vbox = QVBoxLayout()
        _vbox_widget = QWidget()
        _vbox.setSpacing(10)
        _title = QHelpLabel('Method')
        _title.setStyleSheet("border-bottom-width: 1px; border-bottom-style: solid; border-radius: 0px;")
        _title.connect(self.help_method)
        _title.setFont(self.font_gui)

        _vbox.addWidget(_title)
        for name in _parameternames:
            self.input[name] = QCheckBox()
            self.input[name].setFixedHeight(_widget_height)
            self.input[name].setChecked(True)
            self.input[name].setText(name)
            self.input[name].setFont(self.font_gui)
            _vbox.addWidget(self.input[name])
        
        _vbox.setAlignment(Qt.AlignTop)
        _vbox_widget.setLayout(_vbox)
        
        self.input['auto-correlation'].setChecked(False)
        
        #group_box.setContentsMargins(0, 0, 0, 0)
        #group_box.setFixedWidth(400)
        #group_box.setFixedHeight(140)
        #group_box.setFont(self.font_gui)
        #group_box.setLayout(_vbox)

        return _vbox_widget

    def analysis_settings(self):
        #group_box = QGroupBox("Analysis settings")
        #group_box.setObjectName("OUTERLEFT")
        _widget_height = 30

        _parameternames = ['parameter optimization', 'pulse duration correction']

        _vbox = QVBoxLayout()
        _vbox_widget = QWidget()
        _vbox.setSpacing(10)
        _title = QHelpLabel('Analysis settings')
        _title.setStyleSheet("border-bottom-width: 1px; border-bottom-style: solid; border-radius: 0px;")
        _title.connect(self.help_analysis)
        _title.setFont(self.font_gui)

        _vbox.addWidget(_title)
        for name in _parameternames:
            self.input[name] = QCheckBox()
            self.input[name].setFixedHeight(_widget_height)
            self.input[name].setChecked(True)
            self.input[name].setText(name)
            self.input[name].setFont(self.font_gui)
            _vbox.addWidget(self.input[name])
        
        _vbox.setAlignment(Qt.AlignTop)
        _vbox_widget.setLayout(_vbox)

        #group_box.setContentsMargins(0, 0, 0, 0)
        #group_box.setFixedWidth(400)
        #group_box.setFixedHeight(140)
        #group_box.setFont(self.font_gui)
        #group_box.setLayout(_vbox)

        return _vbox_widget

    def help_method(self):
        QMessageBox.about(self, 'Method', 
"""<html>
<b>Description</b>
<hr>
You can analyse a dataset with the following methods.
<br><br>
<b>spike-width:</b>
<br> 
Determines the pulse duration with the weighted average spectral spike width.
<br><br>
<b>auto-correlation:</b>
<br>
Determines the pulse duration by fitting of the auto-correlation function.
<br><br>
<b>Hint</b>
<hr>
When you disable both methods only the average spectrum gets analysed, from which one can obtain the average photon energy, the bandwidth and the intensity fluctuations. Furthermore, you can have a look at the raw spectra under plot and then examples.
""" )

    def help_analysis(self):
        QMessageBox.about(self, 'Analysis settings', 
"""<html>
<b>Description</b>
<hr>
The analysis setting determines how to spectra are analysed. You can control the following parameters.
<br><br>
<b>parameter optimization:</b>
<br> 
Turns on the feedback loop in the beginning of analysis to determine optimal roughness and frequency cutoff parameter. Note for the proper behaviour of the lowpass filter the frequency cutoff is also optimized for every shot individually, so the parameter optimization only determines the initial value (decreases time consumption).
<br><br>
<b>pulse duration correction:</b>
<br>
Correction of the average pulse duration by multiplying it with a factor obtained from the benchmarking of the algorithm with genesis simulations. The factor depends on the number of spikes.
</html>""" )

    def get_plot_tab_widget(self):
        #self.plot widgets and self.layouts for later use
        plot_widget = QWidget() #self.central_widget)
        self.plot_qtab_widget = QTabWidget()  #plot_widget)
        self.plot1_widget = QWidget()#self.plot_qtab_widget)
        self.plot2_widget = QWidget()#self.plot_qtab_widget)
        self.plot3_widget = QWidget()#self.plot_qtab_widget)
        self.plot4_widget = QWidget()#self.plot_qtab_widget)
        self.plot1_layout = QHBoxLayout() #self.plot1_widget)
        self.plot1_widget.setLayout(self.plot1_layout)
        self.plot2_layout = QHBoxLayout() #self.plot2_widget)
        self.plot2_widget.setLayout(self.plot2_layout)
        self.plot3_layout = QHBoxLayout() #self.plot3_widget)
        self.plot3_widget.setLayout(self.plot3_layout)
        self.plot4_layout = QHBoxLayout() #self.plot3_widget)
        self.plot4_widget.setLayout(self.plot4_layout)

        self.plot_qtab_widget.addTab(self.plot1_widget, "Average spectra")
        self.plot_qtab_widget.addTab(self.plot2_widget, "Histogram")
        self.plot_qtab_widget.addTab(self.plot4_widget, 'Auto-correlation')
        self.plot_qtab_widget.addTab(self.plot3_widget, "Examples")
        self.plot_qtab_widget.setFont(self.font_gui) 
        plot_layout = QHBoxLayout() #plot_widget)
        plot_layout.addWidget(self.plot_qtab_widget)
        plot_widget.setLayout(plot_layout)
        plot_layout.setContentsMargins(9, 9, 9, 9)

        return plot_widget

    def plain_figure(self, xlabel, ylabel):
        fig = matplotlib.figure.Figure()
        subplot = fig.add_subplot(111)
        subplot.set_xlabel(xlabel)
        subplot.set_ylabel(ylabel)
        fig.set_tight_layout(True)
        return fig

    def onToggleLeftAction(self, is_checked):
        print(self.h5_filename)
        self.update_plot(-1)

    def onToggleRightAction(self, is_checked):
        self.update_plot(1)

    def update_plot(self, idx):

        if self.figures_dict is None:
            return

        self.current_canvas_idx += idx
        if self.current_canvas_idx >= self.canvas_nplts:
            self.current_canvas_idx = 0
        elif self.current_canvas_idx < 0:
            self.current_canvas_idx = self.canvas_nplts-1

        if self.canvas3 is not None:
            self.plot3_layout.removeWidget(self.canvas3)
            self.canvas3.deleteLater()

        if self.figures_dict[self.spectrum_range[self.current_canvas_idx]] is not None:            
            self.canvas3 = FigureCanvasQTAgg(self.figures_dict[self.spectrum_range[self.current_canvas_idx]])
            
            self.canvas3.draw()
            self.canvas3.setFixedWidth(self.plotwidth)
            self.canvas3.setFixedHeight(self.plotheight)

            self.plot3_layout.addWidget(self.canvas3)
            self.plot3_layout.setContentsMargins(10, 10, 10, 10)
            self.plot3_widget.setLayout(self.plot3_layout)
            self.plot3_widget.show()

        else:
            self.canvas3 = None
        if self.canvas3 is not None:
            self.addGraphToImage(gtitle="Example",
                                 graph=self.plot3_widget)
           

    def spectrum_finder(self):
        group_box = QGroupBox("Spectrum finder")
        group_box.setObjectName("OUTERLEFT")
        _widget_height = 30
        _widget_width = 300
        _parameternames = ['spectrum number', 'peak number']
        _initialvalues = ['', '']

        _vbox = QVBoxLayout()
        _vbox_widget = QWidget()
        
        for name, init in zip(_parameternames, _initialvalues):
            _widget = QWidget()

            _label = QLabel(name + ':')
            _label.setFixedWidth(150)
            _label.setFont(self.font_gui)
            self.finder_input[name] = QLineEdit("")
            self.finder_input[name].setFixedHeight(30)
            self.finder_input[name].setFixedWidth(100)
            self.finder_input[name].setProperty("Write", True)
            self.finder_input[name].setText(init)
            _hbox = QHBoxLayout()
            _hbox.addWidget(_label)
            _hbox.addWidget(self.finder_input[name])
            _vbox.addLayout(_hbox)
        
        button_widget = QWidget()
        button_layout = QHBoxLayout()
        search_button = QPushButton("Search")
        search_button.setFont(self.font_gui)

        search_button.setObjectName("Search")
        search_button.setFixedWidth(175)
        search_button.setFixedHeight(40)
        search_button.setToolTip("Searches for spectra with given constraints")
        search_button.clicked.connect(self.search_spectra) 
        button_layout.addWidget(search_button)
        button_layout.setAlignment(Qt.AlignCenter)
        
        _vbox.addLayout(button_layout)
        
        arrow_button_widget = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 20, 0, 0)
        layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        button = QToolButton()
        button.setArrowType(QtCore.Qt.LeftArrow)
        button.clicked.connect(self.onToggleLeftAction)
        layout.addWidget(button)

        button = QToolButton()
        button.setArrowType(QtCore.Qt.RightArrow)
        button.clicked.connect(self.onToggleRightAction)
        layout.addWidget(button)
        layout.setAlignment(Qt.AlignCenter)

        _vbox.addLayout(layout)
        _vbox.setAlignment(Qt.AlignHCenter)
        _vbox.addStretch()
        _vbox.setSpacing(20)
        _vbox_widget.setLayout(_vbox)
        
        group_box.setContentsMargins(0, 0, 0, 0)
        group_box.setFixedWidth(500)
        group_box.setFixedHeight(self.plotheight)
        group_box.setFont(self.font_gui)
        group_box.setLayout(_vbox)
        return group_box
    
    def search_spectra(self):
        if self.results_dict is None:
            return

        edit1 = self.finder_input['spectrum number'].text()
        edit2 = self.finder_input['peak number'].text()
        
        srtemp = self.readout(edit1, len(self.figures_dict))
        peak_range = self.readout(edit2, int(np.max(self.results_dict['number_of_peaks'])))
                
        found_indices = []
        for i in srtemp:
            if self.results_dict['number_of_peaks'][i] in peak_range:
                found_indices.append(i)
        if len(found_indices) == 0:
            self.spetrum_range = np.arange(len(self.figures_dict))
        else:
            self.spectrum_range = np.array(found_indices)
                               
        self.current_canvas_idx = 0
        self.canvas_nplts = len(self.spectrum_range)
        self.update_plot(0)

    def readout(self, entrystr, maxnumber):
        
        search1 = re.search('-', entrystr)
        search2 = list(re.finditer(',', entrystr))
        if entrystr == '':
            spectrum_range = range(maxnumber)
        elif search1 is not None:
            spectrum_range = range(int(entrystr[:search1.start()]), min(int(entrystr[search1.end():]) + 1, maxnumber))
        elif len(search2) != 0:
            left = 0
            spectrum_range = []
            for match in search2:
                right = match.start()
                try:
                    if int(entrystr[left:right]) < maxnumber:
                        spectrum_range.append(int(entrystr[left:right]))
                except ValueError:
                    pass
                left = match.end()
            try:
                if int(entrystr[left:]) < maxnumber:
                    spectrum_range.append(int(entrystr[left:]))
            except ValueError:
                pass
            if len(spectrum_range) == 0:
                spectrum_range = range(maxnumber)
        else:
            try:
                spectrum_range = [int(entrystr)]
            except ValueError:
                #add message box
                spectrum_range = range(maxnumber)

        return array(spectrum_range)

    def get_measurement_tab_widget(self):
        measurement_widget = QWidget() #self.central_widget)
        measurement_tab_widget = QTabWidget() #measurement_widget)
        self.operator_widget = QWidget() #measurement_tab_widget)
        self.expert_widget = QWidget()  #measurement_tab_widget)
        #operator_layout = QHBoxLayout(self.operator_widget)
        #self.operator_widget.setLayout(operator_layout)
        #expert_layout = QHBoxLayout(self.expert_widget)
        #self.expert_widget.setLayout(expert_layout)
        measurement_tab_widget.addTab(self.operator_widget, "Operator")
        measurement_tab_widget.addTab(self.expert_widget, "Expert")  
        measurement_tab_widget.setFont(self.font_gui)
        measurement_layout = QHBoxLayout() #measurement_widget)
        measurement_layout.addWidget(measurement_tab_widget)
        measurement_widget.setLayout(measurement_layout)
        measurement_layout.setContentsMargins(9, 9, 9, 9)
        
        return measurement_widget

    def start_analysis_thread(self):
        if self.analysis_thread is not None:
            if self.analysis_thread.isRunning():
                qm = QMessageBox()
                qm.setText("Measurement already in progress...")
                qm.exec()
                return
        if self.daqState in (DAQState.BS, DAQState.CA):
            qm = QMessageBox()
            qm.setText(("Data Acquistion in progress ...\n" +
                        "Please wait until DAQ has completed, else\n" +
                        "Stop or Pause DAQ to enable analysis procedure"))
            qm.exec()
            return 


        self.analysis_thread = self.AnalysisThread(self)
        self.analysis_thread.trigger_results_event.connect(
            self.receive_analysis_results)
        self.analysis_thread.trigger_error_event.connect(self.show_analysis_error_message)
        self.analysis_thread.trigger_figures_event.connect(self.receive_figures)
        self.analysis_thread.started.connect(self.analysis_thread_started)
        self.analysis_thread.finished.connect(self.analysis_thread_finished)
        self.progressbar.setVisible(True)
 
        #Un comment this section for when hd5 file is to be loaded 

        
        _file_read_error = False

        if self.filename is not None:
            self.h5_filename = os.path.expanduser(self.filename) #
            if not os.path.isfile(self.h5_filename):
                _file_read_error = True
        else:
            _file_read_error = True

        if _file_read_error:                
            qm = QMessageBox()
            qm.setText( ("No HDF5 file loaded. \nSelect 'Open File' " + 
                         "from the File Menu to load file for analysis"))
            qm.exec()
            return



   
        self.analysis_thread.start()
        QApplication.processEvents()
        self.abort_analysis_button.setVisible(True)
        self.abort_analysis_button.setEnabled(True)
        self.abort_analysis_button.setToolTip(
            'Aborts measurement procedure. Data are not saved')
        self.start_analysis_button.setEnabled(False)
        self.start_analysis_button.setText('Measuring...')

        self.start_analysis_button.setToolTip(
            'Measurement in progress...')



    def analysis_thread_started(self):              
        if self.progressbar is not None: 
            self.reset_progress_bar()
            self.progressbar.setValue(10)
            QApplication.processEvents()
        _mess = "Analysis measurement started"    
        self.statusbar.showMessage(_mess)
        self.showMessage("INFO", _pymodule, __LINE__(), _mess)
     
    def analysis_thread_finished(self):         
        if self.progressbar is not None:
            self.reset_progress_bar()
            self.progressbar.setVisible(False)
       
        self.abort_analysis_button.setVisible(False)
        self.start_analysis_button.setEnabled(True)
        self.start_analysis_button.setText('Start Analysis')
        self.start_analysis_button.setToolTip('Start analysis measurement')                    
        QApplication.processEvents()


    def post_daq_action(self, spectrum_dict):
        #save to hdf5
       
        #http://cafe.psi.ch/cython.html#scn
        #scroll to Retrieving structured data
        #print first event in list
        print("No of events analysed", spectrum_dict['icount']) # len(spectrum_dict['X']))
        #_pvdata = spectrum_dict['X'][0]
        #_pvdata.show()
        
        now = datetime.datetime.now()
        datestr = now.strftime("%Y-%m-%d")
        timestr = now.strftime("%H:%M:%S")
        
        if self.daq_filename_le.toPlainText()+'.h5' in os.listdir(self.rawdataDest+str(now.year).zfill(2)+'/'+str(now.month).zfill(2)+'/'+str(now.day).zfill(2)+'/'):
            filename = str(now.year).zfill(2)+'/'+str(now.month).zfill(2)+'/'+str(now.day).zfill(2)+'/'+self.daq_filename_le.toPlainText()+'_'+now.strftime("%Y-%m-%d_%H:%M:%S")
        else:
            filename = str(now.year).zfill(2)+'/'+str(now.month).zfill(2)+'/'+str(now.day).zfill(2)+'/'+self.daq_filename_le.toPlainText()
                                                
        with h5py.File(self.rawdataDest+filename+'.h5', 'w') as f:
            _icount = spectrum_dict['icount']
            f.create_dataset('x-axis', data=np.array(spectrum_dict['X'][0]))
            f.create_dataset('y-axis', data=np.array(spectrum_dict['Y']))
            f.create_dataset('date', data=datestr)
            f.create_dataset('time', data=timestr)
            f.create_dataset('memo', data=self.daq_memo_le.toPlainText())
            f.create_dataset('pulse ids', data=np.array(spectrum_dict['pulse ids']))
        print(spectrum_dict['pulse ids'])
        self.filename = self.rawdataDest+filename+'.h5'
        #print(self.h5_filename)
        #(spectrum_dict['Y'][0]).show()

    def create_daq_thread(self):
        if self.daq_thread is not None:
            if self.daq_thread.isRunning():
                qm = QMessageBox()
                qm.setText("Daq already in progress...")
                qm.exec()
                return
        self.daq_thread = self.DAQThread(self)
        self.daq_thread.trigger_daq_thread_event.connect(self.post_daq_action)
        #self.daq_thread.started.connect(self.daq_thread_started)
        #self.daq_thread.finished.connect(self.daq_thread_finished)
        #self.progressbar.setVisible(True)
 
      
        #self.daq_thread.start()
        #QApplication.processEvents()spectrum_dict['X']

    '''    
    def daq_thread_started(self):
        print("DAQ thread started")
        if self.progressbar is not None: 
            self.reset_progress_bar()
            self.progressbar.setValue(10)
            QApplication.processEvents()

    def daq_thread_finished(self):  
        print("DAQ thread ended")
        if self.progressbar is not None:
            self.reset_progress_bar()
            self.progressbar.setVisible(False)
        _mess = "DAQ completed"    
        self.statusbar.showMessage(_mess)
        self.showMessage("INFO", _pymodule, __LINE__(), _mess)
        #Optionally stop DAQ
        self.daqStop()       
    '''    

    def reset_progress_bar(self):
        if self.progressbar is not None:             
            self.progressbar_color =  self.progressbar_standard
            self.progressbar.setObjectName(self.progressbar_color)
            self.progressbar.style().polish(self.progressbar)
            

################################################################################
#   Receive signals
################################################################################


    @Slot()
    def receive_abort_analysis(self):
        
        if self.analysis_procedure.abort:
            return
        self.abort_analysis_button.setEnabled(False)
        self.analysis_procedure.trigger_abort.emit()
        self.trigger_progressbar.emit(_PROGRESS_BAR_THREAD_ABORTING)
        _mess = 'Aborting analysis measurement...'
       
        self.abort_analysis_button.setToolTip(_mess)
        self.start_analysis_button.setToolTip(
            'Start button will be re-enabled soon')

    @Slot(str, str, int, str, dict)
    def log_message_update(self, sev, mod, line, mess, options):
        '''Receive message from thread for routing to log window'''   
        if not options:
            self.showMessage(sev, mod, line, mess)
        else:
            self.showMessage(sev, mod, line, mess, options)
    
        self.statusbar.showMessage(mess)

    @Slot(int)
    def progress_update(self,  int_value):
       
        if int_value == _PROGRESS_BAR_THREAD_START:
            self.progressbar.setFormat("Measurement started")
            self.progressbar.setValue(int_value)
        elif int_value == _PROGRESS_BAR_THREAD_ABORTING:
            self.progressbar.setFormat(
                "Aborting procedure at the next available break point")
            self.progressbar.setObjectName("abort")
            self.progressbar.style().polish(self.progressbar)  
            QApplication.processEvents()
            time.sleep(2)
            _mess =  "Measurement procedure aborted"
            self.showMessage("INFO", _pymodule, __LINE__(), _mess)
            self.statusbar.showMessage(_mess)
        elif int_value == _PROGRESS_BAR_THREAD_ABORTED:
            self.progressbar.setFormat("Procedure aborted")
            self.progressbar.setObjectName("abort")
            self.progressbar.style().polish(self.progressbar)
            _mess =  "Measurement procedure aborted"
            self.showMessage("INFO", _pymodule, __LINE__(), _mess)
            self.statusbar.showMessage(_mess)
        elif int_value == _PROGRESS_BAR_THREAD_ERROR:
            self.progressbar.setFormat("No data returned!")
            self.progressbar.setObjectName("abort")
            self.progressbar.style().polish(self.progressbar)
        elif int_value == _PROGRESS_BAR_THREAD_END:
            self.progressbar.setFormat("Measurement completed")
            self.progressbar.setValue(int_value)
            self.daq_completed = True
            time.sleep(2)
        else:
            self.progressbar.setFormat("Measurement in progress...")
            QApplication.processEvents()
            self.progressbar.setValue(int_value)

        QApplication.processEvents()
        time.sleep(0.01)


    @Slot(bool, str, str)
    def receive_elog_notification(self, is_accepted, logbook_url, elog_message):
        '''Receive notification from ELOG, and report to log window'''
        if is_accepted:
            _yes_no = "made"
        else:
            _yes_no = "failed"

        _mess = "Entry into ELOG: {0} {1}".format(logbook_url, _yes_no)

        if is_accepted: 
            self.showMessage(MsgSeverity.INFO, self.pymodule, __LINE__(), _mess)
        else:
            _mess = _mess + ".\n" + elog_message 
            self.showMessage(MsgSeverity.WARN, self.pymodule, __LINE__(), _mess) 
            
        self.statusbar.showMessage(_mess)

    @Slot(dict)
    def receive_analysis_results(self, results_dict):
        '''Receive results from Analysis thread'''
        self.results_dict = results_dict
        with h5py.File(self.h5_filename, 'r') as f:
            if 'memo' in f.keys():
                memo = str(np.array(f['memo']))
            else:
                memo = ''

        important_des = f"Analysed dataset: {os.path.basename(self.h5_filename)} \n\nMemo: "+ memo
        important_dict = {}
        important_dict['number of spectra:'] = str(results_dict['number_of_spectra'])
        important_dict['pulse duration FWHM [as]:'] = str(int(results_dict['minimal_pulse_duration'][0])) + ' \u00B1 ' + str(int(results_dict['minimal_pulse_duration'][1]))
        #important_dict['average spike width FWHM [eV]:'] = str(round(results_dict['average_spike_width'][0], 2)) + ' \u00B1 ' + str(round(results_dict['average_spike_width'][1], 2))
        important_dict['average photon energy [eV]:'] = str(int(results_dict['average_photon_energy'][0])) + ' \u00B1 ' +  str(int(results_dict['average_photon_energy'][1]))

        important_dict['analysed spectra:'] = str(results_dict['number_of_spectra']-results_dict['noisy_spectra'])
        #important_dict['maximal pulse duration [as]:'] = str(int(results_dict['maximal_pulse_duration'][0])) + ' \u00B1 ' + str(int(results_dict['maximal_pulse_duration'][1]))
        important_dict['average spike number:'] = str(round(results_dict['average_spike_number'], 1))
        important_dict['intensity fluctuations [%]:'] = str(int(results_dict['intensity_fluctuations']*100))
        
        average_dict = {}
        average_dict['number of spectra:'] = str(results_dict['number_of_spectra'])
        average_dict['analysed spectra:'] = str(results_dict['number_of_spectra']-results_dict['noisy_spectra']) # + '/' + str(results_dict['number_of_spectra'])
        average_dict['average photon energy [eV]:'] = str(int(results_dict['average_photon_energy'][0]))+ ' \u00B1 ' + str(int(results_dict['average_photon_energy'][1]))
        average_dict['intensity fluctuations [%]:'] = str(int(results_dict['intensity_fluctuations']*100))
        average_dict['bandwidth FWHM [eV]:'] = str(round(results_dict['average_spectrum_width'], 2))
        
        plots = {}
        subplot = results_dict['Figure3'].get_axes()[0]
        subplot.set_xlabel('photon energy [eV]', fontsize=13)
        subplot.set_ylabel('intensity [arb. units]', fontsize=13)
        subplot.set_title('average spectrum', fontsize=13)
        subplot.tick_params(labelsize=10)
        
        subplot = results_dict['Figure4'].get_axes()[0]
        subplot.set_xlabel('number of spikes', fontsize=13)
        subplot.set_ylabel('fraction of shots [%]', fontsize=13)
        subplot.tick_params(labelsize=10)

        plots['avg spectra'] = results_dict['Figure3']
        plots['histogram'] = results_dict['Figure4']
        

        spike_summary_dict = {}
        spike_summary_dict['average spike width FWHM [eV]:'] = str(round(results_dict['average_spike_width'][0], 2)) + ' \u00B1 ' + str(round(results_dict['average_spike_width'][1], 2))
        spike_summary_dict['average spike number:'] = str(round(results_dict['average_spike_number'], 1))

        spike_summary_dict['noisy spectra:'] = str(results_dict['noisy_spectra'])
        
        spike_width_table_dict = {}
        spike_width_table_dict['peak number'] = 'average spike width FWHM [eV]'

        for i in range(len(results_dict['spike_widths'][0])):
            spike_width_table_dict[str(i + 1)] = str(round(results_dict['spike_widths'][0][i], 2)) + ' \u00B1 ' + str(round(results_dict['spike_widths'][1][i], 2))

        auto_dict = {}
        auto_dict['pulse duration FWHM [as]'] = str(round(int(results_dict['tauto'])))
        auto_dict['resolution [eV]'] = str(round(float(results_dict['resolution']), 4))
        
        if self.input['spike-width'] == 0.0:
            important_dict['analysed spectra:'] = '0'
            
        self.result_update(important_des, important_dict, plots, average_dict, spike_summary_dict, spike_width_table_dict, auto_dict)
        self.canvas_update(results_dict)
        
    @Slot(list)
    def receive_figures(self, figures_dict):
        self.canvas_nplts = len(figures_dict)
        self.figures_dict = figures_dict
        self.spectrum_range = array(range(self.canvas_nplts))
        self.current_canvas_idx = 0
        self.update_plot(0)
        
    @Slot(str, str)
    def show_analysis_error_message(self, title, message):
            QMessageBox.critical(self, title, message) 
################################################################################
#  Progagated from Slots
################################################################################
    def result_update(self, des, imp, plots, avg, sp_sum, sp_tab, auto_dict):
        
        if self.important_results_box is not None:
            self.analysis_inner_outer_hbox.removeWidget(self.important_results_box)
            self.important_results_box.deleteLater()

        self.important_results = QResultsWidget(description=des, summary_dict=imp, table_dict={}, figures_dict=plots)
        self.important_results.setBoldFace('minimal pulse duration [as]:')
        self.important_results_box = self.important_results.group_box('Results', maxlabel=3, varvalmin=125, labelwidth=300)
        self.analysis_inner_outer_hbox.addWidget(self.important_results_box)
        
        if self.average_results_box is not None:
            self.plot1_layout.removeWidget(self.average_results_box)
            self.average_results_box.deleteLater()

        self.average_results = QResultsWidget(summary_dict=avg, table_dict={})
        self.average_results_box = self.average_results.group_box('General information', varvalmin=50)
        self.average_results_box.setFixedWidth(500)
        self.average_results_box.setFixedHeight(self.plotheight)
        self.plot1_layout.addWidget(self.average_results_box)

        if self.detailed_spike_results_box is not None:
            self.plot2_layout.removeWidget(self.detailed_spike_results_box)
            self.detailed_spike_results_box.deleteLater()

        self.detailed_spike_results = QResultsWidget(summary_dict=sp_sum, table_dict=sp_tab)
        self.detailed_spike_results_box = self.detailed_spike_results.group_box('Spike statistics', varvalmin=95)
        self.detailed_spike_results_box.setFixedWidth(500)
        self.detailed_spike_results_box.setFixedHeight(self.plotheight)
        self.plot2_layout.addWidget(self.detailed_spike_results_box)

        if self.auto_correlation_results_box is not None:
            self.plot4_layout.removeWidget(self.auto_correlation_results_box)
            self.auto_correlation_results_box.deleteLater()

        self.auto_correlation_results = QResultsWidget(summary_dict=auto_dict)
        self.auto_correlation_results_box = self.auto_correlation_results.group_box('Auto-correlation', varvalmin=50)
        self.auto_correlation_results_box.setFixedWidth(self.infowidth)
        self.auto_correlation_results_box.setFixedHeight(self.plotheight)
        self.plot4_layout.addWidget(self.auto_correlation_results_box)
        
    def canvas_update(self, results_dict): 
        #Average spectrum
        #self.clearGraphsFromImage()
        
        subplot = results_dict['Figure1'].get_axes()[0]
        #subplot.set_xlabel('hi', fontsize=30)
        if self.canvas1 is not None:
            self.plot1_layout.removeWidget(self.canvas1)
            self.canvas1.deleteLater()
            
        if results_dict['Figure1'] is not None:            
            self.canvas1 = FigureCanvasQTAgg(results_dict['Figure1'])
            self.canvas1.draw()
            self.canvas1.setFixedWidth(self.plotwidth)
            self.canvas1.setFixedHeight(self.plotheight)

            self.plot1_layout.addWidget(self.canvas1)
            self.plot1_layout.setContentsMargins(10, 10, 10, 10)
            self.plot1_widget.setLayout(self.plot1_layout)   
            self.plot1_widget.show()
        else:
            self.canvas1 = None
      
        #for plotting    
        if self.canvas1 is not None:
            self.addGraphToImage(gtitle="Average",
                                 graph=self.plot1_widget)
           

        #Histograms
        if self.canvas2 is not None:
            self.plot2_layout.removeWidget(self.canvas2)
            self.canvas2.deleteLater()

        if results_dict['Figure2'] is not None:            
            self.canvas2 = FigureCanvasQTAgg(results_dict['Figure2'])
            self.canvas2.setFixedWidth(self.plotwidth)
            self.canvas2.setFixedHeight(self.plotheight)
            self.canvas2.draw()
            self.plot2_layout.addWidget(self.canvas2)
            self.plot2_widget.setLayout(self.plot2_layout)
            self.plot2_widget.show()
        else:
            self.canvas2 = None

        # for plotting
        if self.canvas2 is not None:
            self.addGraphToImage(gtitle="Histogram",
                                 graph=self.plot2_widget)
        
        #auto-correlation
        if self.canvas4 is not None:
            self.plot4_layout.removeWidget(self.canvas4)
            self.canvas4.deleteLater()

        if results_dict['Figure5'] is not None:            
            self.canvas4 = FigureCanvasQTAgg(results_dict['Figure5'])
            self.canvas4.setFixedWidth(self.plotwidth)
            self.canvas4.setFixedHeight(self.plotheight)
            self.canvas4.draw()
            self.plot4_layout.addWidget(self.canvas4)
            self.plot4_widget.setLayout(self.plot4_layout)
            self.plot4_widget.show()
        else:
            self.canvas4 = None

        # for plotting
        if self.canvas4 is not None:
            self.addGraphToImage(gtitle="Auto-correlation",
                                 graph=self.plot4_widget)
        QApplication.processEvents()    
       


################################################################################
#   Override base methods
################################################################################

    def closeEvent(self, event):
        '''Application must close cleanly!
        '''           
        if self.analysis_thread is not None:
            if self.analysis_thread.isRunning():
                qm = QMessageBox()
                qm.setText(("Measurement in progress." +
                            "Please try again in a few seconds."))
                qm.exec()
                return event.ignore()

        #Are we sure we can close the application cleanly?    
        #Say why it may not be safe to close
        if self.is_dirty:      
            qm = QMessageBox()
            _ret= qm.question(self,'', (
                "Example message: " +
                "Optics have not been restored to their initial values." + 
               "\nAre you sure you want to Quit ?"), qm.Yes | qm.No)
            if _ret == qm.No:
                return event.ignore()
          
        #Close channel access          
        self.closeBaseProc(event)


    def helpAbout(self):
        '''Overrides base class method, providing "About" application specific 
        details"
        '''
        QMessageBox.about(
            self, "About",
            """<b>{0}</b> v {1}
            <p>Copyright &copy; Paul Scherrer Institut (PSI). 
            All rights reserved.</p>
            <p>Author: A. Rutschmann </p> 
            <p>1st Responsible: Eduard Prat Costa,
            <br>Email: eduard.prat@psi.ch </p> 
            <p>2nd Responsible: Alexander Malyzhenkov, 
            <br>Email: alexander.malyzhenkov@psi.ch </p>
            <p>Application reconstructs pulse durations of FEL pulses from the photon spectra. </p>
            <p>Python {2} - Qt {3} - PyQt {4} on {5}""".format(
                 _pymodule, _appversion, platform.python_version(),
                 QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))
   

    #Optional; otherwise use base class instance    
    def sendToELOG(self):
        '''Opens menu to send to ELOG. Overides base class method to tailor 
        user needs.
        ''' 
        _message = None

        if self.analysis_thread is not None:
            if self.analysis_thread.isRunning():
                qm = QMessageBox()
                qm.setText(("Measurement in progress." +
                            "Please try again in a few seconds."))
                qm.exec()
                return


        elif not self.daq_completed:
            qm = QMessageBox()
            qm.setText(
                "Please note that no measurements have yet been undertaken")
            qm.exec()

        _attach_files = []
        _folder_name = self.settings.elogDest + self.appname + "/" 

        if not os.path.exists(_folder_name):
            os.makedirs(_folder_name)

        ''' Save any figure to _folder_name and append them to _attach_files
        '''
      
        if not _attach_files: 
            _attach_files = None

        elogSF=ElogSwissFEL()

        #Optional- Can overide default, else select from Dialog window
        _logbook = None 
        _categoryIdx = elogSF.category.MEASUREMENT
        _systemIdx = elogSF.system.BEAMDYNAMICS
    
        QSendToELOG(self, logbook=_logbook, categoryIdx=_categoryIdx,
                    systemIdx=_systemIdx, title=_elog_title,
                    message=_message, attachFile=_attach_files)
        QApplication.processEvents()

#    def loadFile(self, fname):
#        super(StartMain).loadFile(self, fname)
#        self.
    def saveashdf5(self):    
        input_options=OrderedDict()
        #QLineEdit
        input_options['Phase setting'] = None
        try:
            with h5py.File(self.h5_filename, 'r') as f:
                input_options['Date'] = str(np.array(f['date']))
        except:    
            input_options['Date'] = None
        
        try:
            with h5py.File(self.h5_filename, 'r') as f:
                input_options['Time'] = str(np.array(f['time']))
                print(input_options['Time'])
        except: 
            input_options['Time'] = None

        input_options['Key'] = None
        
        #QCombobox if list
        input_options['Stage setting'] = ['', 'twostage', 'threestage']
        input_options['Comment'] = ''

        QSaveHDF5Analysed(self, input_options=input_options)

    def screenshot(self):
        QScreenshotPlus(self, self.appname)
######################################################################
if __name__ == "__main__":

     
  
    _app = QApplication(sys.argv)
    _splash = SFMainWindow.initApp(_app, appName=_appname, delay=10)
   
    _myapp = StartMain()   
    _myapp.show()

    if _splash is not None:
        _splash.finish(_myapp)

    _app.exec_()
