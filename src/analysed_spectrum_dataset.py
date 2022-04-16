import h5py
import numpy as np
from matplotlib.figure import Figure 
import matplotlib.pyplot as plt

class analysed_spectrum_dataset:

    def __init__(self, filename):
        self.filename = filename
        self.categories = ['metadata', 'data', 'examples']
        self.elements = {}
        self.elements['metadata'] = ['stage_setting', 'phase_setting', 'date', 'time', 'key', 'comment']
        self.elements['data'] = ['number_of_spectra', 'average_photon_energy', 'photon_energy', 'average_spectrum',
                    'average_spike_width', 'minimal_pulse_duration', 'maximal_pulse_duration', 'spike_widths',
                    'average_spike_number', 'noisy_spectra', 'peak_histogram', 'average_spectrum_width', 'pulse_energy',
                    'relative_energy']
        self.elements['examples'] = []
        
    def check_h5(self):
        try:
            file = h5py.File(self.filename, 'r')
        except:
            return False
            file.close()
        return True
        file.close()
        
    def check_format(self):
        file = h5py.File(self.filename, 'r')
        missinglist = []
        for key, items in self.elements.items():
            if not key in file:
                missinglist.append(key)
            else:
                for item in items:
                    if not key + '/' + item in file:
                        missinglist.append(item)
        return len(missinglist) == 0, missinglist
        file.close()
        
    def plot_histogram(self, a):
        file = h5py.File(self.filename, 'r')
        list = np.arange(1, len(file['data/peak_histogram']) + 2)
        a.clear()
        a.bar(list[:-1], file['data/peak_histogram'], np.diff(list), edgecolor='black')
        a.set_xlabel('number of spikes')
        a.set_ylabel('fraction of shots (%)')
        file.close()
        
    def get_histogram(self):
        file = h5py.File(self.filename, 'r')
        fig = Figure()
        subplot = fig.add_subplot(111)
        
        list = np.arange(1, len(file['data/peak_histogram']) + 2)
        subplot.bar(list[:-1], file['data/peak_histogram'], np.diff(list), edgecolor='black')
        subplot.set_xlabel('number of spikes')
        subplot.set_ylabel('fraction of shots (%)')
        fig.set_tight_layout(True)
        file.close()
        return fig
        
    def plot_average_spectra(self, a):
        file = h5py.File(self.filename, 'r')
        a.clear()
        a.plot(file['data/photon_energy'], file['data/average_spectrum'][0], color='black')
        a.plot(file['data/photon_energy'], file['data/average_spectrum'][1], color='red')
        a.set_xlabel('photon energy [eV]')
        a.set_ylabel('intensity')
        a.grid(True)
        file.close()

    def plot_examples(self, a, i):
        file = h5py.File(self.filename, 'r')
        a.clear()
        a.plot(file['data/photon_energy'], file['examples'][2*i], color='black')
        a.plot(file['data/photon_energy'], file['examples'][2*i+1], color='red')
        a.set_xlabel('photon energy [eV]')
        a.set_ylabel('intensity')
        a.grid(True)
        file.close()

    def plot_single_value(self, a, xkey, ykey, fmt='b.', color='blue'):            
        if isinstance(xkey,int):
            try:
                a.errorbar(xkey, self.__getitem__(ykey)[0], yerr=self.__getitem__(ykey)[1], capsize=2, fmt=fmt)
                return self.__getitem__(ykey)[0]+self.__getitem__(ykey)[1]
            except:
                a.plot(xkey, self.__getitem__(ykey), color=color, marker='.')
                return self.__getitem__(ykey)
        else:
            try:
                a.errorbar(self.__getitem__(xkey), self.__getitem__(ykey)[0], yerr=self.__getitem__(ykey)[1], capsize=2, fmt=fmt)
                return self.__getitem__(ykey)[0]+self.__getitem__(ykey)[1]
            except:
                a.plot(self.__getitem__(xkey), self.__getitem__(ykey), color=color, marker='.')
                return self.__getitem__(ykey)

    def __getitem__(self, key):
        file = h5py.File(self.filename, 'r')

        if key in ['stage_setting', 'date', 'time', 'key', 'comment']:
            return str(np.array(file['metadata/'+key]))
        elif key == 'phase_setting':
            if np.array(file['metadata/'+key]) == '':
                return ''
            else:
                return int(np.array(file['metadata/'+key]))
        elif key == 'examples':
            return np.array(file['examples'])
        elif key in file['data'].keys():
            return np.array(file['data/'+key])
        else: 
            return None
        file.close()
