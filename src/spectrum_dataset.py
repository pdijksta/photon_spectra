import h5py
import numpy as np
from .single_spectra import single_spectra
from .slices import Slice
import multiprocessing as mp
from .parallel_call import analyse_single_spectra, feedback
import os
from matplotlib.figure import Figure
from lmfit import Model
from scipy.constants import e, hbar
from uncertainties import ufloat, unumpy
import uncertainties.umath as um
#from qtpy.QtWidgets import QMessageBox

class SpectrumDataset:
    def __init__(self, input_parameters, parallel, pool=None):
        #set parameters

        self.params = input_parameters
        self.datasetname = self.params['h5_filename']
        self.parallel = parallel

        #general information
        self.f = None
        self.spectrum_number = None
        self.windowlength = None
        self.fact = None
        self.all_energy = None
        self.average_energy = None
        self.intensity_fluctuations = None
        self.all_average_photon_energy = None
        self.average_photon_energy = None
        self.average_spectra = None

        self.fit = None
        self.bandwidth = None
        self.noisebool = None

        #rawdata
        self.number_of_peaks = None
        self.FWHM = None
        self.tmin = None
        self.examples = None
        self.spectralist = None

        #processed data
        self.overall_average_FWHM = None
        self.FWHM_std = None
        self.tmin_avg = None
        self.tmax_avg = None
        self.tmin_std = None
        self.tmax_std = None

        self.FWHM_std_per_spike = None
        self.FWHM_per_spike = None

        self.average_peak_number = None
        self.analysed_spectra = None
        self.noisy_spectra = None
        self.counts = None
        self.spectra_with_peaks = None
        self.histogram = None

        self.tauto = None
        self.resolution = None
        self.deltaE = None
        self.G2 = None
        self.autoparams = []
        #setup subprocesses for parallel computation
        if self.parallel:
            if pool is not None:
                self.expool = True
                self.pool = pool
            else:
                self.pool = mp.Pool(8)
                self.expool = False

    def open_dataset(self, photEcutoff=None):
        if self.datasetname.endswith('.npz'):
            if photEcutoff is not None:
                raise ValueError('Parameter photEcutoff not supported!')
            status = self.open_dataset3(self.datasetname)
        elif self.datasetname.endswith('.h5') or self.datasetname.endswith('.hdf5'):
            try:
                with h5py.File(self.datasetname, 'r') as f:
                    status = self.open_dataset2(f, photEcutoff)

            except Exception as e:
                print(e)
                self.send_message('Wrong file format', f'The file {self.datasetname} could not be opened as an h5 file.')
                return False
        _, indices = np.unique(self.intense, axis=0, return_index=True)
        self.intense = self.intense[np.sort(indices)]
        self.intense = self.intense[np.sum(self.intense, axis=1) != 0]
        return status

    def open_dataset3(self, datasetname):
        data = np.load(self.datasetname)
        self.photE = data['e_axis']
        self.intense = (data['map']).T
        return True

    def open_dataset2(self, f, photEcutoff):
        keys = list(f.keys())
        if 'y-axis' in keys and 'x-axis' in keys:
            self.intense = f['y-axis']
            self.photE = np.array(f['x-axis'])
            return True

        success1, success2 = True, True
        for key in keys:
            if key.endswith('SPECTRUM_X'):
                photE = np.array(f[key]['data'][0], np.float64)
                success1 = True
            elif key.endswith('SPECTRUM_Y'):
                intense = np.array(f[key]['data'], np.float64)
                success2 = True
        if success1 and success2:
            self.photE = photE
            self.intense = intense
            return True

        if 'scan 1/data/SARFE10-PSSS059/SPECTRUM_Y' in f and 'scan 1/data/SARFE10-PSSS059/SPECTRUM_X' in f:
            #self.photE = np.array(f['scan 1/data/SARFE10-PSSS059/SPECTRUM_X'][0, 0,:2560])
            #self.intense = f['scan 1/data/SARFE10-PSSS059/SPECTRUM_Y'][:,0,:2560]

            xx = np.array(f['scan 1/data/SARFE10-PSSS059/SPECTRUM_X'])
            yy = np.array(f['scan 1/data/SARFE10-PSSS059/SPECTRUM_Y'])
            shape = xx.shape
            xx1 = xx.reshape([shape[0]*shape[1], shape[2]])
            yy1 = yy.reshape([shape[0]*shape[1], shape[2]])
            mask0 = xx1[0] != 0
            self.photE = xx1[:,mask0][0]
            self.intense = yy1[:,mask0]
            return True

        if 'energy' in keys and 'powerE' in keys:
            if photEcutoff is not None:
                self.photE = np.array(f['energy'][0][photEcutoff:-photEcutoff])
                self.intense = f['powerE'][:,photEcutoff:-photEcutoff]
                return True
            self.photE = np.array(f['energy'][0])
            self.intense = f['powerE']
            return True
        self.send_message('Could not find data','''No photon spectra data could be found in the opened h5 file. We were looking for data under the following keys <ol> <li> x-axis, y-axis </li> <li> SARFE10-PSSS059:SPECTRUM_X, SARFE10-PSSS059:SPECTRUM_Y </li><li>scan 1/data/SARFE10-PSSS059/SPECTRUM_X, scan 1/data/SARFE10-PSSS059/SPECTRUM_Y</li> </ol>''')

        return False

    def get_general_information(self):

        #general information
        #catch some cases of wrong input data
        if len(self.intense.shape) != 2:
            self.send_message('Wrong data shape', f'''We expected that the intensity data has two axes, but we got {len(self.intense.shape)} axes instead.''')
            return

        elif len(self.photE.shape) != 1:
            self.send_message('Wrong data shape', f'''We expected that the photon energy data has one axis, but we got {len(self.photE.shape)} axes instead.''')
            return

        elif len(self.photE) != len(self.intense[0,:]):
            self.send_message('Wrong data shape', f'''We expected that the photon energy data and the intensity data have the same length, but we got {len(self.photE)} data points for the photon energy and {len(self.intense[0, :])} data points for the intensity per spectrum.''')
            return

        else:
            self.spectrum_number = self.intense.shape[0]
            self.windowlength = len(self.photE)
            self.fact = len(self.photE) / (np.max(self.photE) - np.min(self.photE))

        #integrated intensity and intensity fluctuations
        minimas = np.quantile(self.intense, 0.2, axis=1)
        self.all_energy = np.sum(self.intense-minimas[:,np.newaxis], axis=1)/self.fact
        self.average_energy = np.average(self.all_energy)
        self.intensity_fluctuations = np.std(self.all_energy)/self.average_energy

        #average photon energy
        self.all_average_photon_energy = np.sum(self.photE * (self.intense-minimas.reshape(-1, 1)), axis=1) / self.fact/self.all_energy
        self.average_photon_energy = np.average(self.all_average_photon_energy)
        self.std_photon_energy = np.std(self.all_average_photon_energy)

        #average spectra
        self.average_spectra = np.average(self.intense, axis=0)
        self.average_spectra = self.average_spectra - np.min(self.average_spectra)
        self.average_spectra /= self.average_energy

        #determine the cutoffs for average spectrum fit
        maxind = np.argmax(self.average_spectra)

        #catch the case where the maximum is at the edge of the spectra to prevent min1 = min2
        if maxind != 0 and maxind != self.windowlength-1:
            min1 = np.argmin(self.average_spectra[:maxind])
            min2 = np.argmin(self.average_spectra[maxind:])+maxind
        else:
            min1 = 0
            min2 = self.windowlength

        #average spectra fit
        self.fit = Slice(self.average_spectra, min1, min2)
        self.bandwidth = self.fit.sigma*2.355/self.fact

        #sort out spectra with integrated intensity below treshold
        self.noisebool = (self.all_energy > self.params['intensity_thresh'] * self.all_energy.max())

        self.spectralist = self.spectrum_number*[1]

        #import pdb; pdb.set_trace()
        return True

    def analyse(self):
        if self.parallel:
            bo = self.analyse_parallel()
        else:
            bo = self.analyse_serial()

        return bo

    def analyse_serial(self):
        #initialize rawdata containers
        self.number_of_peaks = np.zeros(self.spectrum_number)
        self.FWHM = self.number_of_peaks.copy()
        self.tmin = self.number_of_peaks.copy()
        self.fit_functions = np.zeros([self.spectrum_number, len(self.photE)])
        self.filtered_spectra = self.fit_functions.copy()
        self.examples = []

        #feedbackloop
        if self.params['parameter optimization']:

            self.params['roughness'] = 0
            counter1 = 0
            total_weight = 0
            lpcutoffs = []
            alphas = np.array([5, 20, 50, 100, 200, 300]) #test values for the parameter optimization

            while counter1 < min(self.spectrum_number, 10):

                #find best alpha for spectrum counter1
                spectra = single_spectra(self.intense[counter1, :] / self.average_energy, self.photE, self.fact, counter1)

                noise, lpcutoff, alpha, weight = spectra.findbestalpha(alphas, self.params)
                lpcutoffs.append(lpcutoff) #save the cutoffs for later analysis
                total_weight += weight
                self.params['roughness'] += weight*alpha
                counter1 += 1

            #catch case where the optimization isn't conclusive
            if total_weight == 0:
                self.params['roughness'] = 250
            else:
                self.params['roughness'] /= total_weight

            #determine the corresponding cutoff for the optimized alpha value
            i = np.argmin(np.abs(alphas-self.params['roughness']))
            lpcutoffs = np.array(lpcutoffs)
            self.params['frequency_cutoff'] = np.average(lpcutoffs[:, i][lpcutoffs[:, i] != 0])

        count1 = 0

        for i in range(self.spectrum_number):
            #send signal to progressbar
            if self.send_progress(i/self.spectrum_number*100):
                return

            #noisebool[i] is True if the spectrum lies above intensity threshold
            if self.noisebool[i]:

                #lowpassfilter
                spectra = single_spectra(self.intense[i, :]/self.average_energy, self.photE, self.fact, i)
                lpcutoff = spectra.adaptive_butterworth_filter(self.params['frequency_cutoff'], self.params['snr'], self.params['sbr'])
                #check if the spectrum is noisy or not
                if not(spectra.noisybool):

                    #finds peaks
                    spectra.finding_peak(self.params['prominence'], self.params['height'], self.params['background_prominence'])
                    #slices and fits the spectrum
                    spectra.first_estimate(self.params['cutoff_prominence'], self.params['cutoff_height'], globalminima=False)

                    #gather results
                    if spectra.n >= 1:
                        self.number_of_peaks[i] = spectra.n
                        self.FWHM[i] = spectra.FWHM
                        self.tmin[i] = spectra.tmin
                        self.fit_functions[i] = spectra.total_fit
                        self.filtered_spectra[i] = spectra.spec

                        if count1 <= 1:
                            self.examples.append(spectra.shiftraw)
                            self.examples.append(spectra.total_fit)
                            count1 += 1
                else:
                    #add noisy spectra to the noisebool list
                    self.noisebool[i] = False
                self.spectralist[i] = spectra
        self.data_analysis()

        return True

    def analyse_parallel(self):

        #initialize values
        self.examples = []

        #feedbackloop
        if self.params['parameter optimization']:

            lpcutoffs = []
            self.params['roughness'] = 0
            total_weight = 0

            #calls findbestalpha via feedback function from parallel_call file
            result = self.pool.starmap(feedback, [(i, self.params, self.intense[i, :]/self.average_energy, self.photE, self.fact) for i in range(min(self.spectrum_number, 10))])

            #evaluate results
            for spectra in result:
                total_weight += spectra.weight
                self.params['roughness'] += spectra.weight * spectra.alpha
                lpcutoffs.append(spectra.lpcutoff)

            #catch case where the feedback isn't conclusive
            if total_weight == 0:
                self.params['roughness'] = 250
            else:
                self.params['roughness'] /= total_weight

            #find corresponding cutoff for optimized alpha
            alphas = np.array([5, 20, 50, 100, 200, 300])
            i = np.argmin(np.abs(alphas - self.params['roughness']))
            lpcutoffs = np.array(lpcutoffs)
            self.params['frequency_cutoff'] = np.average(lpcutoffs[:, i][lpcutoffs[:, i] != 0])

        #send to progressbar return if analysis is aborted
        if self.send_progress(30):
            return
        #calls analyse_single_spectra from parallel_call file for parallel spectrum analysis
        resulttmp = self.pool.starmap_async(analyse_single_spectra, [(i, self.params, self.intense[i, :]/self.average_energy, self.photE, self.fact) for i in np.arange(0, self.spectrum_number)[self.noisebool]])

        #if the analysis takes longer than three minutes abort it, one could instead implement a system where the analysis can be aborted by the abort button
        try:
            result = resulttmp.get(timeout=180)
        except Exception as e:
            print(e)
            self.pool.terminate()
            self.pool.join()
            print('timeout')
            self.send_message('Timeout', '''The analysis took over 180 seconds and was terminated.''')
            return False

        if not self.expool:
            self.pool.close()
            self.pool.join()
        count1 = 0
        self.number_of_peaks = np.zeros(self.spectrum_number)
        self.FWHM = self.number_of_peaks.copy()
        self.tmin = self.number_of_peaks.copy()
        self.fit_functions = np.zeros([self.spectrum_number, len(self.photE)])
        self.filtered_spectra = self.fit_functions.copy()

        #readout analysed data
        self.all_spike_widths = {}
        self.all_spike_amplitudes = {}
        self.all_spike_centers = {}
        for num, spectra in zip(np.arange(0, self.spectrum_number)[self.noisebool], result):
            self.spectralist[spectra.spectrum_label] = spectra
            if spectra.noisybool:
                self.noisebool[spectra.spectrum_label] = False
            else:
                #print(num, spectra.n)
                if spectra.n >= 1:

                    self.number_of_peaks[spectra.spectrum_label] = spectra.n
                    self.FWHM[spectra.spectrum_label] = spectra.FWHM
                    self.tmin[spectra.spectrum_label] = spectra.tmin
                    self.fit_functions[spectra.spectrum_label] = spectra.total_fit
                    self.filtered_spectra[spectra.spectrum_label] = spectra.spec
                    self.all_spike_widths[str(num)] = np.array(spectra.all_fwhm)
                    self.all_spike_amplitudes[str(num)] = np.array(spectra.all_amplitude)
                    self.all_spike_centers[str(num)] = np.array(spectra.all_center)

                    if count1 <= 1:
                        self.examples.append(spectra.shiftraw)
                        self.examples.append(spectra.total_fit)
                        count1 += 1
        self.data_analysis()

        return True

    def data_analysis(self):

        #get indices of non-noisy spectra
        self.noisy_spectra = len(np.where(self.noisebool == 0)[0])
        self.analysed_spectra = self.spectrum_number - self.noisy_spectra

        #get the number of spectra with peaks
        self.spectra_with_peaks = len(np.where(self.number_of_peaks != 0)[0])

        if self.spectra_with_peaks > 0:

            #calculate average values of pulse duration, spike number and spike width
            maxpeak = np.max(self.number_of_peaks)
            self.FWHM_per_spike = np.zeros(int(maxpeak))
            self.counts = np.zeros(int(maxpeak))
            self.FWHM_std_per_spike = np.zeros(int(maxpeak))

            #sort data by total spike number
            for i in range(self.spectrum_number):
                if int(self.number_of_peaks[i]) > 0:
                    index = int(self.number_of_peaks[i]) - 1
                    self.counts[index] += 1
                    self.FWHM_per_spike[index] += self.FWHM[i]

            #calculate average
            A = self.counts != 0 #exclude the numbers, where there are no spectrum with this number
            self.FWHM_per_spike[np.logical_not(A)] = 0
            self.FWHM_per_spike[A] = self.FWHM_per_spike[A] / self.counts[A]

            #calculate the standard deviation of the spike width for each number of spikes
            for i in range(self.spectrum_number):
                if int(self.number_of_peaks[i]) > 0:
                    index = int(self.number_of_peaks[i]) - 1
                    self.FWHM_std_per_spike[index] += (self.FWHM[i]-self.FWHM_per_spike[index])**2

            #calculate averages and standard deviations
            self.FWHM_std_per_spike[A] = np.sqrt(self.FWHM_std_per_spike[A])/(np.sqrt(self.counts[A]))

            self.overall_average_FWHM = np.sum(self.FWHM_per_spike * self.counts) / self.spectra_with_peaks
            self.average_peak_number = np.sum(self.counts * np.arange(1, int(np.max(self.number_of_peaks) + 1))) / self.spectra_with_peaks

            average2 = np.average(self.number_of_peaks[np.where(self.number_of_peaks != 0)[0]])
            assert average2 == self.average_peak_number
            self.average_peak_std = np.std(self.number_of_peaks[np.where(self.number_of_peaks != 0)[0]])

            self.FWHM_std = np.std(self.FWHM[np.where(self.FWHM != 0)[0]])

            def rational(x, q0, q1, q2, p0, p1, p2):
                q = [q0, q1, q2]
                p = [p0, p1, p2]
                return np.polyval(q, x)/np.polyval(p, x)

            self.tmin_avg = np.average(self.tmin[np.where(self.tmin != 0)[0]])
            self.tmin_std = np.std(self.tmin[np.where(self.tmin != 0)[0]])

            if self.params['pulse duration correction']:
                with h5py.File(os.path.dirname(__file__)+'/correction_coefficient.h5', 'r') as f:
                    fact = rational(self.average_peak_number, *np.array(f['coef']))
                    print('Correction factor', fact)
                    self.tmin_avg /= fact
                    self.tmin_std /= fact

            self.tmax_avg = np.sqrt(2)*self.tmin_avg
            self.tmax_std = np.sqrt(2)*self.tmin_std

            #calculate number for histogram
            self.histogram = self.counts / self.spectra_with_peaks * 100

    def second_order_correlation(self, number=None):

        if number is not None:
            shiftintense = self.intense[:number,:] - np.min(self.intense[:number,:], axis=1).reshape(-1,1)
            shiftintense = shiftintense.astype(float)
            shiftintense /= self.all_energy[:number].reshape(-1,1)
        else:
            shiftintense = self.intense - np.min(self.intense, axis=1).reshape(-1,1)
            shiftintense = shiftintense.astype(float)
            shiftintense /= self.all_energy.reshape(-1,1)
        avgspectralpower = np.average(shiftintense, axis=0)
        avgborder = np.zeros(3*self.windowlength)
        avgborder[self.windowlength:2*self.windowlength] = avgspectralpower*1.0

        try:
            assert len(avgspectralpower) == self.windowlength
        except AssertionError:
            print(len(avgspectralpower))

        corint1 = np.zeros(self.windowlength)

        for i in range(self.windowlength):
            interval = min(i, self.windowlength-1-i)
            corint1[i] = np.sum(avgspectralpower[i-interval:i+1]*avgspectralpower[i:i+interval+1])/self.fact

        corint2 = np.sum(corint1)/self.fact

        weight = corint1/corint2

        deltaEmax = int(min(self.bandwidth*self.fact, self.windowlength-1))
        deltaE = np.arange(-deltaEmax*2, deltaEmax*2+1, 2, dtype=float)
        G2 = np.zeros_like(deltaE, dtype=float)

        secordercor = np.zeros((self.windowlength, deltaEmax+1))

        for i in range(self.windowlength):
            for j in range(deltaEmax+1):
                if i-j>= 0 and i + j < self.windowlength:
                    secordercor[i, j] = np.average(shiftintense[:, i-j]*shiftintense[:, i+j])
                else:
                    secordercor[i, j] = 0

        for i in range(deltaEmax+1):
            A = secordercor[:, i] != 0
            G2[deltaEmax+i] = np.sum(weight[A]*secordercor[A, i]/(avgborder[self.windowlength-i:self.windowlength*2-i][A]*avgborder[self.windowlength+i:self.windowlength*2+i][A]))/self.fact
            G2[deltaEmax-i] = G2[deltaEmax+i]*1.0

        G2 -= 1
        deltaE /= self.average_photon_energy*self.fact

        self.G2 = G2
        self.deltaE = deltaE
        goffsetmodel = Model(gaussian_central_offset)
        result = goffsetmodel.fit(data=G2, x=deltaE, amplitude=0.5, sigma=10**-4, offset=-0.1)


        amp = ufloat(result.params['amplitude'].value, result.params['amplitude'].stderr)
        sigma = ufloat(result.params['sigma'].value, result.params['sigma'].stderr)
        for key in result.params.keys():
            self.autoparams.append(result.params[key])
        if unumpy.nominal_values(amp) >= 1:
            #t_g = ufloat(0,0)
            self.tauto = ufloat(0,0)
            self.resolution = ufloat(0,0)
        else:
            #t_g = 0.5*um.sqrt(1/amp**2-1)
            self.tauto = 1/(np.sqrt(2)*amp*sigma*self.average_photon_energy*e/hbar)*2.355*10**18
            self.resolution = um.sqrt((1-amp**2)/(4*amp**2*(self.tauto/(2.355*10**18))**2))/e*hbar

#        if number is not None:
#            num = number
#        else:
#            num = self.spectrum_number

#        with h5py.File('second_order_correlation/Genesis_analysed/'+os.path.basename(self.datasetname)[:-3]+'_n='+str(num)+'.h5', 'w') as f:
#            f.create_dataset('dE', data=deltaE)
#            f.create_dataset('G2', data=G2)
#            pd = np.array([unumpy.nominal_values(pulse_duration),unumpy.std_devs(pulse_duration)])
#            f.create_dataset('pulse duration', data=pd)
#            f.create_dataset('spectra number', data=num)
#            f.create_dataset('photon energy', data=np.array([self.average_photon_energy, self.std_photon_energy]))
#            tg = np.array([unumpy.nominal_values(t_g), unumpy.std_devs(t_g)])
#            f.create_dataset('tg', data=tg)
#            relrestmp = np.array([unumpy.nominal_values(relres), unumpy.std_devs(relres)])
#            f.create_dataset('relres', data=relrestmp)

    def get_single_spectra(self, i):

        fig = Figure()
        axes = fig.add_subplot(111)
        self.plot_single_spectra(i, axes)
        fig.set_tight_layout(True)

        return fig

    def plot_single_spectra(self, i, axes):

        axes.clear()
        axes.plot(self.photE, self.intense[i, :]/self.average_energy, color='black')
        axes.grid(True)
        axes.set_xlabel('photon energy [eV]')
        axes.set_ylabel('intensity')
        axes.set_title('spectrum '+str(i))

    def get_average_figure(self):
        fig = Figure()
        axes = fig.add_subplot(111)
        self.plot_average_spectra(axes)
        fig.set_tight_layout(True)

        return fig

    def plot_average_spectra(self, axes):

        axes.clear()
        axes.plot(self.photE, self.average_spectra, label="Average spectrum", color='black')
        axes.plot(self.photE, self.fit.gausscurve, label="Average spectrum", color='red')
        axes.set_xlabel('photon energy [eV]')
        axes.set_ylabel('intensity [arb. units]')
        axes.grid(True)

    def get_histogram(self):

        fig = Figure()
        axes = fig.add_subplot(111)
        self.plot_histogram(axes)
        fig.set_tight_layout(True)

        return fig

    def plot_histogram(self, axes):
        if self.counts is not None:
            list = np.arange(1, len(self.counts)+2)
            axes.clear()
            axes.bar(list[:-1], self.counts/self.spectra_with_peaks*100, np.diff(list), edgecolor='black')
        axes.set_xlabel('number of spikes')
        axes.set_ylabel('fraction of shots [%]')
        axes.set_title('Spike histogram')

    def get_auto_correlation_plot(self):

        fig = Figure()
        axes = fig.add_subplot(111)
        self.plot_auto_correlation_function(axes)
        fig.set_tight_layout(True)
        return fig

    def plot_auto_correlation_function(self, axes):
        if self.deltaE is not None and self.G2 is not None:
            axes.clear()
            axes.scatter(self.deltaE, self.G2, color='black')
            x = np.linspace(self.deltaE[0], self.deltaE[-1], 100)
            axes.plot(x, gaussian_central_offset(x, *self.autoparams), color='red')
            axes.set_xlim([self.deltaE[0], self.deltaE[-1]])
            axes.grid(True)
        axes.set_xlabel('relative photon energy width')
        axes.set_ylabel('second order correlation $G_2$')

    def get_figure_list(self):
        #get a list of matplotlib figures of the spectra
        spectra_plot_list = []
        N = min(self.spectrum_number, 100) #maximal number is 100

        for i, sp in enumerate(self.spectralist[:N]):
            if sp == 1: #below intensity treshold, gives raw spectra
                spectra_plot_list.append(self.get_single_spectra(i))

            elif self.noisebool[i]: #too noisy, gives raw spectra + filtered signal
                spectra_plot_list.append(sp.get_gaussian_plot())

            else: #analysed spectrum, gives raw spectra + fit
                spectra_plot_list.append(sp.get_lowpass_signal_plot())

        return spectra_plot_list

    def saveh5py(self, filename, stage='', date='', phase='', key='', comment='', time=''):
        dictkeys = ['number_of_spectra', 'average_photon_energy', 'photon_energy', 'average_spectrum',
                    'average_spike_width', 'minimal_pulse_duration', 'maximal_pulse_duration', 'spike_widths',
                    'average_spike_number', 'noisy_spectra', 'peak_histogram', 'average_spectrum_width', 'pulse_energy',
                    'relative_energy']

        dictvalues = [self.spectrum_number, self.average_photon_energy, self.photE,
                      [self.average_spectra, self.fit.gausscurve],
                      [self.overall_average_FWHM, self.FWHM_std],
                      [self.tmin_avg, self.tmin_std], [self.tmax_avg, self.tmax_std],
                      [self.FWHM_per_spike, self.FWHM_std_per_spike],
                      [self.average_peak_number, self.average_peak_std], self.noisy_spectra,
                      self.histogram, self.bandwidth, 0,
                      self.all_energy/self.average_energy]

        datadict = dict(zip(dictkeys, dictvalues))

        newfile = h5py.File(filename, 'w')
        newfile.create_dataset('metadata/phase_setting', data=phase)
        newfile.create_dataset('metadata/stage_setting', data=stage)
        newfile.create_dataset('metadata/date', data=date)
        newfile.create_dataset('metadata/time', data=time)
        newfile.create_dataset('metadata/key', data=key)
        newfile.create_dataset('metadata/comment', data=comment)

        for key in datadict:
            newfile.create_dataset('data/' + key, data=np.array(datadict[key]))

        newfile.create_dataset('examples', data=np.array(self.examples))
        newfile.close()

    def send_progress(self, percent):
        #override for implementation in GUI
        return False

    def send_message(self, title, message):
        #override for implementation in GUI
        print(title, message)
        pass

def gaussian_central_offset(x, amplitude, sigma, offset):
    return amplitude* np.exp(-x**2 / (2*sigma**2)) + offset

class GUISpectrumDataset(SpectrumDataset):

    def __init__(self, parent, input_parameters, parallel, pool=None):
        super(GUISpectrumDataset, self).__init__(input_parameters, parallel, pool=pool)
        self.parent = parent

    def send_progress(self, percent):
        if percent > 10:
            if self.parent.parent is not None:
                self.parent.parent.trigger_progressbar.emit(int(percent))
        return self.parent.abort

    def send_message(self, title, message):
        if self.parent.parent is not None:
            self.parent.parent.trigger_analysis_message.emit(title, message)
        else:
            print(title, message)

