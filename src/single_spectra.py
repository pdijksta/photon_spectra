import numpy as np
from numpy import array, std, zeros_like, zeros, logical_and, empty, where, all
from numpy import sum as npsum
from numpy import max as npmax
from numpy import min as npmin
from numpy import sqrt as npsqrt
from numpy import average as npaverage
import time
from scipy.signal import find_peaks, butter, filtfilt
from .slices import Slice

from matplotlib.figure import Figure

class single_spectra:
    def __init__(self, intensity, photE, fact, i):
        self.spectrum_label = i
        self.intensity = intensity
        self.indivtime = zeros(4)
        self.photE = photE
        self.fact = fact

        #Lowpassfilter
        self.deg = 5
        self.spec = None
        self.shiftraw = None
        self.filtered_intensity = None
        self.r2filter = None
        self.max_intensity = None
        self.noise = None
        self.noise_mag = None
        self.background = None
        self.noisybool = None

        #Peakfinder and fit
        self.xpeaks = None
        self.ypeaks = None
        self.oldpeaks = None
        self.n = None

        self.minima_indices = None
        self.minima = None
        self.leftcutoff = None
        self.rightcutoff = None
        self.Slices = None

        self.average_sigma = None
        self.total_amplitude = None
        self.total_fit = None
        self.FWHM = None
        self.tmin = None
        self.redchi = None

    def adaptive_butterworth_filter(self, lowpasscutoff, alpha, snr, sbr):
        tlowstart = time.time()
        lpcutoffinit = lowpasscutoff*1.0
        lpcutoff = lowpasscutoff*1.0
        rd = np.inf
        switched = False
        while True:

            if lpcutoff >= 1 or lpcutoff <= 0:
                b, a = butter(self.deg, lpcutoffinit, 'low')
                self.spec = filtfilt(b, a, self.intensity)
            else:
                b, a = butter(self.deg, lpcutoff, 'low')
                #print('r2filter:', rd, 'cutoff:', lpcutoff)
                self.spec = filtfilt(b, a, self.intensity)

            r2 = (npsum((self.intensity - self.spec) ** 2)+ alpha * npsum((self.spec[:-1] - self.spec[1:]) ** 2))

            old_rd = rd
            rd = npsqrt(r2 / len(self.photE))
            if rd > old_rd:
                if switched:
                    lpcutoff += 0.001
                    b, a = butter(self.deg, lpcutoff, 'low')
                    self.spec = filtfilt(b, a, self.intensity)
                    break
                else:
                    lpcutoff -= 0.001
                    switched = True
            else:
                if switched:
                    lpcutoff -= 0.001
                else:
                    lpcutoff += 0.001

        #if np.any(self.spec):
        #    print('Filtered data contains something')
        #    #import pdb; pdb.set_trace()
        #else:
        #    print('Filtered data contains nothing')





        self.r2filter = rd
        self.filtered_intensity = self.spec -npmin(self.spec)
        self.shiftraw = self.intensity - npmin(self.spec)

        #if np.any(self.shiftraw):
        #    print('shiftraw contains something')
        #else:
        #    print('shiftraw contains nothing')

        self.average_intensity = npaverage(self.shiftraw)

        self.max_intensity = npmax(self.shiftraw)

        self.noise = self.intensity - self.spec

        A = (self.shiftraw < npmax(self.shiftraw)*0.2)
        self.noise_mag = std(self.noise[A])
        self.background = 0.5*(npaverage(self.shiftraw[:10]) + npaverage(self.shiftraw[-10:]))


        self.noisybool = False

        if self.noise_mag > snr*self.max_intensity or self.background > sbr*self.max_intensity:
            self.noisybool = True


        tlowend = time.time()
        self.indivtime[0] = tlowend - tlowstart

        return lpcutoff

    def findbestalpha(self, alphas, params):

        N = len(alphas)

        redchi = zeros(N)
        peaknumbers = zeros(N)
        lpcutoffs = zeros(N)

        for i, alpha in enumerate(alphas):

            lpcutoffs[i] = self.adaptive_butterworth_filter(params['frequency_cutoff'], alpha, params['snr'], params['sbr'])

            if self.noisybool:
                self.alpha = alphas[-1]
                self.weight = 3
                lpcutoffs[-1] = self.adaptive_butterworth_filter(params['frequency_cutoff'], alphas[-1], params['snr'], params['sbr'])
                self.lpcutoff = lpcutoffs
                return True, lpcutoffs, self.alpha, self.weight
            else:
                self.finding_peak(params['prominence'], params['height'], 0)

                if self.n >= 1:
                    if i >= 1:
                        self.update_slices()
                    else:
                        self.first_estimate(0, 100, globalminima=True)
                    redchi[i] = self.redchi
                    peaknumbers[i] = self.n
                else:
                    self.alpha = alphas[-1]
                    self.weight = 3
                    lpcutoffs[-1] = self.adaptive_butterworth_filter(params['frequency_cutoff'], alphas[-1], params['snr'], params['sbr'])
                    self.lpcutoff = lpcutoffs
                    return True, lpcutoffs, self.alpha, self.weight

        #print(redchi, peaknumbers)

        redchi1 = npmin(redchi)

        minima = where(redchi < redchi1*1.3)[0]
        self.alpha = npaverage(alphas[minima])
        self.weight = N - len(minima)
        self.lpcutoff = lpcutoffs

        return False, lpcutoffs, self.alpha, self.weight

    def finding_peak(self, prominence, absheight, bgprominence):
        tpeakstart = time.time()

        if self.xpeaks is not None:
            self.oldpeaks = self.xpeaks*1
        minheight = max(absheight*npmax(self.filtered_intensity)/100, bgprominence*self.background)
        #print(minheight, prominence, print(npmax(self.filtered_intensity)))
        #print(find_peaks(self.filtered_intensity, height=minheight, prominence=prominence*self.max_intensity/100))
        self.xpeaks = find_peaks(self.filtered_intensity, height=minheight, prominence=prominence*self.max_intensity/100)[0]

        self.n = len(self.xpeaks)
        self.filterpeak()
        self.ypeaks = self.filtered_intensity[self.xpeaks]

        tpeakend = time.time()
        self.indivtime[1] = tpeakend - tpeakstart

    def filterpeak(self):
        A = logical_and(self.xpeaks > len(self.photE)/25, self.xpeaks < len(self.photE)*(1-1/25))
        self.xpeaks = self.xpeaks[A]
        self.n = len(self.xpeaks)

    def first_estimate(self,cutoff_prom, cutoff_height, globalminima=False):
        tslicestart = time.time()

        self.total_fit = zeros_like(self.filtered_intensity)
        self.minima_indices = empty(self.n+1, dtype=int)
        self.Slices = [1]*self.n
        self.average_sigma = 0
        self.total_amplitude = 0

        if self.n > 0:
            self.getnumber = True
            peak1 = 0
            for j in range(self.n):

                peak2 = self.xpeaks[j]

                if j == 0:

                    a = where(self.filtered_intensity[:peak2] < self.background*1.5)[0]
                    if len(a) > 0:
                        leftborder = a[-1]
                    else:
                        #print(self.spectrum_label, np.where(self.filtered_intensity[:peak2] < self.background*1.5), np.where(self.filtered_intensity[:peak2] < self.background*1.5)[0])
                        leftborder = 0

                    if globalminima:
                        Min_y= npmin(self.filtered_intensity[:peak2])
                        Min_x = where(self.filtered_intensity[:peak2] == Min_y)[0][-1] + peak1
                        Min_x = int(max(Min_x, self.xpeaks[j] / 2))
                    else:
                        try:
                            prom = cutoff_prom*self.max_intensity/100
                            minima = find_peaks(-self.filtered_intensity[leftborder:peak2], height=-cutoff_height/100 * self.ypeaks[0], prominence=prom)[0]
                            Min_x = minima[-1] + leftborder

                        except IndexError:
                            Min_x = 0

                        Min_x = int(max(Min_x, self.xpeaks[j] / 2, leftborder))

                    self.minima_indices[0] = Min_x*1.0
                else:
                    Min_y = npmin(self.filtered_intensity[peak1:peak2])
                    Min_x = where(self.filtered_intensity[peak1:peak2] == Min_y)[0] + peak1
                    self.minima_indices[j] = int(npaverage(Min_x))
                peak1 = peak2

            peak2old = len(self.photE)

            a = where(self.filtered_intensity[peak1:] > self.background*1.5)[0]
            if len(a) > 0:
                peak2 = int(a[-1]) +peak1
            else:
                peak2 = peak2old*1

            if globalminima:
                Min_y = npmin(self.filtered_intensity[peak1:peak2])
                Min_x = where(self.filtered_intensity[peak1:peak2] == Min_y)[0][0] + peak1

            else:

                try:
                    prom = cutoff_prom * self.max_intensity / 100 # promfactor*self.noise_mag*0.5)

                    minima = find_peaks(-self.filtered_intensity[peak1:peak2], height=-cutoff_height/100 * self.ypeaks[-1], prominence=prom)[0]
                    Min_x = minima[0] + peak1

                except IndexError:
                    Min_x = peak2old

                Min_x = int(min(Min_x, self.xpeaks[self.n - 1] + (len(self.filtered_intensity) - self.xpeaks[self.n - 1]) / 2, peak2))

            self.minima_indices[self.n] = Min_x*1.0

            self.minima = self.minima_indices[1:-1]
            self.leftcutoff = self.minima_indices[0]
            self.rightcutoff = self.minima_indices[-1]
            for k in range(self.n):

                index1 = self.minima_indices[k]
                index2 = self.minima_indices[k + 1]

                self.Slices[k] = Slice(self.shiftraw, int(index1), int(index2), model='gaussian')

                self.average_sigma += self.Slices[k].sigma*self.Slices[k].amplitude
                self.total_amplitude += self.Slices[k].amplitude
                self.total_fit += self.Slices[k].gausscurve

            self.average_sigma = self.average_sigma/self.total_amplitude

            averagescaled = self.average_sigma / self.fact

            self.FWHM = array(2.355 * averagescaled)
            self.tmin = 4 * 6.582 * 10 ** (-16) * 0.69314718056 * 10 ** (18) / self.FWHM

            self.redchi = npsum((self.shiftraw[self.leftcutoff:self.rightcutoff] - self.total_fit[self.leftcutoff:self.rightcutoff]) ** 2)/(len(self.intensity[self.leftcutoff:self.rightcutoff]))/self.average_intensity**2
        tsliceend = time.time()
        self.indivtime[2] = tsliceend - tslicestart

    def update_slices(self):
        if len(self.xpeaks) == len(self.oldpeaks) and all(self.xpeaks < self.minima_indices[1:]) and all(self.xpeaks > self.minima_indices[:-1]):
            pass
        else:
            self.first_estimate(0, 0, globalminima=True)

    def update_total_fit(self):
        total_fit = zeros(len(self.total_fit))

        for slice in self.Slices:
            total_fit += slice.gausscurve
        self.total_fit = total_fit

    def update_avg_sigma(self):
        sum = 0
        total_amplitude = 0
        for slice in self.Slices:
            sum += slice.sigma * slice.amplitude
            total_amplitude += slice.amplitude
        self.average_sigma = sum/total_amplitude
        average = self.average_sigma / self.fact
        self.FWHM = array(2.355 * average)

    def get_lowpass_signal_plot(self):

        fig = Figure()
        axes = fig.add_subplot(111)
        self.plot_lowpass_signal(axes)
        fig.set_tight_layout(True)

        return fig

    def plot_lowpass_signal(self, axes):

        axes.clear()
        axes.plot(self.photE, self.shiftraw, c='black')
        axes.plot(self.photE,self.filtered_intensity, c='orange')

        axes.set_xlabel('photon energy [eV]')
        axes.set_ylabel('intensity')
        axes.set_title('spectrum ' + str(self.spectrum_label))
        axes.grid(True)

    def get_peak_plot(self):
        fig = Figure()
        axes = fig.add_subplot(111)

        axes.plot(self.photE, self.shiftraw, c='black')
        axes.plot(self.photE, self.filtered_intensity, c='orange')
        allpeaks = find_peaks(self.filtered_intensity)[0]
        allypeaks = self.filtered_intensity[allpeaks]

        axes.scatter(self.photE[allpeaks], allypeaks, c='red', marker='x', zorder=5)
        axes.scatter(self.photE[self.xpeaks], self.ypeaks, c='green', marker='x', zorder=5, s=80)
        axes.set_xlabel('photon energy [eV]')
        axes.set_ylabel('intensity')
        axes.set_title('spectrum ' + str(self.spectrum_label))
        axes.grid(True)
        fig.set_tight_layout(True)

        return fig

    def plot_gaussian(self, axes):

        axes.clear()
        axes.plot(self.photE, self.shiftraw, color='black')
        axes.plot(self.photE, self.total_fit, color='red')
        axes.set_title('spectrum' + str(self.spectrum_label))
        axes.set_xlabel('photon energy [eV]')
        axes.set_ylabel('intensity')
        axes.grid(True)

    def get_gaussian_plot(self):

        fig = Figure()
        axes = fig.add_subplot(111)
        self.plot_gaussian(axes)
        fig.set_tight_layout(True)

        return fig
