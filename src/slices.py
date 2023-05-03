#import time
from numpy import arange, zeros_like, pi, exp
from numpy import min as npmin
from numpy import max as npmax
from numpy import sqrt as npsqrt
from lmfit.models import GaussianModel
from lmfit import Parameters, Model

def gaussian_central_offset(x, amplitude, sigma, offset):
    return amplitude/(npsqrt(2*pi)*sigma) * exp(-x**2 / (2*sigma**2)) + offset
    
def gaussian_offset(x, center, amplitude, sigma, offset):
    return amplitude/(npsqrt(2*pi)*sigma) * exp(-(x-center)**2 / (2*sigma**2)) + offset
    
class Slice:
    def __init__(self, intensity, leftindex, rightindex, model='gaussian'):
        x = arange(len(intensity))
        self.leftindex = leftindex 
        self.rightindex = rightindex
        self.intensity = zeros_like(intensity)
        self.intensity[leftindex:rightindex] = intensity[leftindex:rightindex]#-npmin(intensity[leftindex:rightindex])

        #print(self.intensity, intensity)
        
        if model == 'gaussian':
            mod = GaussianModel()
            pars = mod.guess(self.intensity, x=x)

            #pars = Parameters()
            #sigmaest = (rightindex-leftindex)/4
            #pars.add('sigma',value = sigmaest)
            #pars.add('center',value = (rightindex + leftindex)/2)
            #pars.add('amplitude',value =np.max(self.intensity)*(2*np.pi)**0.5*(rightindex-leftindex)/2)
            
            out = mod.fit(self.intensity, pars, x=x)
            self.sigma = out.params['sigma'].value
            self.amplitude = out.params['amplitude'].value
            self.center = out.params['center'].value

            self.gausscurve = out.best_fit
            #self.gausscurve[leftindex:rightindex] += npmin(intensity[leftindex:rightindex])
            
        elif model == 'central offset gaussian':
            goffsetmodel = Model(gaussian_offset)
            sigmaest = (rightindex-leftindex)/4
            ampest = npmax(self.intensity)*(2*pi)**0.5*(rightindex-leftindex)/4
            offsetest = npmin(self.intensity)
            result = goffsetmodel.fit(data=self.intensity, x=x, amp=ampest, sigma=sigmaest, offset=offsetest)
            self.sigma = result.params['sigma'].value
            self.amplitude = result.params['amplitude'].value
            self.center = 0

            self.gausscurve = result.best_fit
        elif model == 'offset gaussian':
            goffsetmodel = Model(gaussian_offset)
            sigmaest = (rightindex-leftindex)/4
            ampest = npmax(self.intensity)*(2*pi)**0.5*sigmaest
            offsetest = npmin(self.intensity[leftindex:rightindex])
            result = goffsetmodel.fit(data=self.intensity, x=x, centerest=(rightindex + leftindex)/2, amp=ampest, sigma=sigmaest, offset=offsetest)
            self.sigma = result.params['sigma'].value
            self.amplitude = result.params['amplitude'].value
            self.center = result.params['center'].value

            self.gausscurve = result.best_fit-result.params['offset'].value
            
        #self.redchi = out.redchi/self.max**2*10**3
