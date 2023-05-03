import time
from .single_spectra import single_spectra
from numpy import array

def analyse_single_spectra(i, params, intense, photE, fact):
    spectra = single_spectra(intense, photE, fact, i)

    spectra.adaptive_butterworth_filter(params['frequency_cutoff'], params['roughness'], params['snr'], params['sbr'])

    if not (spectra.noisybool):

        spectra.finding_peak(params['prominence'], params['height'], params['background_prominence'])
        if spectra.n > 0:
            spectra.first_estimate(params['cutoff_prominence'], params['cutoff_height'], globalminima=False)
    
    return spectra

def feedback(i, params, intense, photE, fact):
    spectra = single_spectra(intense, photE, fact, i)
    alphas = array([5, 20, 50, 100, 200, 300])
    spectra.findbestalpha(alphas, params)
    return spectra
