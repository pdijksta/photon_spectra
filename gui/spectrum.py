from photon_spectra import analysis_procedure

default_input_parameters_xfel = {
    'height': 20,
    'prominence': 5,
    'spike-width': 1,
    'auto-correlation': 0,
    'parameter optimization': 0,
    'pulse duration correction': 0,
    'snr': 1,
    'sbr': 0.1,
    'intensity_thresh': 0.25,
    'frequency_cutoff': 0.7,
    'background_prominence': 6,
    'cutoff_prominence': 5,
    'cutoff_height': 20,
    }

default_input_parameters_swissfel = {
        'height': 20,
        'prominence': 10,
        'spike-width': 1,
        'auto-correlation': 0,
        'parameter optimization': 0,
        'pulse duration correction': 0,
        'snr': 0.05,
        'sbr': 0.15,
        'intensity_thresh': 0.1,
        'frequency_cutoff': 0.1,
        'background_prominence': 1.5,
        'cutoff_prominence': 10.0,
        'cutoff_height': 20.0,
        }

def analyze_spectrum(file_, parameters):
    parameters['h5_filename'] = file_
    analyzer = analysis_procedure.AnalysisProcedure()
    results_dict = analyzer.PerformMeas(parameters)
    return analyzer, results_dict

