from photon_spectra import analysis_procedure

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
    'roughness': 0.75,
    'frequency_cutoff': 0.05,
    'background_prominence': 6,
    'cutoff_prominence': 5,
    'cutoff_height': 20,
    }


def analyze_spectrum(file_, parameters):
    parameters['h5_filename'] = file_
    analyzer = analysis_procedure.AnalysisProcedure()
    results_dict = analyzer.PerformMeas(parameters)
    #import pdb; pdb.set_trace()
    return analyzer, results_dict

