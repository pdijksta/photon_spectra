from photon_spectra import analysis_procedure

def get_input_parameters():
    return {
        'height': 20,
        'prominence': 10,
        'spike-width': 1,
        'auto-correlation': 0,
        'parameter optimization': 1,
        'pulse duration correction': 0,
        'snr': 0.05,
        'sbr': 0.15,
        'intensity_thresh': 0.4,
        'roughness': 200,
        'frequency_cutoff': 0.0001,
        'background_prominence': 1.5,
        'cutoff_prominence': 5.0,
        'cutoff_height': 70.0,
        }.copy()

def analyze_spectrum(file_, parameters):
    parameters['h5_filename'] = file_
    analyzer = analysis_procedure.AnalysisProcedure()
    results_dict = analyzer.PerformMeas(parameters)
    return analyzer, results_dict

