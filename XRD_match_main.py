from pattern_matching import SimExp_pattern_matching
from experimental_patterns import experimental_pattern
from simulated_patterns import simulated_pattern
import pattern_matching as pm
import argparse
import numpy as np
import time
import configuration as config

GA = config.general_arguments
OA = config.optimizer_arguments

def run_matching():

	parser = argparse.ArgumentParser(description='Optional arguments for XRD matching')
	parser.add_argument('-E', '--exp_pattern', action='store', dest='exp_pattern', type=str, required=True,
		help='the experimental pattern to match,')
	parser.add_argument('-S', '--sim_cif', action='store', dest='sim_cif', type=str, required=True,
		help='the cif to match to the experimental pattern')
	parser.add_argument('-O', '--optimizer', action='store', dest='optimizer', type=str, required=True,
		help='the global optimizer to use')
	parser.add_argument('-F', '--error_function', action='store', dest='error_function', type=str, required=False, default=None,
		help='the error metric to minimize')

	ef_dict = {
	'city_block': pm.city_block,
	'MAE': pm.MAE,
	'Euclidean': pm.Euclidean,
	'squared_Euclidean': pm.squared_Euclidean,
	'MSE': pm.MSE,
	'Chebyshev': pm.Chebyshev,
	'Clark': pm.Clark,
	'probabilistic_symmetric': pm.probabilistic_symmetric,
	'divergence': pm.divergence,
	'Kullback_Leibler': pm.Kullback_Leibler
	}

	args = parser.parse_args()

	exp = experimental_pattern(args.exp_pattern, max_2theta=GA['max_2theta'], delimiter=GA['exp_pattern_delimiter'])
	sim = simulated_pattern(args.sim_cif, max_2theta=GA['max_2theta'] + 0.5)

	if args.error_function == None:
		error_function = GA['error_function']
		print('error function specified in configuration.py', str(error_function).split()[1])
	else:
		error_function = ef_dict[args.error_function]
		print('error function specified in terminal as', str(error_function).split()[1])

	if GA['randomly_perturb'][0]:
		print('randomly perturbing starting lattice constants...')
		sim.randomly_perturb_lattice_constants(GA['randomly_perturb'][1], GA['randomly_perturb'][2])

	xvals = np.linspace(GA['min_2theta'], GA['max_2theta'], GA['num_points'])
			   
	fname = args.sim_cif.split('/')[-1].split('.')[0]

	matching = SimExp_pattern_matching(exp, sim, xvals, 
									   mode=GA['constraint_mode'], 
									   fixed=GA['fixed_parameters'], 
									   error_function=error_function, 
	                                   find_peaks_kwargs={'width':GA['peak_find_width'], 'prominence':0.005},
	                                   outname=fname, wavelength=GA['wavelength'],
	                                   exp_model_method=GA['experimental_model'])

	matching.run(GA['outer_iterations'], GA['max_fraction_change'], args.optimizer, OA[args.optimizer], GA['halting_criteria'])

if __name__ == '__main__':

	start_time = time.time()
	run_matching()
	print('Normal termination of XRD_match after')
	print('--- %s seconds ---' % (time.time() - start_time))
	print()
