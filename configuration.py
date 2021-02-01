import pattern_matching as pm

general_arguments = {
	'exp_pattern_delimiter':'',
	'outer_iterations': 1,
	'num_points': 1000,
	'min_2theta': 0.5,
	'max_2theta': 3.5,
	'constraint_mode': 'aniso',
	'fixed_parameters': (0,0,0,1,1,0),
	'error_function': pm.city_block,
	'peak_find_width': 3,
	'max_fraction_change': 0.06,
	'experimental_model': 'InterpolationModel',
	'randomly_perturb': (False, 0.08, 0.12),
	#'halting_criteria': (pm.MAE, 0.0012),
	'halting_criteria': None,
	'wavelength': 0.45256
}

optimizer_arguments = {
	'differential_evolution': {
	# a larger recombination, popsize and/or smaller mutation value may speed convergence but will 
	# result in a more narrow search. For more difficult cases parameters geared towards a broad 
	# search may work better.
		'maxiter': 30,
		'popsize': 15, 
		'recombination': 0.80,
		'mutation': (0.5,1.0),
		'strategy': 'randtobest1bin'
	},
	'basinhopping': {
		'niter': 100,
		'T': 1.0,
		'stepsize':0.5,
		'interval': 50,
	},
	'dual_annealing': {
		'maxiter': 500,
		'initial_temp': 5000.0,
		'restart_temp_ratio': 2e-05,
		'visit': 2.0,
		'accept': -0.1,
		'local_search_options': {'method':'CG'}
	}
}