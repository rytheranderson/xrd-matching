from pymatgen.io.cif import CifParser
import os
import glob

def get_lattice_parameters(cif):

	struct = CifParser(cif).get_structures(primitive=False)[0]
	uc_params = struct.lattice.abc + struct.lattice.angles

	return uc_params

def print_lattice_parameter_matching_data(directory):

	cifs = glob.glob(directory + os.sep + '*.cif')
	MOFs = set([f.split('/')[-1].split('.')[0].split('_')[0] for f in cifs])
	results_dict = dict((m,{}) for m in MOFs)

	for mc in sorted(glob.glob(directory + os.sep + '*matched.cif')):

		MOF = mc.split('/')[-1].split('.')[0].split('_')[0]
		
		ef = None
		if 'Chebyshev' in mc:
			ef = 'Chebyshev'
		elif 'Clark' in mc:
			ef = 'Clark'
		elif 'Euclidean' in mc and 'squared' not in mc:
			ef = 'Euclidean'
		elif 'MAE' in mc:
			ef = 'MAE'
		elif 'MSE' in mc:
			ef = 'MSE'
		elif 'city_block' in mc:
			ef = 'city_block'
		elif 'divergence' in mc:
			ef = 'divergence'
		elif 'probabilistic_symmetric' in mc:
			ef = 'probabilistic_symmetric'
		elif 'squared' in mc:
			ef = 'squared_Euclidean'
		elif 'Kullback_Leibler' in mc:
			ef = 'Kullback_Leibler'
		if ef == None:
			raise ValueError('error function not recognized for', mc)

		results_dict[MOF][ef] = get_lattice_parameters(mc)

	for MOF,ef_dict in results_dict.items():

		correct_uc_params = get_lattice_parameters(directory + os.sep + MOF + '.cif')

		for ef,uc_params in ef_dict.items():

			for pn, pm, pc in zip(['a','b','c','alpha','beta','gamma'], uc_params, correct_uc_params):
			
				print(MOF, ef, pn, pm, pc)

print_lattice_parameter_matching_data('8-12P_matched_cifs_and_patterns')
#print_lattice_parameter_matching_data('4-6P_matched_cifs_and_patterns')

