import os
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import differential_evolution, basinhopping, dual_annealing, minimize
from multiprocessing import cpu_count
from pymatgen.io.cif import CifWriter

Nfeval = 1

def city_block(x, y):

    return np.sum(np.abs(x-y))

def MAE(x, y):

    return np.average(np.abs(x-y))

def Euclidean(x, y):

    return np.sqrt(np.sum((x-y)*(x-y)))

def squared_Euclidean(x, y):

    return np.sum((x-y)*(x-y))

def MSE(x, y):

    return np.average((x-y)*(x-y))

def Chebyshev(x ,y):

    return np.max(np.abs(x-y))

def Clark(x, y):

    num = np.abs(x-y)
    den = x+y
    term = num/den
    term = term[np.logical_not(np.isnan(term))]
    term = term*term

    return np.sqrt(np.sum(term))

def probabilistic_symmetric(x, y):

    num = (x-y)*(x-y)
    den = x+y
    term = num/den
    term = term[np.logical_not(np.isnan(term))]

    return 2*np.sum(term)

def divergence(x, y):

    num = (x-y)*(x-y)
    den = (x+y)*(x+y)
    term = num/den
    term = term[np.logical_not(np.isnan(term))]

    return 2*np.sum(term)

def Kullback_Leibler(x, y):

    term = x/y
    term = x*np.log(term)
    term = term[np.logical_not(np.isnan(term))]
    term = term[np.logical_not(np.isinf(term))]

    return np.sum(term)

class SimExp_pattern_matching(object):

    def __init__(self, exp, sim, xvals, mode='aniso', fixed=[0,0,0,1,1,1], 
                 exp_model_method='GaussianModel', error_function=city_block,
                 find_peaks_kwargs={'width':3, 'prominence':0.005}, outname='output'):

        if not os.path.exists(outname):
            os.makedirs(outname)
            
        self.outname = outname
        self.mode = mode
        self.sim = sim
        self.exp = exp
        self.xvals = xvals
        self.find_peaks_kwargs = find_peaks_kwargs

        EMP = exp.model_pattern_predictions(xvals, method=exp_model_method, find_peaks_kwargs=find_peaks_kwargs)
        
        self.EMP = EMP
        self.error_function = error_function

        fixed_dict = {'aniso':fixed,
                      'iso':[fixed[0],fixed[3],fixed[4],fixed[5]],
                      'ab':[fixed[0],fixed[2],fixed[3],fixed[4],fixed[5]],
                      'ac':[fixed[0],fixed[1],fixed[3],fixed[4],fixed[5]],
                      'bc':[fixed[0],fixed[1],fixed[3],fixed[4],fixed[5]]}

        self.fixed = list(map(bool,fixed_dict[mode]))

    def objective(self, X):

        SMP = self.sim.model_pattern_predictions(self.xvals, X, mode=self.mode)
        sim_y = SMP[:,1]
        exp_y = self.EMP[:,1]
        error = self.error_function(sim_y, exp_y)

        return error

    def match_patterns(self, max_fraction_change, method, method_kwargs, halting_criteria):

        UC0 = self.sim.struct.lattice.abc + self.sim.struct.lattice.angles

        mode_vary_params = {'aniso':UC0,
                            'iso':[UC0[0],UC0[3],UC0[4],UC0[5]],
                            'ab' :[UC0[0],UC0[2],UC0[3],UC0[4],UC0[5]],
                            'ac' :[UC0[0],UC0[1],UC0[3],UC0[4],UC0[5]],
                            'bc' :[UC0[0],UC0[1],UC0[3],UC0[4],UC0[5]]}

        MFC = max_fraction_change
        vary_params = mode_vary_params[self.mode]
        bounds = [(p,p) if f else (p-MFC*p, p+MFC*p) for f,p in zip(self.fixed, vary_params)]

        print('Starting unit cell parameters:', [round(p, 4) for p in UC0])
        print('Bounds:', [(round(b[0],4), round(b[1],4)) for b in bounds])
        print()

        if method == 'differential_evolution':

            if halting_criteria != None:
                
                def callbackF(xk, convergence=100.0):
    
                    SMP = self.sim.model_pattern_predictions(self.xvals, xk)[:,1]
                    EMP = self.EMP[:,1]
                    crit = halting_criteria[0](SMP, EMP)
                    print('halting criteria value =', np.round(crit, 5))
    
                    if crit < halting_criteria[1]:
                        print('halting inner optimization due to halting criteria, val less than', halting_criteria[1])
                        return True

            else:
                callbackF = None

            res = differential_evolution(self.objective, bounds, **method_kwargs, 
                                         updating='deferred', polish=True, workers=cpu_count(), disp=True, callback=callbackF)
        elif method == 'basinhopping':
            res = basinhopping(self.objective, self.uc_params0, **method_kwargs, 
                               disp=True)
        elif method == 'dual_annealing':

            def callbackF(X, f, c):
                global Nfeval
                print(Nfeval, self.objective(X))
                Nfeval += 1

            res = dual_annealing(self.objective, bounds, **method_kwargs, 
                                 no_local_search=True, callback=callbackF)
            res = minimize(self.objective, res.x, method='L-BFGS-B', bounds=bounds)
        else:
            raise ValueError('optimizer', method, 'is not implemented.')

        UC = res.x
        UC = list(UC) + list(np.zeros(2))

        mode_uc_params = {'aniso':UC[0:6],
                          'iso':[UC[0],UC[0],UC[0],UC[1],UC[2],UC[3]],
                          'ab': [UC[0],UC[0],UC[1],UC[2],UC[3],UC[4]],
                          'ac': [UC[0],UC[1],UC[0],UC[2],UC[3],UC[4]],
                          'bc': [UC[0],UC[1],UC[1],UC[2],UC[3],UC[4]]}

        matched_uc_params = mode_uc_params[self.mode]
        self.sim.change_lattice_constants(matched_uc_params)
        self.matched_uc_params = matched_uc_params

        print('Final unit cell parameters:', [round(p, 4) for p in self.matched_uc_params])
        final_error = res.fun
        print('Final error:', round(final_error, 4))
        print()

        return matched_uc_params, final_error, res

    def plot_matching_resuts(self, SMP0):

        SMP1 = self.sim.model_pattern_predictions(self.xvals, self.matched_uc_params)

        fig, axes = plt.subplots(2, 1, figsize=(7.5, 5.0))
        axes[0].plot(self.xvals, SMP0[:,1] + 0.1, color='orange', label='Initial Simulated Pattern', linewidth=0.6)
        axes[0].plot(self.xvals, self.EMP[:,1], color='blue', label='Experimental Pattern Model', linewidth=0.5)
        axes[0].legend(loc='upper right')
        axes[0].set_ylabel('Normalized Intensity')
        axes[1].plot(self.xvals, SMP1[:,1] + 0.1, color='orange', label='Matched Simulated Pattern', linewidth=0.6)
        axes[1].plot(self.xvals, self.EMP[:,1], color='blue', label='Experimental Pattern Model', linewidth=0.5)
        axes[1].legend(loc='upper right')
        axes[1].set_ylabel('Normalized Intensity')
        axes[1].set_xlabel(r'$2\Theta$')

        fig.savefig(self.outname + os.sep + self.outname + '_matched_XRD.png', bbox_inches='tight', dpi=300)

    def write_matched_pattern(self, xvals, uc_params):

        pattern = self.sim.model_pattern_predictions(self.xvals, self.matched_uc_params)
        np.savetxt(self.outname + os.sep + self.outname + '_matched_pattern.txt', pattern, delimiter=' ')

    def write_cif(self):

        writer = CifWriter(self.sim.struct, symprec=None)
        writer.write_file(self.outname + os.sep + self.outname + '_matched.cif')

    def run(self, outer_iterations, max_fraction_change, method, method_kwargs, halting_criteria):

        UC0 = self.sim.struct.lattice.abc + self.sim.struct.lattice.angles
        SMP0 = self.sim.model_pattern_predictions(self.xvals, UC0)
        EMP = self.EMP[:,1]

        for i in range(outer_iterations):

            print('outer iteration:', i+1, 'out of', outer_iterations)
            print('max fraction change:', max_fraction_change/(i+1.0))
            uc_params, error, res = self.match_patterns(max_fraction_change/(i+1.0), method, method_kwargs, halting_criteria)

            if res.message == 'callback function requested stop early by returning True':
                print('halting outer iterations due to halting criteria')
                break

        self.plot_matching_resuts(SMP0)
        self.write_matched_pattern(self.xvals, self.matched_uc_params)
        self.write_cif()

        FSMP = self.sim.model_pattern_predictions(self.xvals, self.matched_uc_params)[:,1]
        mae = MAE(EMP, FSMP)
        mse = MSE(EMP, FSMP)
        clark = Clark(EMP, FSMP)
        ps = probabilistic_symmetric(EMP, FSMP)

        error_line = [mae, mse, clark, ps]
        efname = str(self.error_function).split()[1]

        mode_lines = {'aniso':'The cell was allowed to change anisotropically', 
                      'iso':'Patterns were matched with the constraint a = b = c',
                      'ab':'Patterns were matched with the constraint a = b',
                      'ac':'Patterns were matched with the constraint a = c',
                      'bc':'Patterns were matched with the constraint b = c'}

        with open(self.outname + os.sep + self.outname + '_summary.txt', 'w') as out:
            out.write(str(res))
            out.write('\n')
            out.write(mode_lines[self.mode] + '\n')
            out.write('Error function: ' + efname + '\n')
            out.write('UC param max fraction change: ' + str(max_fraction_change) + '\n')
            out.write('Final error: ' + str(error) + '\n')
            out.write('Final uc params: ' + ' '.join(list(map(str, self.matched_uc_params))) + '\n')
            out.write('{:<8} {:<8} {:<8} {:<8}'.format('MAE', 'MSE', 'Clark', 'PS'))
            out.write('\n')
            out.write('{:<8.6f} {:<8.6f} {:<8.2f} {:<8.4f}'.format(*error_line))
