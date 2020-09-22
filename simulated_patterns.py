import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.io.cif import CifParser
from model_functions import GaussianSum
from random import uniform, choice

class simulated_pattern(object):

    def __init__(self, cif, max_2theta=40.0):

        cif = CifParser(cif)
        self.struct = cif.get_structures(primitive=False)[0]
        self.max_2theta = max_2theta

    def change_lattice_constants(self, uc_params):

        a,b,c,alpha,beta,gamma = uc_params
        pi = np.pi
        ax = a
        ay = 0.0
        az = 0.0
        bx = b * np.cos(gamma * pi / 180.0)
        by = b * np.sin(gamma * pi / 180.0)
        bz = 0.0
        cx = c * np.cos(beta * pi / 180.0)
        cy = (c * b * np.cos(alpha * pi / 180.0) - bx * cx) / by
        cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
        unit_cell = np.array([[ax,ay,az],[bx,by,bz],[cx,cy,cz]])
        
        new_lattice = Lattice(unit_cell)
        self.struct.lattice = new_lattice

    def randomly_perturb_lattice_constants(self, pmin, pmax):

        uc_params0 = self.struct.lattice.abc + self.struct.lattice.angles
       
        perturbed_params_len = [p + choice((-1,1))*uniform(pmin,pmax)*p for p in uc_params0[0:3]]
        #perturbed_params_ang = [p + choice((-1,1))*uniform(pmin/2.0,pmax/2.0)*p for p in uc_params0[3:]]
        perturbed_params_ang = list(uc_params0[3:])
        perturbed_params = perturbed_params_len + perturbed_params_ang

        self.change_lattice_constants(perturbed_params)
            
    def find_peaks(self):

        XRD = XRDCalculator(symprec=0, wavelength='CuKa').get_pattern(self.struct, scaled=True, two_theta_range=(0,self.max_2theta))
        
        return np.array([XRD.x, XRD.y]).T

    def fit_GaussianModel(self, uc_params, mode='aniso', sd=0.05):

        if mode == 'iso':
            uc_params = np.array([uc_params[0],uc_params[0],uc_params[0],uc_params[1],uc_params[2],uc_params[3]])
        elif mode == 'ab':
            uc_params = np.array([uc_params[0],uc_params[0],uc_params[1],uc_params[2],uc_params[3],uc_params[4]])
        elif mode == 'ac':
            uc_params = np.array([uc_params[0],uc_params[1],uc_params[0],uc_params[2],uc_params[3],uc_params[4]])
        elif mode == 'bc':
            uc_params = np.array([uc_params[0],uc_params[1],uc_params[1],uc_params[2],uc_params[3],uc_params[4]])

        self.change_lattice_constants(uc_params)
        peaks = self.find_peaks()
        data = np.c_[peaks, np.full((len(peaks),), fill_value=sd)]
        model = lambda x: GaussianSum(x, data)

        return model

    def model_pattern_predictions(self, xvals, uc_params, mode='aniso'):

        model = self.fit_GaussianModel(uc_params, mode=mode)
        yvals = model(xvals)
        yvals /= max(yvals)

        return np.array([xvals,yvals]).T
