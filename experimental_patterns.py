import numpy as np
from scipy import signal, interpolate
from model_functions import GaussianSum

class experimental_pattern(object):

    def __init__(self, filename, max_2theta=40.0, delimiter=''):

        exp_pattern = np.genfromtxt(filename, delimiter=delimiter)
        exp_pattern = exp_pattern[:,0:2]
        # normalize max peak height to unity
        exp_pattern[:,1] /= max(exp_pattern[:,1])
        exp_pattern[:,1] -= min(exp_pattern[:,1])
        exp_pattern = exp_pattern[exp_pattern[:,0] <= max_2theta]

        self.experimental_pattern = exp_pattern

    def find_peaks(self, width=5, prominence=0.005):

        """
            Width and prominence should be sufficient to find all peaks.
            The algorithm is highly dependent on these values, the defaults
            here will not work in all cases: always check the peak finding with a plot.
        """

        exp_pattern = self.experimental_pattern
        peak_inds = signal.find_peaks(exp_pattern[:,1], width=width, prominence=prominence)[0]

        return peak_inds

    def fit_GaussianModel(self, sd=0.05, find_peaks_kwargs={'width':5, 'prominence':0.005}):

        """
            Creates a model of Gaussian functions from peak locations. The Gaussian functions have the same height
            as the peaks, with a constant standard deviation.
        """
        
        exp_pattern = self.experimental_pattern
        peak_inds = self.find_peaks(**find_peaks_kwargs)
        peaks = exp_pattern[peak_inds]
        data = np.c_[peaks, np.full((len(peaks),), fill_value=sd)]
        model = lambda x: GaussianSum(x, data)
        
        return model

    def fit_InterpolationModel(self):

        """
            An linear interpolation function is a faster and often more reliable way to model 
            an experimental PXRD, especially if the  pattern has sharp peaks without much noise.
        """

        exp_pattern = self.experimental_pattern
        x = exp_pattern[:,0]
        y = exp_pattern[:,1]
        model = interpolate.interp1d(x, y, kind='linear', fill_value='extrapolate')

        return model

    def model_pattern_predictions(self, xvals, method='GaussianModel', find_peaks_kwargs={'width':3, 'prominence':0.005}):

        """
            Returns values for a pattern model for new values of 2theta (xvals). You can either use
            an interpolation model of fit a sum-of-Gaussians model.
        """

        xvals = np.array(xvals)
        if method == 'InterpolationModel':
            model = self.fit_InterpolationModel()
        elif method == 'GaussianModel':
            model = self.fit_GaussianModel(find_peaks_kwargs=find_peaks_kwargs)
        else:
            raise ValueError('Experimental model type should be GaussianModel or InterpolationModel')

        yvals = model(xvals)

        return np.array([xvals,yvals]).T
