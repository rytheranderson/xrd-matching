import numpy as np

def Gaussian(x, center, height, sd):

    return height*np.exp(-(x-center)**2/(2*sd**2))

def GaussianSum(x, data):

    return sum([Gaussian(x, row[0], row[1], row[2]) for row in data])
