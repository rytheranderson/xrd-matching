# cif2lammps
## Authors

- Ryther Anderson

## Motivation
XRD_matching is a Python 3 program which matches the XRD pattern of a simulated material with an experimental XRD pattern by varying the lattice parameters of the simulated material.
This is what Rietveld refinement seeks to do. However, Rietveld refinement uses least squares, which does not work well for complex materials such as MOFs, or when significant peak splitting is occuring.
The implementation here uses differential evolution (from scipy) to find a global minimum of an arbitrary similarity metric (i.e. error function).

## Current Status
The code has been tested on 10 materials using 10 error functions by perturbing the unit cell parameters of each material by random amounts and seeing if they can be returned to their
original values. Perturbations of up to 12% of the original parameter value were used (this is a very large perturbation in crystallography). The results show that this code is capable of 
quickly refining structures even with these large perturbations (that often result in significant peak splitting). An experimental XRD pattern can be represented as a linear interpolation function,
which is appropriate when the pattern has distinct peaks without much background. Alternatively, each peak in the experimental pattern can be represented as a Gaussian function (the peaks are 
automatically identified) of corresponding height and center (currently a constant standard deviation is used). This procedure removes any background and peak broadening.

## Usage
Generally speaking, just run:
```
python XRD_match_main.py -E experimental_pattern_file -S simulated_cif -O differential_evolution
```
where "experimental_pattern_file" is the experimental pattern (the first two columns should be 2theta and intensity, respectively), "simulated_cif" is the CIF of the structure to be refined 
(it can be in a variety of formats readable by pymatgen), and -O is the optimization method to use. Differential evolution seems to be the most robust optimization method (I will do more testing), 
but basin hopping and dual annealing can also be used. The configuration.py can also be used to alter other variables that may influence fitting, see the comments in this file for a brief description of 
each of these variables. The values currently in configuration.py are good defaults, but may not be appropriate for all refinements.

## Requirements
The Anaconda distribution of Python 3.X with pymatgen will work, otherwise have numpy, scipy, matplotlib, and pymatgen.
Install instructions for pymatgen can be found here: https://pymatgen.org/installation.html.

