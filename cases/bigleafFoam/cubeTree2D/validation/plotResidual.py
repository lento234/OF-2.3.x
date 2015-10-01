#!/usr/local/bin/python
# -*- coding: utf-8 -*-
"""
PLOTRESIDUALS
=============

Plot convergence of residuals

References:
[1]_ : Bourdin, Patrick, and John D. Wilson. 2008.
       “Windbreak Aerodynamics: Is Computational Fluid Dynamics Reliable?”
       Boundary-Layer Meteorology 126 (2): 181–208.
           Fig. 2

[2]_ : Santiago, J. L., F. Martín, a. Cuerva, N. Bezdenejnykh, and A. Sanz-Andrés. 2007.
       “Experimental and Numerical Study of Wind Flow behind Windbreaks.”
       Atmospheric Environment 41 (30): 6406–20.
           Fig. 1

"""

# ----------------------------------------------------------
# ----------------------------------------------------------
# -------------- L O A D I N G   M O D U L E S   -----------
# ----------------------------------------------------------
# ----------------------------------------------------------

# Standard/computational modules
import sys
import os
import glob
import numpy as np
import pyFlowStat.LineContainer as LC
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

# Plotting modules
import matplotlib as mpl
from matplotlib import pyplot as plt
from pyFlowStat import plotEnv

# ----------------------------------------------------------
# ----------------------------------------------------------
# -- S E T U P   P L O T T I N G   E N V I R O N M E N T  --
# ----------------------------------------------------------
# ----------------------------------------------------------

if 'print' in sys.argv[1:]:
    interactive = False
    saveFig     = True
else:
    interactive = True
    saveFig     = False

# Custom plotting environment
palette = plotEnv.set(plotType='line',
                      numColors=5,
                      interactive=interactive)

# Additional modifications
# mpl.rcParams.update({'legend.fontsize': 8})

# ----------------------------------------------------------
# ----------------------------------------------------------
# -------------- L O A D I N G   P A R A M E T E R S -------
# ----------------------------------------------------------
# ----------------------------------------------------------

if os.path.exists('../logs'):
    dataDir = glob.glob('../logs/*_0')
else:
    os.system('cd .. && foamLog log.solver')
    dataDir = glob.glob('../logs/*_0')


data = {}

for dataPath in dataDir:
	data[dataPath[8:-2]] = np.loadtxt(dataPath).T


# Plotting

plt.figure()
plt.xlabel(r'Iterations')
plt.ylabel(r'Rel. Tolerance')
plt.semilogy(data['p'][0], 	data['p'][1], 	label=r'p')
plt.semilogy(data['Ux'][0], data['Ux'][1], 	label=r'Ux')
plt.semilogy(data['Uz'][0], data['Uz'][1], 	label=r'Uz')
plt.semilogy(data['k'][0], 	data['k'][1], 	label=r'k')
plt.semilogy(data['T'][0], 	data['T'][1], 	label=r'T')
plt.semilogy(data['epsilon'][0], data['epsilon'][1], label=r'epsilon')

plt.legend(loc=0)
plotEnv.cleanupFigure()
if saveFig: plt.savefig('residual_initial')


plt.figure()
plt.title('Final Res')
plt.semilogy(data['pFinalRes'][0], 	data['pFinalRes'][1], 	label=r'p')
plt.semilogy(data['UxFinalRes'][0], data['UxFinalRes'][1], 	label=r'Ux')
plt.semilogy(data['UzFinalRes'][0], data['UzFinalRes'][1], 	label=r'Uz')
plt.semilogy(data['kFinalRes'][0], 	data['kFinalRes'][1], 	label=r'k')
plt.semilogy(data['epsilonFinalRes'][0], data['epsilonFinalRes'][1], label=r'epsilon')

plt.xlabel(r'Iterations')
plt.ylabel(r'Rel. Tolerance')
plt.legend(loc=0)
plotEnv.cleanupFigure()
if saveFig: plt.savefig('residual_final')
