#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
PLOTHORIZPROFILE
================

Plot horizontal variation of quantities

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
    saveFig = True
else:
    interactive = True
    saveFig = False

# Custom plotting environment
palette = plotEnv.set(plotType='line',
                      numColors=1,
                      interactive=interactive)

# Additional modifications
mpl.rcParams.update({'legend.fontsize': 8})

# ----------------------------------------------------------
# ----------------------------------------------------------
# -------------- L O A D I N G   P A R A M E T E R S -------
# ----------------------------------------------------------
# ----------------------------------------------------------

# Load flow properties (initial conditions)
flowProp = ParsedParameterFile('../constant/flowProperties')
ustar = flowProp['ustar']
kappa = flowProp['kappa']
z0    = flowProp['z0']

Uref = 2.0

# Locate sampled data
dataDir = glob.glob('../postProcessing/sets/*/')[0]
varNames= ['U']
lnNames = ['xMidVeg']
labels  = ['midVeg']

# Load sampled data
lc = LC.LineContainer.createFromFoamFolder(dataDir, names=lnNames, underscoreHeaders=varNames)

# ----------------------------------------------------------
# ----------------------------------------------------------
# --------------P L O T T I N G ----------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------

# --------------------------
# P L O T   V E L O C I T Y
# --------------------------

# line = lc.lines[lnNames[0] + '_' + 'U']
# z    = line.xyz[:,2]
# U    = line.vx
# U0  = (ustar/kappa)*np.log(z/z0)
#
# plt.figure('U')
# plt.plot(U/Uref, z, '-', c=palette[0])
# plt.plot(U0/Uref, z ,'--', c=palette[0])
# curax = plt.axis()
# plt.xlabel(r'$\overline{u}/U_{\mathit{ref}}$')
# plt.ylabel(r'$z$ [m]')
# # plt.legend(loc=0, ncol=2)
# plotEnv.cleanupFigure()


# ----------------
# P L O T   L A I
# ----------------

z = lc.lines[lnNames[0] + '_' + 'LAI'].xyz[:,2]
LAI = lc.lines[lnNames[0] + '_' + 'LAI'].s

plt.figure('LAI')
plt.plot(LAI, z, '-', c=palette[0])
curax = plt.axis()
plt.xlabel(r'$LAI$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()

# ---------------------------
# P L O T   R A D I A T I O N
# ---------------------------

Rg = lc.lines[lnNames[0] + '_' + 'Rg'].vz
Rg2 = 100*np.exp(-0.75*LAI)
Rg3 = (100*np.exp(-0.75*1*(10.0-z)))[:71]
plt.figure('Rg')
plt.plot(Rg, z, '-', c=palette[0])
plt.plot(100*np.exp(-0.75*LAI), z, '.--', c='r')
plt.plot((100*np.exp(-0.75*1*(10.0-z)))[:71], z[:71], '.--', c='g')
curax = plt.axis()
plt.xlabel(r'$Rg(z)$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()

Rn = lc.lines[lnNames[0] + '_' + 'Rn'].s
dz = z[1:]-z[:-1]
zMid = 0.5*(z[1:]+z[:-1])

plt.figure('Rn')
plt.plot(Rn, z, '-', c=palette[0])
plt.semilogx((Rg[1:]-Rg[:-1])/dz, zMid, '--', c='r')
plt.semilogx((Rg2[1:]-Rg2[:-1])/dz, zMid, ':', c='g')
plt.semilogx((Rg3[1:]-Rg3[:-1])/dz[:70], zMid[:70], '--', c='b')
curax = plt.axis()
plt.xlabel(r'$Rn$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()

# ----------------
# temperature
# ----------------

T = lc.lines[lnNames[0]+'_'+'T'].s
plt.figure('T')
plt.plot(T, z, c=palette[0])
plt.xlabel(r'$T$')
plt.ylabel(r'$z$')
plotEnv.cleanupFigure()


rhosat = lc.lines[lnNames[0]+'_'+'rhosat'].s
plt.figure('rhosat')
plt.plot(rhosat, z, c=palette[0])
plt.plot(np.exp(), z, c='r')
plt.xlabel(r'$\rho*$')
plt.ylabel(r'$z$')
plotEnv.cleanupFigure()

"""
plt.semilogy(data['k'][0], 	data['k'][1], 	label=r'k')
plt.figure('$u$ vs. Z/Zref')
for i,label in enumerate(labels):

    # Load data
    line = lc.lines[lnNames[i] + '_' + 'U']
    xoHB = line.xyz[:,0]/HB
    UoU04 = np.linalg.norm(line.rawVars(), axis=1)/U04

    # Plot numerical results
    plt.plot(xoHB, UoU04, '-', c=palette[i], label=label + ' Present Study')

    # Plot validation data
    plt.plot(valData['BnW_rKE'][lnNames[i]][:,0],
             valData['BnW_rKE'][lnNames[i]][:,1], '.--', c=palette[i],
             label=label + ' B&W 2007, Num')

    plt.plot(valData['santiagoetal_rKE'][lnNames[i]][:,0],
             valData['santiagoetal_rKE'][lnNames[i]][:,1], '*--', c=palette[i],
             label=label + ' S. et al. 2007, Num')

    plt.plot(valData['BnM'][lnNames[i]][:,0],
             valData['BnM'][lnNames[i]][:,1], 'o', c='k',
             label=label + ' B&M 1983, Exp')

# Modify Figure
plt.axis([-5, 35, 0.1, 1])
curax = plt.axis()
plt.xlabel(r'$x/{h_b}$')
plt.ylabel(r'$u/u_{04}$')
plt.legend(loc=0, ncol=2)
plotEnv.cleanupFigure()
if saveFig: plt.savefig('plots/horizProfile')
"""
