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
import os
import glob
import numpy as np
import pyFlowStat.old.LineContainer as LC
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

# Plotting modules
import matplotlib as mpl
from matplotlib import pyplot as plt
#import plotenv

# ----------------------------------------------------------
# ----------------------------------------------------------
# -- S E T U P   P L O T T I N G   E N V I R O N M E N T  --
# ----------------------------------------------------------
# ----------------------------------------------------------

if 'print' in sys.argv[1:]:
    interactive = False
    saveFig = True
elif 'printSave' in sys.argv[1:]:
    interactive = False
    saveFig = True
    saveData = True
else:
    interactive = True
    saveFig = False
    saveData = False

# # Custom plotting environment
# palette = plotenv.set(plotType='line',
#                       numColors=9,
#                       interactive=interactive)

# Additional modifications
#mpl.rcParams.update({'legend.fontsize': 16})
mpl.rcParams.update({'axes.labelsize': 'large'})
plt.ion()

# ----------------------------------------------------------
# ----------------------------------------------------------
# -------------- L O A D I N G   P A R A M E T E R S -------
# ----------------------------------------------------------
# ----------------------------------------------------------

vegProp  = ParsedParameterFile('../constant/vegetationProperties')
H = vegProp['H'][-1]
C = vegProp['C'][-1]
l = vegProp['l'][-1]
kc = vegProp['kc'][-1]
Rg0 = vegProp['Rg0'][-1]
Rl0 = vegProp['Rl0'][-1]
cpa = vegProp['cpa'][-1]
rhoa = vegProp['rhoa'][-1]
rsMin = vegProp['rsMin'][-1]
a1 = vegProp['a1'][-1]
a2 = vegProp['a2'][-1]
a3 = vegProp['a3'][-1]
D0 = vegProp['D0'][-1]

transProp = ParsedParameterFile('../constant/transportProperties')
nu   = transProp['nu'][-1]
TRef = transProp['TRef'][-1]
qRef = transProp['qRef'][-1]

# Locate sampled data
dataDir = glob.glob('../postProcessing/sets/*/')[0]

# Load sampled data
lc = LC.LineContainer.createFromFoamFolder(dataDir, names=['xMidVegV'])

x,y,z = lc.lines['xMidVegV' + '_' + 'T'].xyz.T[:,:-2]
T = lc.lines['xMidVegV' + '_' + 'T'].rawVars()[:-2]-273.85
Tl = lc.lines['xMidVegV' + '_' + 'Tl'].rawVars()[:-2]-273.85
Rn = lc.lines['xMidVegV' + '_' + 'Rn'].rawVars()[:-2]
Ql = lc.lines['xMidVegV' + '_' + 'Ql'].rawVars()[:-2]
Qs = lc.lines['xMidVegV' + '_' + 'Qs'].rawVars()[:-2]


# # Plot leaf temperature
# plt.figure('T')
# plt.plot(T,z-0.5,'-',c=palette[0], label='$T_a$')
# plt.plot(Tl,z-0.5,'-',c=palette[1], label='$T_l$')
# plt.xlabel(r'$T ^{\circ}$C ')
# plt.ylabel(r'$\tilde{z}/H$ [m]')
# plt.legend()
# plotenv.cleanupFigure()
# if saveFig: plt.savefig('T.png',dpi=400)
#
# plt.figure('dT')
# plt.plot(Tl-T,z-0.5,'-',c=palette[0])
# plt.xlabel(r'$\Delta T ^{\circ}$C ')
# plt.ylabel(r'$\tilde{z}/H$ [m]')
# plt.legend()
# plt.plot(z*0.,z-0.5,'-',c='k')
# #plotenv.cleanupFigure()
# if saveFig: plt.savefig('dT.png',dpi=400)

plt.figure('Energy balance',figsize=(8,5))
plt.plot(-Rn,z-0.5,'-',c='0',label='$R_{n}$')
plt.plot(Ql,z-0.5,'--',c='0',label='$\mathrm{LAD}{\cdot}Q_{l}$')
plt.plot(Qs,z-0.5,':',c='0',label='$\mathrm{LAD}{\cdot}Q_{s}$')
plt.xlabel(r'$W/m{^3}$')
plt.ylabel(r'$\tilde{z}/H$ [m]')
plt.legend(loc=0)
plt.plot(z*0.,z-0.5,'-',c='k')
#plotenv.cleanupFigure()
if saveFig: plt.savefig('EnergyBalance.png',dpi=400)


fig, ax1 = plt.subplots(figsize=(8,5))

ax2 = ax1.twiny()
ax1.plot(T,z-0.5,'--',c='0', label='$T$')
ax1.plot(Tl,z-0.5,':',c='0', label='$T_l$')
ax1.set_xlabel(r'$T ^{\circ}$C ')
ax1.set_ylabel(r'$\tilde{z}/H$ [m]')

ax2.plot(Tl-T,z-0.5,'-',c='k',label='$\Delta T$')
ax2.set_xlabel(r'$\Delta T ^{\circ}$C ')
ax2.set_xlim(-5,5)

ax1.legend(loc=4)
ax2.legend(loc=2)
#plotenv.cleanupFigure()
plt.tight_layout()



# Store data
if saveData:
    curwd = os.getcwd().split('/')[-2]
    np.savez('../'+'data_'+curwd+'.npz',Tl=Tl,T=T,z=z,Rn=Rn,Ql=Ql,Qs=Qs)
