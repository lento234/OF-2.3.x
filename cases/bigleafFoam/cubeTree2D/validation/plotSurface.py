#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
PLOTSURFACE
===========

Plot fields

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
from pyFlowStat import TriSurfaceContainer as TC

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
palette = plotEnv.set(plotType='surface',
                      interactive=interactive)

# ----------------------------------------------------------
# ----------------------------------------------------------
# ----------- L O A D I N G   P A R A M E T E R S ----------
# ----------------------------------------------------------
# ----------------------------------------------------------

# Load flow properties (initial conditions)
flowProp = ParsedParameterFile('../constant/flowProperties')
URef = 2.0
ustar = flowProp['ustar']
kappa = flowProp['kappa']
z0    = flowProp['z0']
# URef  = (ustar/kappa)*np.log(10.0/z0)

# Load surface data
dataDir = glob.glob('../postProcessing/surfaces/*/*/')[0]
tsc = TC.TriSurfaceContainer.createFromFoamFolder(dataDir, (1,0,0))


# ----------------------------------------------------------
# ----------------------------------------------------------
# ------------- P L O T T I N G  F I E L D S  --------------
# ----------------------------------------------------------
# ----------------------------------------------------------

# ---------------------------------------------------------- Mean velocity
levels = np.arange(0,1.25,0.02) #np.hstack((, 0.99,1.01, np.arange(1.1,2.1,0.1)))
ticks  = np.arange(0,1.4,0.2) #np.hstack((np.arange(0,1,0.2), 1.0, np.arange(1.2,2.1,0.2)))

plt.figure(r'Normalized Mean velocity field - U/U_HB')
plt.tricontourf(tsc.x, tsc.y, tsc.triangles,tsc.fields['U'].Umag()/URef,
		        levels,cmap=palette['DIV'])
plt.plot([12,18,18,12,12],[4,4,10,10,4],'k',lw=0.5)
plt.axis([0,40,0,20])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

plotEnv.colorbar(ticks,orientation='h')
plotEnv.cleanupFigure(despine=False)
# if saveFig: plt.savefig('surface_UField')


"""
# ---------------------------------------------------------- Turb. Disspation rate

plt.figure(r'Turbulent dissipation rate field (logscale)')

# levelsMax = [np.floor(np.log10(np.abs(tsc.fields['epsilon'].s).min())),
# 		     np.ceil(np.log10(np.abs(tsc.fields['epsilon'].s).max()))]
# levels = np.logspace(levelsMax[0],levelsMax[1],np.abs(levelsMax).sum()*4+1)
# ticks  = np.logspace(levelsMax[0],levelsMax[1],np.abs(levelsMax).sum()+1)

plt.tricontourf(tsc.x/HB,tsc.y/HB,tsc.triangles,tsc.fields['epsilon'].s,
				levels,locator=mpl.ticker.LogLocator(), cmap=palette['SEQ_COLD'])#SEQCMAPCOLD)

plt.plot([-HB/200,HB/200,HB/200,-HB/200,-HB/200],[0,0,1,1,0],'k',lw=0.5)

plt.axis([-5,10,0,2])
plt.xlabel(r'$x/h_b$')
plt.ylabel(r'$y/h_b$')
plotEnv.colorbar(ticks,orientation='h',format='exp')
plotEnv.cleanupFigure(despine=False)
if saveFig: plt.savefig('surface_epsilonField')


# ---------------------------------------------------------- TKE

plt.figure(r'Normalized TKE field')

levelsMin = np.floor(tsc.fields['k'].s.min()/(ustar**2))
levelsMax = np.ceil(tsc.fields['k'].s.max()/(ustar**2))
levels = np.arange(levelsMin,levelsMax+np.spacing(1e3))
ticks = np.arange(levelsMin,levelsMax+np.spacing(1e3),2)

plt.tricontourf(tsc.x/HB, tsc.y/HB, tsc.triangles,tsc.fields['k'].s/(ustar**2),
			    levels, cmap=palette['SEQ_HOT'])
plt.plot([-HB/200,HB/200,HB/200,-HB/200,-HB/200],[0,0,1,1,0],'k',lw=0.5)
plt.axis([-5,10,0,2])
plt.xlabel(r'$x/h_b$')
plt.ylabel(r'$y/h_b$')
plotEnv.colorbar(ticks,orientation='h')
plotEnv.cleanupFigure(despine=False)
if saveFig: plt.savefig('surface_kField')

# ---------------------------------------------------------- Plot epsilon field

plt.figure(r'Turbulent viscosity field (logscale)')

levelsMax = [np.floor(np.log10(np.abs(tsc.fields['nut'].s).min())),
		     np.ceil(np.log10(np.abs(tsc.fields['nut'].s).max()))]
levels = np.logspace(levelsMax[0],levelsMax[1],np.abs(levelsMax).sum()*4+1)
ticks  = np.logspace(levelsMax[0],levelsMax[1],np.abs(levelsMax).sum()+1)

plt.tricontourf(tsc.x/HB, tsc.y/HB, tsc.triangles,tsc.fields['nut'].s,
				levels,locator=mpl.ticker.LogLocator(),cmap=palette['SEQ_COLD'])

plt.plot([-HB/200,HB/200,HB/200,-HB/200,-HB/200],[0,0,1,1,0],'k',lw=0.5)
plt.axis([-5,10,0,2])
plt.xlabel(r'$x/h_b$')
plt.ylabel(r'$y/h_b$')
plotEnv.colorbar(ticks,orientation='h',format='exp')
plotEnv.cleanupFigure(despine=False)
if saveFig: plt.savefig('surface_nutField')

# ---------------------------------------------------------- Mean velocity

levels    = np.hstack((np.mgrid[-30:0:3], -0.3, 0.3, np.mgrid[3:30+np.spacing(1e3):3]))
ticks     = np.hstack((np.mgrid[-30:0:6], 0, np.mgrid[6:30+np.spacing(1e3):6]))

plt.figure(r'Mean vorticity field')
plt.tricontourf(tsc.x/HB, tsc.y/HB, tsc.triangles,tsc.fields['vorticity'].vz,#*HB/U0HB,
                levels,cmap=palette['DIV'],extend='both')

plt.axis([-5,10,0,2])
plt.plot([-HB/200,HB/200,HB/200,-HB/200,-HB/200],[0,0,1,1,0],'k',lw=0.5)
plt.xlabel(r'$x/h_b$')
plt.ylabel(r'$y/h_b$')
plotEnv.colorbar(ticks,orientation='h')
plotEnv.cleanupFigure(despine=False)
if saveFig: plt.savefig('surface_meanVorticity')

"""

# ---------------------------------------------------------- Pressure
# plt.figure(r'Normalized Pressure field - P/(U_h*h)**2')
#
# levelsMax = np.ceil(np.abs([tsc.fields['p'].s.min(),tsc.fields['p'].s.max()]).min()/((U0HB**2)*(HB**2))*100)/100
# levels    = np.hstack((np.mgrid[-0.22:0:0.02], -0.002, 0.002, np.mgrid[0.02:0.22+np.spacing(1e3):0.02]))
# ticks     = np.hstack((np.mgrid[-0.22:0:0.06], 0, np.mgrid[0.04:0.22+np.spacing(1e3):0.06]))
#
# plt.tricontourf(tsc.x/HB, tsc.y/HB, tsc.triangles,tsc.fields['p'].s/((U0HB**2)*(HB**2)),levels,
#                 cmap=palette['DIV'], extend='both')
#
# plt.plot([-HB/200,HB/200,HB/200,-HB/200,-HB/200],[0,0,1,1,0],'k',lw=0.5)
# plt.axis([-5,10,0,2])
# plt.xlabel(r'$x/h_b$')
# plt.ylabel(r'$y/h_b$')
#
# plotEnv.colorbar(ticks,orientation='h')
# plotEnv.cleanupFigure(despine=False)
# if saveFig: plt.savefig('surface_pField')
