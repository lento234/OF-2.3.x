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
import pyFlowStat.old.LineContainer as LC
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

# Plotting modules
import matplotlib as mpl
from matplotlib import pyplot as plt
from pyFlowStat.old import TriSurfaceContainer as TC
import plotenv

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
palette = plotenv.set(plotType='surface',
                      interactive=interactive)

# ----------------------------------------------------------
# ----------------------------------------------------------
# ----------- L O A D I N G   P A R A M E T E R S ----------
# ----------------------------------------------------------
# ----------------------------------------------------------

# Load flow properties (initial conditions)
flowProp = ParsedParameterFile('../constant/flowProperties')
# URef = 2.0
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
levels = np.linspace(0,0.08,84) #np.hstack((, 0.99,1.01, np.arange(1.1,2.1,0.1)))
ticks  = np.linspace(0,0.08,11) #np.hstack((np.arange(0,1,0.2), 1.0, np.arange(1.2,2.1,0.2)))

plt.figure(r'Normalized Mean velocity field - U/U_HB')
plt.tricontourf(tsc.x, tsc.y, tsc.triangles,tsc.fields['U'].Umag(),
		        levels,cmap=palette['SPECTRAL']['DIV'])
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.85,0.85,1.3,1.3,0.85],'g',lw=0.5)
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.8,0.8,0.85,0.85,0.8],'r',lw=0.5)

plt.axis('scaled')
plt.axis([-3,8,0,3])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

plotenv.colorbar(ticks,orientation='h',format='%g')
plotenv.cleanupFigure(despine=False)
if saveFig: plt.savefig('surface_UField.png')

# ---------------------------------------------------------- Temperature
levels = np.linspace(298,303,21)-273.15 #np.hstack((, 0.99,1.01, np.arange(1.1,2.1,0.1)))
ticks  = np.linspace(298,303,11)-273.15 #np.hstack((np.arange(0,1,0.2), 1.0, np.arange(1.2,2.1,0.2)))

plt.figure(r'Temperature - T - 2')
plt.tricontourf(tsc.x, tsc.y, tsc.triangles,tsc.fields['T'].s-273.15,
		        levels,cmap=palette['SPECTRAL']['COLD_R'])
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.85,0.85,1.3,1.3,0.85],'g',lw=0.5)
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.8,0.8,0.85,0.85,0.8],'r',lw=0.5)

plt.axis('scaled')
plt.axis([-3,8,0,3])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

plotenv.colorbar(ticks,orientation='h',format='%.1f')
plotenv.cleanupFigure(despine=False)
if saveFig: plt.savefig('surface_TField.png')


# ---------------------------------------------------------- Temperature
levels = np.linspace(0.005,0.014,19)*1000 #np.hstack((, 0.99,1.01, np.arange(1.1,2.1,0.1)))
ticks  = np.linspace(0.005,0.014,10)*1000 #np.hstack((np.arange(0,1,0.2), 1.0, np.arange(1.2,2.1,0.2)))

plt.figure(r'Temperature - q - 2')
plt.tricontourf(tsc.x, tsc.y, tsc.triangles,tsc.fields['q'].s*1000,
		        levels,cmap='Blues')
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.85,0.85,1.3,1.3,0.85],'g',lw=0.5)
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.8,0.8,0.85,0.85,0.8],'r',lw=0.5)
plt.axis('scaled')
plt.axis([-3,8,0,3])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

# cb = plotEnv.colorbar(ticks,orientation='h')
plotenv.colorbar(ticks,orientation='h',format='%g')
# cb.formatter.set_powerlimits((0, 0))
plotenv.cleanupFigure(despine=False)
if saveFig: plt.savefig('surface_qField.png')

# e	Kref Tek	TE	WCI	WCT	H	WBGT	mR	Mrt	UTCI
bioklima = np.loadtxt('exportdata_bioklima.dat',skiprows=1,delimiter='\t')

levels = np.linspace(298,303,21)-273.15 #np.hstack((, 0.99,1.01, np.arange(1.1,2.1,0.1)))
ticks  = np.linspace(298,303,11)-273.15 #np.hstack((np.arange(0,1,0.2), 1.0, np.arange(1.2,2.1,0.2)))

plt.figure(r'bioklima - WBGT')
plt.tricontourf(tsc.x, tsc.y, bioklima[:,0],
		        levels,cmap=palette['SPECTRAL']['COLD_R'])
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.85,0.85,1.3,1.3,0.85],'g',lw=0.5)
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.8,0.8,0.85,0.85,0.8],'r',lw=0.5)
plt.axis('scaled')
plt.axis([-3,8,0,3])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plotenv.colorbar(ticks,orientation='h',format='%.1f')
plotenv.cleanupFigure(despine=False)



levels = np.linspace(20,40,21) #np.hstack((, 0.99,1.01, np.arange(1.1,2.1,0.1)))
ticks  = np.linspace(20,40,11) #np.hstack((np.arange(0,1,0.2), 1.0, np.arange(1.2,2.1,0.2)))

plt.figure(r'bioklima - UTCI')
plt.tricontourf(tsc.x, tsc.y, bioklima[:,1],
		        levels,cmap='jet',extend='both')
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.85,0.85,1.3,1.3,0.85],'g',lw=0.5)
plt.plot([-1.57,1.57,1.57,-1.57,-1.57],[0.8,0.8,0.85,0.85,0.8],'r',lw=0.5)
plt.axis('scaled')
plt.axis([-3,8,0,3])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plotenv.colorbar(ticks,orientation='h',format='%.1f')
plotenv.cleanupFigure(despine=False)
