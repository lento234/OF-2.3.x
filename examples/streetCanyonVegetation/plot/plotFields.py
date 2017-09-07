#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Plots fields
"""

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
import scipy.interpolate as spinterp

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


# ------------- P L O T T I N G  F I E L D S  --------------

#levels = np.arange(12)*0.25/2.0
#ticks = [0.,0.25,0.5,0.75,1.0,1.25,1.5]
levels = np.linspace(29,30.2,25)
ticks = np.linspace(29,30.2,7)

theta = np.linspace(0,2*np.pi,100)
xc = 0.5*np.cos(theta)
yc = 0.5*np.sin(theta)+1

mask = np.where((tsc.x>=-2.5) & (tsc.x<=6.5) & (tsc.y<=3.5))[0]

plt.figure(r'Normalized Mean velocity field - U/U_HB')
#plt.tricontourf(tsc.x, tsc.y, tsc.fields['U'].Umag(),levels, cmap='viridis')
plt.tricontourf(tsc.x, tsc.y, tsc.fields['T'].s-273.15,levels, cmap='viridis')

plt.plot(xc,yc,'--w',lw=2.5)

plt.axis('scaled')
plt.axis([-2,6,0,3])
plt.xlabel(r'$x/H$')
plt.ylabel(r'$z/H$')

cb = plotenv.colorbar(ticks,orientation='h',format='%g')
plotenv.add_text('(a)', color='w')
cb.set_label('$T$ ($^{\circ}$C)')

x,y = np.meshgrid(np.linspace(-2.5,6.5,180),np.linspace(0,3.5,70))
u = spinterp.griddata(np.vstack((tsc.x, tsc.y)).T,
                      tsc.fields['U'].vx, (x,y))
v = spinterp.griddata(np.vstack((tsc.x, tsc.y)).T,
                      tsc.fields['U'].vz, (x,y))
plt.streamplot(x,y,u,v, density=2,color='k',linewidth=0.5)



levels = np.linspace(0,800,17)
ticks = np.array([0, 200, 400, 600, 800])


plt.figure(r'Solar radiation')
#plt.tricontourf(tsc.x, tsc.y, tsc.fields['U'].Umag(),levels, cmap='viridis')
plt.tricontourf(tsc.x, tsc.y, tsc.fields['Rg'].Umag(),levels, cmap='viridis')
plt.plot(xc,yc,'--w',lw=2.5)

plt.axis('scaled')
plt.axis([-2,6,0,3])
plt.xlabel(r'$x/H$')
plt.ylabel(r'$z/H$')

cb = plotenv.colorbar(ticks,orientation='h',format='%g')
plotenv.add_text('(b)', color='k')
cb.set_label('$q_{r,sw} (W/m^2)$')

u = spinterp.griddata(np.vstack((tsc.x, tsc.y)).T,
                      tsc.fields['Rg'].vx, (x,y))
v = spinterp.griddata(np.vstack((tsc.x, tsc.y)).T,
                      tsc.fields['Rg'].vz, (x,y))

plt.quiver(x[::5,::5] , y[::5,::5], u[::5,::5], v[::5,::5], color='k')



"""


# ---------------------------------------------------------- Temperature
# levels = np.linspace(299,305,21)-273 #np.hstack((, 0.99,1.01, np.arange(1.1,2.1,0.1)))
# ticks  = np.linspace(299,305,11)-273 #np.hstack((np.arange(0,1,0.2), 1.0, np.arange(1.2,2.1,0.2)))
#levels = np.linspace(303,306,21)-273 #np.hstack((, 0.99,1.01, np.arange(1.1,2.1,0.1)))
#ticks  = np.linspace(303,306,11)-273 #np.hstack((np.arange(0,1,0.2), 1.0, np.arange(1.2,2.1,0.2)))

#levels = np.linspace(26,32,21)
#ticks = np.linspace(26,32,11)
levels = np.hstack((-3,np.arange(20)*0.3-2.85,3))
ticks  = np.hstack((-3,np.arange(9)[1::2]*0.3-2.85,0,2.85-np.arange(9)[1::2][::-1]*0.3,3))
#ticks = np.linspace(-3,3,11)

plt.figure(r'Temperature - T')
plt.tricontourf(tsc.x, tsc.y,tsc.triangles,tsc.fields['T'].s-303.85,#-273.15,
		        levels,cmap=palette['CMB']['DIV'])
plt.plot([-0.5,0.5,0.5,-.5,-.5],[0.5,0.5,1.5,1.5,0.5],'g',lw=1)
plt.axis('scaled')
plt.axis([-2,8,0,3])
plt.xlabel(r'$x/H$')
plt.ylabel(r'$y/H$')

plotenv.colorbar(ticks,orientation='h',format='%.1f')
plotenv.cleanupFigure(despine=False)
if saveFig: plt.savefig('surface_TField.png')

"""
