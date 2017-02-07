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
from pyFlowStat.old import TriSurfaceContainer as TC

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

x,y = tsc.x, tsc.y
size = x.shape

# data
v = tsc.fields['U'].Umag() # velocity,  m/s
t = tsc.fields['T'].s - 273.15 # temperature, celsius
f = tsc.fields['RH'].s # relative humidity, %
Kglob = tsc.fields['Rg'].Umag() # global radiation
hSl = np.ones(size)*90.0 # solar altitude, angle

# save data
np.savetxt('importdata_bioklima.dat',  np.vstack((v,t,f,Kglob,hSl)).T,
            delimiter='\t', header='v\tt\tf\tKglob\thSl',comments='')
