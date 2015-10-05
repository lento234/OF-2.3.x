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
z0 = flowProp['z0']
RH = flowProp['RH']
Tref = flowProp['RH']
qRef = flowProp['qRef']
ustar = flowProp['ustar']
kappa = flowProp['kappa']
kappa = flowProp['kappa']

vegProp  = ParsedParameterFile('../constant/vegetationProperties')
H = vegProp['H'][-1]
C = vegProp['C'][-1]
l = vegProp['l'][-1]
kc = vegProp['kc'][-1]
Cdf = vegProp['Cdf'][-1]
Rg0 = vegProp['Rg0'][-1]
cpa = vegProp['cpa'][-1]
Uref = 2.0
rhoa = vegProp['rhoa'][-1]
rsMin = vegProp['rsMin'][-1]
betaD = vegProp['betaD'][-1]
betaP = vegProp['betaP'][-1]
lmbda = vegProp['lambda'][-1]

transProp = ParsedParameterFile('../constant/transportProperties')
nu = transProp['nu'][-1]

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

# Extract coordinates
x,y,z = lc.lines[lnNames[0] + '_' + 'U'].xyz.T

# Extract data
U, V, W = lc.lines[lnNames[0] + '_' + 'U'].rawVars().T
Umid, Vmid, Wmid = 0.5*(U[1:]+U[:-1]), 0.5*(V[1:]+V[:-1]), 0.5*(W[1:]+W[:-1])
magU = np.sqrt(Umid**2 + Vmid**2 + Wmid**2)
LAD = lc.lines[lnNames[0] + '_' + 'LAD'].rawVars()
LAI = lc.lines[lnNames[0] + '_' + 'LAI'].rawVars()
Rg = lc.lines[lnNames[0] + '_' + 'Rg'].rawVars()[:,2]
Rn = lc.lines[lnNames[0] + '_' + 'Rn'].rawVars()
ra = lc.lines[lnNames[0] + '_' + 'ra'].rawVars()
rs = lc.lines[lnNames[0] + '_' + 'rs'].rawVars()
T =  lc.lines[lnNames[0] + '_' + 'T'].rawVars()
Tl =  lc.lines[lnNames[0] + '_' + 'Tl'].rawVars()
q = lc.lines[lnNames[0] + '_' + 'q'].rawVars()
E = lc.lines[lnNames[0] + '_' + 'E'].rawVars()
rhosat = lc.lines[lnNames[0] + '_' + 'rhosat'].rawVars()
qsat = lc.lines[lnNames[0] + '_' + 'qsat'].rawVars()
Ql = lc.lines[lnNames[0] + '_' + 'Ql'].rawVars()
Qs = lc.lines[lnNames[0] + '_' + 'Qs'].rawVars()

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

# Plot LAD
# plt.figure('LAD')
# plt.plot(LAD, z, '-', c=palette[0])
# plt.xlabel(r'$LAD$')
# plt.ylabel(r'$z$ [m]')
# plotEnv.cleanupFigure()

# Plot LAI
# plt.figure('LAI')
# plt.plot(LAI, z, '-', c=palette[0])
# plt.plot(LAD*(10.0-z), z, '--', c='r')
# plt.xlabel(r'$LAI$')
# plt.ylabel(r'$z$ [m]')
# plotEnv.cleanupFigure()

# Plot Rg
Rg2 = Rg0*np.ones(Rg.shape)
Rg2[z<=H] = Rg0*np.exp(-kc*LAD[z<=H]*(H-z[z<=H]))

# plt.figure('Rg')
# plt.semilogx(Rg, z, '-', c=palette[0])
# plt.semilogx(Rg2, z, '--', c='r')
# plt.xlabel(r'$Rg$')
# plt.ylabel(r'$z$ [m]')
# plotEnv.cleanupFigure()

# Plot Rn
Rn2 = (Rg2[1:]-Rg2[:-1])/(z[1:]-z[:-1])
zMid = 0.5*(z[1:]+z[:-1])

# plt.figure('Rn')
# plt.semilogx(Rn, z, '.-', c=palette[0])
# plt.semilogx(Rn2, zMid, '--', c='r')
# plt.xlabel(r'$Rn$')
# plt.ylabel(r'$z$ [m]')
# plotEnv.cleanupFigure()

# Plot aerodynamic resistance
ra2 = C*np.sqrt(l/magU)
plt.figure('ra')
plt.plot(ra,z,'.-',c=palette[0])
plt.plot(ra2,zMid,'--',c='r')
plt.xlabel(r'$r_a$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()

# Plot stomatal resistance
# rs2 = rs*(1+0.11*np.exp(0.34*(6.107)))
Tmid = 0.5*(T[1:]+T[:-1])
rs2 = rsMin*(31+Rn2)*(1+0.016*(Tmid-16.4-273.15)**2)/(6.7+Rn2)
rs2[(zMid>H) | (zMid<4)] = 0.
# rs3 = rsMin*(1 + 0.11*np.exp(0.34*s(6.10710**(7.5*T/(237.5+T))-1629*)))
# rs3[(z>H) | (z<4)] = 0.
plt.figure('rs')
plt.plot(rs,z,'.-',c=palette[0])
plt.plot(rs2,zMid,'.-',c='r')
# plt.plot(rs3,z,'.-',c='g')
plt.xlabel(r'$r_s$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()


# Temperature of air/leaf
# plt.figure('Tl-T')
# Ttemp = Tl.copy()
# Ttemp[Tl<0.001] = 293.15
# plt.plot(Ttemp-T,z,'.-',c=palette[0])
# plt.xlabel(r'$T_l$')
# plt.ylabel(r'$z$ [m]')
# plotEnv.cleanupFigure()


# Saturated vapour density
# Tc = Tmid-273.15
# rhosat2 = (5.018 + (0.32321)*Tc + (8.1847e-3)*Tc**2 + (3.1243e-4)*Tc**3)/1000.
rhosat2 = 0.0022*np.exp(77.3450+0.0057*Tmid-7235/Tmid)/(Tmid**9.2)
# rhosat2 = (6.335 + (0.6718)*Tc - (2.0887e-2)*Tc**2 + (7.3095e-4)*Tc**3)/1000.
#psat = np.exp(77.3450 + 0.0057*T - 7235/T)/(T**8.2)
#rhosat3 = 0.0022*psat/T
plt.figure('rhosat')
plt.plot(rhosat, z, '.-', c=palette[0])
plt.plot(rhosat2, zMid, '--', c='r')
plt.xlabel(r'$rhosat$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()


# saturated specific humidity
qsat2 = rhosat2/rhoa

plt.figure('qsat')
plt.plot(qsat, z, '.-', c=palette[0])
plt.plot(qsat2, zMid, '--', c='r')
plt.xlabel(r'$rhosat$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()


# Transpiration
LADmid = 0.5*(LAD[1:] + LAD[:-1])
qmid = 0.5*(q[1:] + q[:-1])
E2 = LADmid*rhoa*(qsat2-qmid)/(ra2+rs2)

plt.figure('E')
plt.plot(E, z, '.-', c=palette[0])
plt.plot(E2, zMid, '--', c='r')
plt.xlabel(r'$E$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()

# Latent heat
Ql2 = lmbda*E2

plt.figure('Ql')
plt.plot(Ql, z, '.-', c=palette[0])
plt.plot(Ql2, zMid, '--', c='r')
plt.xlabel(r'$Ql$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()


# leaf temperature
Tl2 = Tmid + (Rn2 - Ql2)*(ra2/(2.0*rhoa*cpa))

plt.figure('Tl')
plt.plot(Tl, z, '.-', c=palette[0])
plt.plot(Tl2, zMid, '.-', c='r')
plt.xlabel(r'$Tl$')
plt.xlabel(r'$z$')
plotEnv.cleanupFigure()

# sensible heat
Qs2 = 2.0*rhoa*cpa*LADmid*(Tl2-Tmid)/ra2

plt.figure('Qs')
plt.plot(Qs, z, '.-', c=palette[0])
plt.plot(Qs2, zMid, '--', c='r')
plt.xlabel(r'$Qs$')
plt.ylabel(r'$z$ [m]')
plotEnv.cleanupFigure()


"""

# ----------------
# P L O T   L A I
# ----------------

# ---------------------------
# P L O T   R A D I A T I O N
# ---------------------------
#
# Rg = lc.lines[lnNames[0] + '_' + 'Rg'].vz
# Rg2 = 100*np.exp(-0.75*LAI)
# Rg3 = (100*np.exp(-0.75*1*(10.0-z)))[:71]
# plt.figure('Rg')
# plt.plot(Rg, z, '-', c=palette[0])
# plt.plot(100*np.exp(-0.75*LAI), z, '.--', c='r')
# plt.plot((100*np.exp(-0.75*1*(10.0-z)))[:71], z[:71], '.--', c='g')
# curax = plt.axis()
# plt.xlabel(r'$Rg(z)$')
# plt.ylabel(r'$z$ [m]')
# plotEnv.cleanupFigure()
#
# Rn = lc.lines[lnNames[0] + '_' + 'Rn'].s
# dz = z[1:]-z[:-1]
# zMid = 0.5*(z[1:]+z[:-1])
#
# plt.figure('Rn')
# plt.plot(Rn, z, '-', c=palette[0])
# plt.semilogx((Rg[1:]-Rg[:-1])/dz, zMid, '--', c='r')
# plt.semilogx((Rg2[1:]-Rg2[:-1])/dz, zMid, ':', c='g')
# plt.semilogx((Rg3[1:]-Rg3[:-1])/dz[:70], zMid[:70], '--', c='b')
# curax = plt.axis()
# plt.xlabel(r'$Rn$')
# plt.ylabel(r'$z$ [m]')
# plotEnv.cleanupFigure()
#
# # ----------------
# # temperature
# # ----------------
#
# T = lc.lines[lnNames[0]+'_'+'T'].s
# plt.figure('T')
# plt.plot(T, z, c=palette[0])
# plt.xlabel(r'$T$')
# plt.ylabel(r'$z$')
# plotEnv.cleanupFigure()
#
#
# rhosat = lc.lines[lnNames[0]+'_'+'rhosat'].s
# plt.figure('rhosat')
# plt.plot(rhosat, z, c=palette[0])
# plt.plot(np.exp(), z, c='r')
# plt.xlabel(r'$\rho*$')
# plt.ylabel(r'$z$')
# plotEnv.cleanupFigure()
"""

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
