#
# Code to calculate Rossby stationary wave number, following Hoskins and Ambrizzi 1993.
# Rossby waves are refracted in latitude towards regions of higher Ks
# As dl/dt = cg dKs/dy
# Ks is defined as (Beta*/Uzm)^0.5
# Or on a Mercator coordinate:
# Ks = (a BetaM/vbar)^0.5
# where a is Earth radius, vbar = Ubar/acos(phi) and phi is latitue
# BetaM = [2*omega - (1/cos(phi) d/dphi)(1/cos(phi) d/dphi(cos2(phi)vbar)]*cos2(phi)/a
#
# Author: Rachel White, rachel.white@cantab.net
# Created: 7-Dec-2016
#
import os, errno
import netCDF4
import numpy as np
import math
import datetime as dt
import xray
from rhwhitepackages.readwrite import xrayOpen

# define level you want to calculate Ks at
ilev = 250

# pre-defined variables, including file names and directories
startyr=1979
endyr=2015
nyears = endyr-startyr+1

Dir = '/home/disk/eos4/rachel/Obs/ERAI/'
filename = ('ERAI_monthlymeans_DJF_' + str(startyr) + '-' + str(endyr) + '.nc')

# constants
a = 6.37122e06  # radius of Earth
a2 = math.pow(a,2)
omega =  7.2921e-5  # rotation rate of earth
g = 9.80616

datain = xrayOpen(Dir + filename)
uIn = datain['u'].sel(level=ilev)   # read in U at level ilev
latsIn = uIn.coords['latitude']
lonsIn = uIn.coords['longitude']
if lonsIn[-1] - lonsIn[1] < 358:
    exit('caution, not a global grid, and this could cause things to break')
print uIn.shape
# Calculate trigonometric functions
phi = np.radians(latsIn)
cphi = np.cos(phi)
sphi = np.sin(phi)
c2phi = np.power(cphi,2)
acphi = a * cphi
asphi = a * sphi 
# Calculate coriolis parameter f
f = 2*omega*sphi
f2 = np.power(f,2)

def ddphi(arrayin,phiIn):
    if len(arrayin.shape) == 3:   # check we still have 3 dimensions, so DT, DX and DY correspond as expected
        if arrayin.shape[1] < arrayin.shape[2]:     # check which way round lats and lons are.
            # caution! this may be wrong if it's not a global grid
            darrayDT,darrayDY, darrayDX = np.gradient(arrayin)
            dphi = np.gradient(phiIn)
            return(darrayDY/dphi[None,:,None])  # using broadcasting to extend dphi array dimensions

# Calculate terms in:
# Ks = (a BetaM/vbar)^0.5
# where a is Earth radius, vbar = Ubar/acos(phi) and phi is latitue
# BetaM = [2*omega - {(1/cos(phi) d/dphi)(1/cos(phi) d/dphi(cos2(phi)vbar))}]*cos2(phi)/a
#
# Using inbuilt numpy broadcasting where possible

# Calculate vbar
vbar = uIn / (a * cphi) 


# Calculate c2phivbar = cos2phi * U / (cos(phi) a) = U cos(phi)/a
c2phivbar = uIn * cphi / a

# calculate gradient wrt latitude
dc2phivbarDPHI = ddphi(c2phivbar,phi)

# Multiply by 1/cphi
dphi = np.gradient(cphi)
workingcalc = dc2phivbarDPHI / cphi.values[None,:,None]

# calcuate gradient wrt latitude phi
workingcalcDPHI = ddphi(workingcalc,phi)

# multiply by 1/cphi
workingcalcfinal = workingcalcDPHI / cphi.values[None,:,None]

# calculate betaM = (2* omega - workingcalcfinal) * c2phi/a
betaM = (2*omega - workingcalcfinal) * c2phi.values[None,:,None]/a
# calculate Ks = (a * betaM/ vbar)^0.5
ks2 = a * betaM /vbar
ks = np.sqrt(ks2)
print np.amax(ks)
# Create xarray dataarray from this

# write out to netcdf

