## import packages
import numpy as np
import sys, platform, os
import matplotlib.pyplot as plt
import seaborn as sns
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck15 as cosmo
import camb
from camb import model, initialpower
from scipy.interpolate import UnivariateSpline

## set constants
h = 0.676 ##BOSS value
# eta_star = cosmo.comoving_distance(1100).value
eta_star = cosmo.comoving_distance(1059.94).value ## z_drag from Planck 2018 cosmology paper Table 2, all Planck alone
# thetamc = 1.0409/100. ##in rad., Planck 2015(?) value.
# rs = thetamc*eta_star
rs = 147.09
# lstar = np.pi/thetamc
lstar = np.pi*eta_star/rs
dklss = np.pi/19. ##width of last scattering

## import all powerspectra
cambk, cambpkz0, cambpkz1100, cambkt, camb_tk, primordialpk = np.loadtxt('../data/camb_pk.dat', unpack=True)
sdssk, __, sdsspk, sdsspk_err  = np.loadtxt('../data/Beutler_2016_BAO/Beutleretal_pk_monopole_DR12_NGC_z1_postrecon_120.dat', unpack=True)
sdssk *= h
sdsspk *= h

planckk, pk_z1100, pk_z1100_l, pk_z1100_u = np.loadtxt('../results/pb_z1100.dat', unpack=True)

## make splines
sdss_spline = UnivariateSpline(sdssk, sdsspk, s=0., ext='zeros')
planck_spline = UnivariateSpline(planckk, np.log10(pk_z1100), s=1., ext='zeros')

## tk
ks = np.linspace((lstar+0.5)/eta_star, np.max(sdssk), 1000)
tk = sdss_spline(ks)/10**planck_spline(ks)
results = np.array([ks, tk]).T
np.savetxt('../results/transfer.dat', results)

## ak
tk_spline = UnivariateSpline(ks, tk, s=0.0, ext=1)
cambtk_spline = UnivariateSpline(cambkt, camb_tk*cambkt**2, s=0.0, ext=1)
primordialpk_spline = UnivariateSpline(cambkt, primordialpk, s=0.0, ext=1)
ak = -1.j*ks*tk_spline(ks)/tk_spline(ks[0])*cambtk_spline(ks[0])*primordialpk_spline(ks)
aGR = -1.j*ks*cambtk_spline(ks)*primordialpk_spline(ks)
results = np.array([ks, ak, aGR]).T
np.savetxt('../results/accel_k.dat', results)
