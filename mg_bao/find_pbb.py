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
from scipy.stats import norm
from scipy.optimize import nnls


from mg_bao.convenience import *

## set constants
h = 0.676 ##BOSS value
eta_star = cosmo.comoving_distance(1059.94).value ## z_drag from Planck 2018 cosmology paper Table 2, all Planck alone
# thetamc = 1.0409/100. ##in rad., Planck 2015(?) value.
# rs = thetamc*eta_star
rs = 147.09 ## try rs=r_drag from Planck 2018 same table as z_drag above
# lstar = np.pi/thetamc
dklss = np.pi/19. ##width of last scattering

## define functions
def tik_lsq_reg(A_,b_,l_,x0=0):
    ## use a vector normalization "constant" because range of y values is huge
    lvec = np.logspace(np.log10(l_), np.log10(l_)+0.7, len(b_))
    C = np.vstack([A_, np.diagflat([lvec])])
    d = np.hstack([b_, np.zeros((len(b_)))])
    sol, res = nnls(C,d)
    return sol, res

def add_thickness_lss_effect(k_, pk_, dklss_):
    return pk_*np.exp(k_/dklss_)

## first, get CAMB spectra
cambk, cambpkz0, cambpkz1100 = np.loadtxt('../data/camb_pk.dat',usecols=(0,1,2), unpack=True)

## load data
planckk, pkee, pkee_u, pkee_l = np.loadtxt('../results/data/pk_ee.dat', unpack=True)
sdssk, __, sdsspk, sdsspk_err  = np.loadtxt('../data/Beutler_2016_BAO/Beutleretal_pk_monopole_DR12_NGC_z1_postrecon_120.dat', unpack=True)
sdssk *= h
sdsspk *= h

## transform pkee -> pkbb (w/ regularization)
lambdaa = 1.8e-2
AA = np.identity((len(planckk)))*(np.sin(planckk*rs))**2
Aksq_tik = tik_lsq_reg(AA, pkee/planckk**2, lambdaa)[0]
pk_z1100_tik = Aksq_tik*np.cos(rs*planckk)**2

Aksq_tik_u = tik_lsq_reg(AA, pkee_u/planckk**2, lambdaa)[0]
pk_z1100_tik_u = Aksq_tik_u*np.cos(rs*planckk)**2

Aksq_tik_l = tik_lsq_reg(AA, pkee_l/planckk**2, lambdaa)[0]
pk_z1100_tik_l = Aksq_tik_l*np.cos(rs*planckk)**2

## add thickness of last scatter effect
pk_z1100 =  add_thickness_lss_effect(planckk, pk_z1100_tik, dklss)
pk_z1100_u =  add_thickness_lss_effect(planckk, pk_z1100_tik_u, dklss)
pk_z1100_l =  add_thickness_lss_effect(planckk, pk_z1100_tik_l, dklss)


pk1100_lerr = pk_z1100 - pk_z1100_l
pk1100_uerr = pk_z1100_u - pk_z1100

results = np.array([planckk, pk_z1100, pk_z1100_l, pk_z1100_u]).T
np.savetxt('../results/pb_z1100.dat', results)
np.savetxt('../results/aksq.dat', Aksq_tik)

# dimvalue = 1.e-8
# pk_z1100*=dimvalue
# pk1100_lerr*=dimvalue
# pk1100_uerr*=dimvalue
#

# cambpkz1100 *= cambpkz1100*1.e8
# ## plot figure
# ks = np.linspace(np.min(sdssk), np.max(sdssk), 1000)
# plt.figure()
# plt.errorbar(sdssk, sdsspk, yerr=sdsspk_err, fmt='o', color='black', label=r'$P_{bb}(k, z=0.38)$')
# # plt.plot(ks, sdss_spline(ks), color='black')
# plt.errorbar(planckk, pk_z1100, yerr=[pk1100_lerr, pk1100_uerr], fmt='o', label=r'$P_{bb}(k, z=1100)$')
# plt.axvline((lstar+0.5)/eta_star, linestyle='dashed', color='black', linewidth=3)
# plt.plot(cambk, cambpkz0, color=cs[1], linestyle='dotted')
# plt.plot(cambk, cambpkz1100, color=cs[1], linestyle='dotted', label='CAMB')
# plt.yscale('log')
# plt.ylim([2e-3, 1e6])
# plt.xlim([0.005, 0.135])
# plt.legend(fontsize=12, loc='upper right')
# plt.xlabel(r'$k~[\rm{Mpc}^{-1}]$')
# plt.ylabel(r'$P_{bb}(k)~[\rm{Mpc}^3]$');
# plt.tight_layout();
# plt.show()
