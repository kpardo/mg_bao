## import packages
import numpy as np
import sys, platform, os
import matplotlib.pyplot as plt
import seaborn as sns
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import binned_statistic

from mg_bao import *

## set plot settings
sns.set_context("paper")
sns.set_style('ticks')
sns.set_palette('colorblind')
figparams = {
        'text.latex.preamble': [r'\usepackage{amsmath}'],
        'text.usetex':True,
        'axes.labelsize':20.,
        'xtick.labelsize':16,
        'ytick.labelsize':16,
        'figure.figsize':[8., 6.],
        'font.family':'DejaVu Sans',
        'legend.fontsize':12}
plt.rcParams.update(figparams)
cs = plt.rcParams['axes.prop_cycle'].by_key()['color']

## set cosmology to Planck 2018 Paper I Table 6
cosmo = FlatLambdaCDM(H0=67.32, Om0=0.3158, Ob0=0.03324)

## load data
planckl, planckdl, lerr, herr = np.loadtxt('../data/planck_ps/COM_PowerSpect_CMB-EE-full_R3.01.txt', unpack=True)
actl, actcl, acterr = np.loadtxt('../data/act_ps/cl_cmb_ee.dat', unpack=True)

## change planck data to cl from dl
planckcl = planckdl/(planckl*(planckl+1))*2*np.pi
lclerr = lerr/(planckl*(planckl+1))*2*np.pi
hclerr = herr/(planckl*(planckl+1))*2*np.pi

## mask data
lmask1 = planckl < 2026.
lmask2 = planckl > 24
lmask = np.array(lmask1 & lmask2)
planckcl = planckcl[lmask]
planckl = planckl[lmask]
hclerr = hclerr[lmask]
lclerr = lclerr[lmask]
actlmask = actl < 2001.
actcl = actcl[actlmask]
acterr = acterr[actlmask]
actl = actl[actlmask]

## bin planck data so it has the same mid bin values as act data
bins = np.arange(25,np.max(actl)+25, 50)
clbin, lbinedges, __ = binned_statistic(planckl, planckcl, bins=bins)
clupp, __, __ = binned_statistic(planckl, planckcl+hclerr, bins=bins)
cldown, __, __ = binned_statistic(planckl, planckcl-lclerr, bins=bins)
lbinmid = np.array([(a + b) /2 for a,b in zip(lbinedges[:-1], lbinedges[1:])])

## combine data via quadrature where they overlap.
newlmask = lbinmid >= np.min(actl)
cls = clbin[~newlmask]
cls = np.append(cls, (clbin[newlmask] + actcl[:-1])/2.)
cl_upp = clupp[~newlmask]
cl_upp = np.append(cl_upp, (clupp[newlmask]+actcl[:-1]+acterr[:-1])/2.)
cl_down = cldown[~newlmask]
cl_down = np.append(cl_down, (cldown[newlmask]+ actcl[:-1] - acterr[:-1])/2.)

## convert to pk and ks
# eta_star = cosmo.comoving_distance(1100).value
eta_star = cosmo.comoving_distance(1059.94).value ## z_drag from Planck 2018 cosmology paper Table 2, all Planck alone
planckk = (lbinmid+0.5)/(eta_star)
pkee = cls*lbinmid**2/planckk**3*np.pi
pkee_u = cl_upp*lbinmid**2/planckk**3*np.pi
pkee_l = cl_down*lbinmid**2/planckk**3*np.pi

## print to file
results = np.array([planckk, pkee, pkee_u, pkee_l]).T
np.savetxt('../results/data/pk_ee.dat', results)

## check by making figures
plt.figure()
plt.errorbar(planckl[::5], planckcl[::5], yerr=[lclerr[::5], hclerr[::5]],fmt='o', c=cs[0], label=r'Planck Cls')
plt.errorbar(actl, actcl, yerr=acterr, c=cs[1],fmt='o',label=r'ACT Cls')
plt.errorbar(lbinmid, cls, yerr=[cls-cl_down, cl_upp-cls],fmt='o', c='black', label=r'Combined', zorder=3)
plt.yscale('log')
plt.legend()
plt.xlabel(r'$l$')
plt.ylabel(r'$C_l$')
plt.tight_layout()
plt.savefig('../results/planck_cls.png')

plt.figure()
plt.errorbar(planckk, pkee, yerr=[pkee-pkee_l, pkee_u-pkee], fmt='o', color='black')
plt.yscale('log')
plt.xlabel(r'$k~[\rm{Mpc}^{-1}]$')
plt.ylabel(r'$P_{EE}(k)~[\rm{Mpc}^3]$');
plt.savefig('../results/planck_pkee.png')
