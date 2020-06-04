'''
plotting functions
'''


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from scipy.ndimage import gaussian_filter1d
from datetime import datetime

from mg_bao.constants import *
from mg_bao.convenience import *

## set fig params
sns.set_context("paper")
sns.set_style('ticks')
sns.set_palette('colorblind')
figparams = {
        'text.latex.preamble': [r'\usepackage{amsmath}'],
        'text.usetex':True,
        'axes.labelsize':20.,
        'xtick.labelsize':16,
        'ytick.labelsize':16,
        'figure.figsize':[10., 8.],
        'font.family':'DejaVu Sans',
        'legend.fontsize':18}
plt.rcParams.update(figparams)
cs = plt.rcParams['axes.prop_cycle'].by_key()['color']

def savefig(fig, figpath, writepdf=False, dpi=450):
    ## stolen from luke
    fig.savefig(figpath, dpi=dpi, bbox_inches='tight')
    print('{}: made {}'.format(datetime.now().isoformat(), figpath))

    if writepdf:
        pdffigpath = figpath.replace('.png','.pdf')
        fig.savefig(pdffigpath, bbox_inches='tight', rasterized=True, dpi=dpi)
        print('{}: made {}'.format(datetime.now().isoformat(), pdffigpath))

    plt.close('all')

def pk_plot(planckk, pkz1100, pkz1100yerr, sdssk, sdsspk, sdsspk_err, cambk,
        cambpkz0, cambpkz1100, figpath):
    '''
    makes power spectrum plot
    '''
    dimvalue = 5.e-8 # random number that scales down Pk,z=1100 from data
    upvalue = 3.e5 # random number that scales up camb Pk, z=1100`
    f = plt.figure()
    plt.errorbar(sdssk, sdsspk, yerr=sdsspk_err, fmt='o', color='black', label=r'$P_{bb}(k, z=0.38)$')
    plt.errorbar(planckk[6:], pkz1100[6:]*dimvalue,
            yerr=np.array(pkz1100yerr)[:, 6:]*dimvalue, fmt='o',
                 c='black', mfc='white', label=r'$P_{bb}(k, z=1100)$')
    plt.axvline(lstar/eta_star, linestyle='dashed', color='black', linewidth=3)
    plt.plot(cambk, cambpkz0, color=cs[0],linewidth=2, linestyle='dotted')
    plt.plot(cambk, cambpkz1100*upvalue, color=cs[0], linewidth=2, linestyle='dotted', label='CAMB, z=1100')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim([2e-3, 1e5])
    plt.xlim([0.009, 0.14])
    plt.legend(fontsize=12, loc='lower left')
    plt.xlabel(r'$k~[\rm{Mpc}^{-1}]$')
    plt.ylabel(r'$P_{bb}(k)~[\rm{Mpc}^3]$');
    savefig(f, figpath, writepdf=True)

def tk_plot(ks, tk, cambk, camb_pkdiv, filepath):
    norm = 1./tk[0] ##first non-zero sdss spline value
    cambnorm = 1.e-7
    f = plt.figure()
    plt.plot(ks, tk*norm, c='black', linewidth=2, label='Data')
    plt.plot(cambk, camb_pkdiv*cambnorm, color=cs[0],linewidth=2, linestyle='dotted',label='CAMB')
    plt.axvline(lstar/eta_star, linestyle='dashed', color='black', linewidth=3)
    plt.legend()
    plt.yscale('log')
    plt.xlim([0.009, 0.105])
    plt.ylim([1e-1, 1e5])
    plt.xlabel(r'$k~[\rm{Mpc}^{-1}]$')
    plt.ylabel(r'$T^2(k)$');
    savefig(f, filepath, writepdf=True)

def greens_plot(rs, Gr,rs2, Gr2, filepath):
    Gr = gaussian_filter1d(Gr/Gr[0], 1) ## normalize to 1 at small scales and smooth
    Gr2 = gaussian_filter1d(Gr2/Gr2[0], 1)
    f = plt.figure()
    plt.plot(rs, Gr, c='black', linewidth=2, label='Zero Extrapolation')
    plt.plot(rs2, Gr2, c='black', linewidth=2, linestyle='dotted', label='Constant Value Extrapolation')
    plt.legend()
    plt.xlabel(r'$r~[\rm{Mpc}]$')
    plt.ylabel(r'$\mathcal{G}(r)$')
    savefig(f, filepath, writepdf=True)
