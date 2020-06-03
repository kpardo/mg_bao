## import packages
import numpy as np
import astropy.units as u
import astropy.constants as const
import pandas as pd
from scipy.interpolate import UnivariateSpline
from datetime import datetime
from hankel import SymmetricFourierTransform

from mg_bao.constants import *
## TODO: add errorbars to this analysis.
def import_powerspectra():
    ## import all powerspectra
    cambk, cambpkz0, cambpkz1100, cambpkz750 = np.loadtxt('../data/camb_pk.dat', usecols=(0,1,2,6),unpack=True)
    camb = pd.read_csv('../results/data_products/camb_pk.dat')
    cambk = camb['k'].to_numpy()
    cambpkz0 = camb['pkz0'].to_numpy()
    cambpkz1100 = camb['pkz1100'].to_numpy()

    sdssk, __, sdsspk, sdsspk_err  = np.loadtxt('../data/Beutler_2016_BAO/Beutleretal_pk_monopole_DR12_NGC_z1_postrecon_120.dat', unpack=True)
    sdssk *= boss_h
    sdsspk *= boss_h

    planck = pd.read_csv('../results/data_products/pb_z1100.dat')
    planckk = planck['k'].to_numpy()
    pk_z1100 = planck['pbz1100'].to_numpy()
    pk_z1100_u = planck['pbz1100_u'].to_numpy()
    pk_z1100_l = planck['pbz1100_l'].to_numpy()
    return planckk, pk_z1100, sdssk, sdsspk

def make_splines( planckk, pk_z1100, sdssk, sdsspk):
    sdss_spline = UnivariateSpline(sdssk, sdsspk, s=0., ext='zeros')
    ## use log10 of planckpk for spline because of large fluctuations
    log10planck_spline = UnivariateSpline(planckk, np.log10(pk_z1100), s=1., ext='zeros')
    return sdss_spline, log10planck_spline

def make_tk():
    planckk, pk_z1100, sdssk, sdsspk  = import_powerspectra()
    sdss_spline, log10planck_spline = make_splines(planckk, pk_z1100, sdssk,
            sdsspk)
    ks = np.linspace((lstar+0.5)/eta_star, np.max(sdssk), 1000)
    tk = sdss_spline(ks)/10**log10planck_spline(ks)
    results = np.array([ks, tk]).T
    table = pd.DataFrame(results, columns=['k', 'Tk'])
    filepath = '../results/data_products/transfer.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))

def create_r_array(ks):
    ws = ks/np.pi
    xmin = 0.5*1./ws[-1]
    xmax =  1./ws[0]
    deltax = 0.5*1./ws[-1]
    numx = 2*len(ws)
    r = np.linspace(xmin, xmax, numx)
    return r

def make_greens():
    ## get powerspectra to create spline again
    planckk, pk_z1100, sdssk, sdsspk  = import_powerspectra()
    ## create spline for tk
    sdss_spline = UnivariateSpline(sdssk, sdsspk, s=0., ext='zeros')
    ## use log10 of planckpk for spline because of large fluctuations
    log10planck_spline = UnivariateSpline(planckk, np.log10(pk_z1100), s=1., ext='zeros')
    ks = np.linspace((lstar+0.5)/eta_star, np.max(sdssk), 1000)
    tk = UnivariateSpline(ks,
            np.sqrt(sdss_spline(ks)/10**log10planck_spline(ks)),s=0.,
            ext='zeros')
    ## do the fourier transform with the help of Hankel
    ft = SymmetricFourierTransform(ndim=3, N = 200, h = 0.03)
    r = create_r_array(ks)
    Gr = ft.transform(tk,r, ret_err=False, inverse=True)
    ## save data
    results = np.array([r, Gr]).T
    table = pd.DataFrame(results, columns=['r', 'Gr'])
    filepath = '../results/data_products/greens.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))
