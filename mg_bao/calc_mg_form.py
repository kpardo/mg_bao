## import packages
import numpy as np
import astropy.units as u
import astropy.constants as const
import pandas as pd
from scipy.interpolate import UnivariateSpline
from datetime import datetime
from hankel import SymmetricFourierTransform, get_h

from mg_bao.constants import *
from mg_bao.convenience import *

## TODO: add errorbars to this analysis.
def import_powerspectra(lerr = False, uerr = False):
    ## import all powerspectra
    cambk, cambpkz0, cambpkz1100, cambpkz750 = np.loadtxt('../data/camb_pk.dat', usecols=(0,1,2,6),unpack=True)
    camb = pd.read_csv('../results/data_products/camb_pk.dat')
    cambk = camb['k'].to_numpy()
    cambpkz0 = camb['pkz0'].to_numpy()
    cambpkz1100 = camb['pkz1100'].to_numpy()

    sdssk, __, sdsspk, sdsspk_err  = np.loadtxt('../data/Beutler_2016_BAO/Beutleretal_pk_monopole_DR12_NGC_z1_postrecon_120.dat', unpack=True)
    sdssk *= boss_h
    sdsspk *= boss_h
    sdsspk *= sdsspk_err

    planck = pd.read_csv('../results/data_products/pb_z1100.dat')
    planckk = planck['k'].to_numpy()
    pk_z1100 = planck['pbz1100'].to_numpy()
    pk_z1100_u = planck['pbz1100_u'].to_numpy()
    pk_z1100_l = planck['pbz1100_l'].to_numpy()

    if lerr:
        return planckk, pk_z1100_l, sdssk, sdsspk-sdsspk_err

    if uerr:
        return planckk, pk_z1100_u, sdssk, sdsspk+sdsspk_err
    return planckk, pk_z1100, sdssk, sdsspk

def make_splines( planckk, pk_z1100, sdssk, sdsspk):
    sdss_spline = UnivariateSpline(sdssk, sdsspk, s=0., ext='zeros')
    ## use log10 of planckpk for spline because of large fluctuations
    log10planck_spline = UnivariateSpline(planckk, np.log10(pk_z1100), s=0., ext='zeros')
    return sdss_spline, log10planck_spline

def make_tk():
    planckk, pk_z1100, sdssk, sdsspk  = import_powerspectra()
    sdss_spline, log10planck_spline = make_splines(planckk, pk_z1100, sdssk,
            sdsspk)
    planckk, pk_z1100_u, sdssk, sdsspk_u = import_powerspectra(uerr=True)
    pk_z1100_u[pk_z1100_u == 0] = 1.e-32 ## make it some tiny number so no nan.
    sdss_spline_u, log10planck_spline_u = make_splines(planckk, pk_z1100_u,
            sdssk, sdsspk_u)
    planckk, pk_z1100_l, sdssk, sdsspk_l = import_powerspectra(lerr=True)
    sdss_spline_l, log10planck_spline_l = make_splines(planckk, pk_z1100_l,
            sdssk, sdsspk_l)
    ks = np.linspace((lstar+0.5)/eta_star, np.max(sdssk), 1000)
    tk = sdss_spline(ks)/10**log10planck_spline(ks)
    tk_l = sdss_spline_l(ks)/10**log10planck_spline_l(ks)
    tk_u = sdss_spline_u(ks)/10**log10planck_spline_u(ks)
    results = np.array([ks, tk, tk_l, tk_u]).T
    table = pd.DataFrame(results, columns=['k', 'Tk', 'Tk_l', 'Tk_u'])
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

def make_greens(ext='zeros'):
    ## get powerspectra to create spline again
    planckk, pk_z1100, sdssk, sdsspk  = import_powerspectra()
    __, pk_z1100_l, ___, sdsspk_l  = import_powerspectra(lerr=True)
    __, pk_z1100_u, ___, sdsspk_u  = import_powerspectra(uerr=True)
    pk_z1100_u[pk_z1100_u == 0] = 1.e-32 ## make it some tiny number so no nan.
    ## create spline for tk
    sdss_spline = UnivariateSpline(sdssk, sdsspk, s=0., ext='zeros')
    sdss_spline_l = UnivariateSpline(sdssk, sdsspk_l, s=0., ext='zeros')
    sdss_spline_u = UnivariateSpline(sdssk, sdsspk_u, s=0., ext='zeros')
    ## use log10 of planckpk for spline because of large fluctuations
    log10planck_spline = UnivariateSpline(planckk, np.log10(pk_z1100), s=0., ext='zeros')
    log10planck_spline_l = UnivariateSpline(planckk, np.log10(pk_z1100_l), s=0., ext='zeros')
    log10planck_spline_u = UnivariateSpline(planckk, np.log10(pk_z1100_u), s=0., ext='zeros')
    ks = np.linspace((lstar+0.5)/eta_star, np.max(sdssk), 1000)
    tk = UnivariateSpline(ks,
            np.sqrt(sdss_spline(ks)/10**log10planck_spline(ks)),s=0.,
            ext=ext)
    tk_l = UnivariateSpline(ks,
            np.sqrt(sdss_spline_l(ks)/10**log10planck_spline_l(ks)),s=0.,
            ext=ext)
    tk_u = UnivariateSpline(ks,
            np.sqrt(sdss_spline_u(ks)/10**log10planck_spline_u(ks)),s=0.,
            ext=ext)
    ## do the fourier transform with the help of Hankel
    ## first create the r array using sdss k -- more conservative
    r = create_r_array(sdssk)
    ## find optimal parameters for hankel to use
    deltah, err, N = get_h(tk, nu=3, K=[np.min(r), np.max(r)],cls=SymmetricFourierTransform, inverse=True)
    if np.any(np.abs(err) > 1.e-2):
        print(err)
        print("The error on the FT is high (> 1\%). You should check this!")
    ft = SymmetricFourierTransform(ndim=3, N = N, h = deltah)
    Gr = ft.transform(tk,r, ret_err=False, inverse=True)
    Gr_l = ft.transform(tk_l,r, ret_err=False, inverse=True)
    Gr_u = ft.transform(tk_u,r, ret_err=False, inverse=True)
    ## save data
    results = np.array([r, Gr, Gr_l, Gr_u]).T
    table = pd.DataFrame(results, columns=['r', 'Gr', 'Gr_l', 'Gr_u'])
    filepath = '../results/data_products/greens_'+ext+'.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))
