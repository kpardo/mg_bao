## import packages
import numpy as np
import astropy.units as u
import astropy.constants as const
import pandas as pd
from scipy.interpolate import UnivariateSpline
from datetime import datetime
from hankel import SymmetricFourierTransform, get_h
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel

from mg_bao.constants import *
from mg_bao.convenience import *
from mg_bao.paths import DATADIR, RESULTSDIR
from mg_bao.plotting import tk_for_ref

def import_powerspectra(lerr = False, uerr = False, err=False):
    ## import all powerspectra
    camb = pd.read_csv(RESULTSDIR+'data_products/camb_pk.dat')
    cambk = camb['k'].to_numpy()
    cambpkz0 = camb['pkz0'].to_numpy()
    cambpkz1100 = camb['pkz1100'].to_numpy()

    sdssk, __, sdsspk, sdsspk_err  = np.loadtxt(DATADIR+'Beutler_2016_BAO/Beutleretal_pk_monopole_DR12_NGC_z1_postrecon_120.dat', unpack=True)
    sdssk *= boss_h
    sdsspk *= boss_h
    sdsspk_err *= boss_h

    planck = pd.read_csv(RESULTSDIR+'data_products/pb_z1100.dat')
    planckk = planck['k'].to_numpy()
    pk_z1100 = planck['pbz1100'].to_numpy()
    pk_z1100_u = planck['pbz1100_u'].to_numpy()
    pk_z1100_l = planck['pbz1100_l'].to_numpy()

    if lerr:
        return planckk, pk_z1100_l, sdssk, sdsspk - sdsspk_err

    if uerr:
        return planckk, pk_z1100_u, sdssk, sdsspk + sdsspk_err
    if err:
        return planckk, pk_z1100, sdssk, sdsspk, sdsspk_err
    return planckk, pk_z1100, sdssk, sdsspk

def make_splines( planckk, pk_z1100, sdssk, sdsspk):
    sdss_spline = UnivariateSpline(sdssk, sdsspk, s=0., ext='zeros')
    ## use log10 of planckpk for spline because of large fluctuations
    log10planck_spline = UnivariateSpline(planckk, np.log10(pk_z1100), s=0., ext='zeros')
    return sdss_spline, log10planck_spline

def get_tk_err(ks, pk, sdsspk, pk_u, pk_l, sdsserr):
    ## first get arrays from splines
    tk = sdsspk(ks)/10**(pk(ks))
    puerr = 10**(pk_u(ks)) -10**(pk(ks))
    plerr = 10**(pk(ks))-10**(pk_l(ks))
    uerr = err_div(tk, 10**(pk(ks)), sdsspk(ks), puerr, sdsserr(ks))
    lerr = err_div(tk,10**(pk(ks)), sdsspk(ks), plerr, sdsserr(ks))
    return tk+uerr, tk-lerr

def err_div(f, a,b,siga, sigb):
    '''
    calculates error for division of two variables
    assumes they are independent
    '''
    aterm = (siga/a)**2
    bterm = (sigb/b)**2
    return np.abs(f)*np.sqrt(aterm + bterm)

def make_tk():
    planckk, pk_z1100, sdssk, sdsspk, sdsserr  = import_powerspectra(err=True)
    sdsserr_spline = UnivariateSpline(sdssk, sdsserr, s=0., ext='zeros')
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
    tk_u, tk_l = get_tk_err(ks,log10planck_spline, sdss_spline,
            log10planck_spline_u, log10planck_spline_l, sdsserr_spline)
    results = np.array([ks, tk, tk_l, tk_u]).T
    table = pd.DataFrame(results, columns=['k', 'Tk', 'Tk_l', 'Tk_u'])
    filepath = RESULTSDIR+'data_products/transfer.dat'
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

def get_gp_cov(k, tk, norm, cambspline):
    X = k[:, np.newaxis]
    y = np.log10(tk*norm) - cambspline(k)
    kernel = RBF() +WhiteKernel()
    gp = GaussianProcessRegressor(kernel=kernel, alpha=2., normalize_y=False,
            n_restarts_optimizer=50).fit(X,y)
    X_ = np.logspace(-5, 1, 2000)
    y_mean, y_cov = gp.predict(X_[:, np.newaxis], return_cov=True)
    return X_, y_mean, y_cov

def splinefortransform(x, logspline, maxk=5.):
    s = 10**logspline(x)
    s[x>maxk] = 0.
    return s

def do_FT(r, k_, logkfunc):
    logsamplespline = UnivariateSpline(k_, logkfunc, s=0.,ext='zeros')
    samplespline = lambda x: np.sqrt(splinefortransform(x, logsamplespline))
    deltah, result, N = get_h(samplespline, nu=3, K=[np.min(r), np.max(r)],
            cls=SymmetricFourierTransform, inverse=True, maxiter=25)
    ft = SymmetricFourierTransform(ndim=3, N = N, h = deltah)
    Gr = ft.transform(samplespline,r, ret_err=False, inverse=True)
    return Gr

def make_greens():
    ## get tk
    try:
        tkdat = pd.read_csv(RESULTSDIR+'data_products/transfer.dat')
        k = tkdat['k']
        tk = tkdat['Tk']
    except:
        make_tk()
    tknorm = 1./tk[0] ##first non-zero tk spline value
    ## get sdssk for rarray later.
    sdssk, __, __, __  = np.loadtxt(DATADIR+'Beutler_2016_BAO/Beutleretal_pk_monopole_DR12_NGC_z1_postrecon_120.dat', unpack=True)
    ## get CAMB tk for gp
    camb = pd.read_csv(RESULTSDIR+'data_products/camb_pk.dat')
    cambk = camb['k'].to_numpy()
    cambpkz0 = camb['pkz0'].to_numpy()
    cambpkz1100 = camb['pkz1100'].to_numpy()
    cambpkdiv = cambpkz0/cambpkz1100
    cambnorm = 1.e-7 ## sets it about equal to tk so we can splice
    cambspline = UnivariateSpline(cambk, np.log10(cambpkdiv*cambnorm), s=0.)
    ## get gp results -- these are in log space
    X_, y_mean, y_cov = get_gp_cov(k, tk, tknorm, cambspline)
    ## make tk plot to check gp results
    tk_for_ref(k[:, np.newaxis],y_mean, y_cov, tk, tknorm, X_, cambspline)
    ## get real space values
    kforgreen = np.linspace(np.min(sdssk), np.max(sdssk), 100)
    r = create_r_array(kforgreen)
    ##thin out because noise.
    ## do transform 
    Gr = do_FT(r, X_, y_mean + cambspline(X_))
    Gr_l = do_FT(r, X_, y_mean - np.sqrt(np.diag(y_cov)) + cambspline(X_))
    Gr_u = do_FT(r, X_, y_mean + np.sqrt(np.diag(y_cov)) + cambspline(X_))

    ## save data
    results = np.array([r, Gr, Gr_l, Gr_u]).T
    table = pd.DataFrame(results, columns=['r', 'Gr', 'Gr_l', 'Gr_u'])
    filepath = RESULTSDIR+'data_products/greens.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))


