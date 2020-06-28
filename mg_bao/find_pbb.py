## import packages
import numpy as np
import astropy.units as u
import astropy.constants as const
import pandas as pd
from scipy.optimize import nnls
from datetime import datetime


from mg_bao.convenience import *
from mg_bao.constants import *
from mg_bao.paths import RESULTSDIR

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

def calc_pbb(k, pk, lss, lambdaa = 1.8e-2):
    ## lambdaa is the regularization constant.
    ## it was found using L-curve optimization
    AA = np.identity(len(k))*(np.sin(k*rs))**2
    Aksq_tik = tik_lsq_reg(AA, pk/k**2, lambdaa)[0]
    pbb0 = Aksq_tik*np.cos(rs*k)**2
    pbb = add_thickness_lss_effect(k, pbb0, lss)
    return pbb

def get_pbb():
    ## load data
    planck = pd.read_csv(RESULTSDIR+'data_products/pkee.dat')
    planckk = planck['k'].to_numpy()
    pkee = planck['pkee'].to_numpy()
    pkee_u = planck['pkee_u'].to_numpy()
    pkee_l = planck['pkee_l'].to_numpy()

    ## transform pkee -> pkbb (w/ regularization)
    pk_z1100 = calc_pbb(planckk, pkee, dklss)
    pk_z1100_u = calc_pbb(planckk, pkee_u, dklss)
    pk_z1100_l = calc_pbb(planckk, pkee_l, dklss)

    ## save results
    results = np.array([planckk, pk_z1100, pk_z1100_l, pk_z1100_u]).T
    table = pd.DataFrame(results, columns=['k', 'pbz1100', 'pbz1100_u', 'pbz1100_l'])
    filepath = RESULTSDIR+'data_products/pb_z1100.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))
