## import packages
import numpy as np
import astropy.units as u
import astropy.constants as const
import pandas as pd
from scipy.stats import binned_statistic
from datetime import datetime

from mg_bao.constants import eta_star
from mg_bao.convenience import prange
from mg_bao.paths import DATADIR, RESULTSDIR

def make_pkee():
    ## load data
    planckl, planckdl, lerr, herr = np.loadtxt(DATADIR+'planck_ps/COM_PowerSpect_CMB-EE-full_R3.01.txt', unpack=True)
    actl, actcl, acterr = np.loadtxt(DATADIR+'act_ps/cl_cmb_ee.dat', unpack=True)

    ## change planck data to cl from dl
    planckcl = planckdl/(planckl*(planckl+1))*2*np.pi
    lclerr = lerr/(planckl*(planckl+1))*2*np.pi
    hclerr = herr/(planckl*(planckl+1))*2*np.pi

    ## mask data
    planckl1, planckcl1 = mask_data(planckl, planckcl, lower=24, upper=2026)
    __, lclerr1 = mask_data(planckl, lclerr, lower=24, upper=2026)
    __, hclerr1 = mask_data(planckl, hclerr, lower=24, upper=2026)
    actl1, actcl1 = mask_data(actl, actcl, upper=2001.)
    __, acterr1 = mask_data(actl, acterr, upper=2001.)

    ## bin planck data so it has the same mid bin values as act data
    bins = np.arange(25,np.max(actl1)+25, 50)
    clbin, lbinedges, __ = binned_statistic(planckl1, planckcl1, bins=bins)
    clupp, __, __ = binned_statistic(planckl1, planckcl1+hclerr1, bins=bins)
    cldown, __, __ = binned_statistic(planckl1, planckcl1-lclerr1, bins=bins)
    lbinmid = np.array([(a + b) /2 for a,b in zip(lbinedges[:-1], lbinedges[1:])])

    ## average data  where they overlap.
    newlmask = lbinmid >= np.min(actl1)
    cls = clbin[~newlmask]
    cls = np.append(cls, (clbin[newlmask] + actcl1[:-1])/2.)
    cl_upp = clupp[~newlmask]
    cl_upp = np.append(cl_upp, (clupp[newlmask]+actcl1[:-1]+acterr1[:-1])/2.)
    cl_down = cldown[~newlmask]
    cl_down = np.append(cl_down, (cldown[newlmask]+ actcl1[:-1] - acterr1[:-1])/2.)

    ## convert to pk and ks
    planckk = (lbinmid+0.5)/(eta_star)
    pkee = cls*lbinmid**2/planckk**3*np.pi
    pkee_u = cl_upp*lbinmid**2/planckk**3*np.pi
    pkee_l = cl_down*lbinmid**2/planckk**3*np.pi

    ## print to file
    results = np.array([planckk, pkee, pkee_u, pkee_l]).T
    table = pd.DataFrame(results, columns=['k', 'pkee', 'pkee_u', 'pkee_l'])
    filepath = RESULTSDIR+'data_products/pkee.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))

def mask_data(ls, data, lower=0, upper=None):
    lmask = (lower < ls) & (ls < upper)
    return ls[lmask], data[lmask]
