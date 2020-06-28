## import packages
import numpy as np
import sys, platform, os
import astropy.units as u
import astropy.constants as const
import pandas as pd
import camb
from datetime import datetime

from mg_bao.convenience import *
from mg_bao.constants import *
from mg_bao.paths import RESULTSDIR

## define functions
def run_camb():
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.36, ombh2=0.02237, omch2=0.1200, mnu=0.06, omk=0,
            tau=0.0544) ## set to Planck 2018, cosmo paper, table 2, all Planck
    pars.InitPower.set_params(As=2.100e-9, ns=0.9649, r=0)
    pars.WantTransfer = True
    pars.set_for_lmax(2500, lens_potential_accuracy=0);
    pars.set_matter_power(redshifts=[0.38, 1100], kmax=10)
    pars.set_accuracy(AccuracyBoost=7.) ## need this to get good resolution
    results= camb.get_results(pars)
    kh,zs,PK = results.get_linear_matter_power_spectrum('delta_baryon', 'delta_baryon',hubble_units=False)
    print(zs)
    return kh*h, PK[0], PK[1]

def run_dmonly_camb():
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.36, ombh2=0., omch2=0.3153*.6736**2, mnu=0.06, omk=0,
            tau=0.0544) ## set omch2 = Omega_m*h**2 
    pars.InitPower.set_params(As=2.100e-9, ns=0.9649, r=0)
    pars.WantTransfer = True
    pars.set_for_lmax(2500, lens_potential_accuracy=0);
    pars.set_matter_power(redshifts=[0.38, 1100], kmax=10)
    pars.set_accuracy(AccuracyBoost=7.) ## need this to get good resolution
    results= camb.get_results(pars)
    kh,zs,PK = results.get_linear_matter_power_spectrum(hubble_units=False)
    print(zs)
    return kh*h, PK[0], PK[1]

def get_camb_spectra():
    ##first, get CAMB spectra
    cambk, cambpkz0, cambpkz1100 = run_camb()
    results = np.array([cambk, cambpkz0, cambpkz1100]).T
    table = pd.DataFrame(results, columns=['k', 'pkz0', 'pkz1100'])
    filepath = RESULTSDIR+'/data_products/camb_pk.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))


def get_camb_dmonly_spectra():
    ##first, get CAMB spectra
    cambk, cambpkz0, cambpkz1100 = run_dmonly_camb()
    results = np.array([cambk, cambpkz0, cambpkz1100]).T
    table = pd.DataFrame(results, columns=['k', 'pkz0', 'pkz1100'])
    filepath = RESULTSDIR+'data_products/cambdm_pk.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))
