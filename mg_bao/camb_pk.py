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

## define functions
def run_camb():
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.32, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0,
            tau=0.06) ## FIXME Need to make consistent with Planck 2018.
    pars.InitPower.set_params(As=2e-9, ns=0.965, r=0) ## FIXME same.
    pars.WantTransfer = True
    pars.set_for_lmax(2500, lens_potential_accuracy=0);
    pars.set_matter_power(redshifts=[0.38, 1100], kmax=10)
    pars.set_accuracy(AccuracyBoost=7.) ## need this to get good resolution
    results= camb.get_results(pars)
    kh,zs,PK = results.get_linear_matter_power_spectrum('delta_baryon', 'delta_baryon',hubble_units=False)
    #transfer = results.get_matter_transfer_data()
    #kh_tr = transfer.transfer_data[0,:,0]
    #camb_tk = transfer.transfer_data[model.Transfer_b-1,:,1]
    print(zs)
    return kh*h, PK[0], PK[1]

def get_camb_spectra():
    ##first, get CAMB spectra
    cambk, cambpkz0, cambpkz1100 = run_camb()
    results = np.array([cambk, cambpkz0, cambpkz1100]).T
    table = pd.DataFrame(results, columns=['k', 'pkz0', 'pkz1100'])
    filepath = '../results/data_products/camb_pk.dat'
    table.to_csv(filepath, index=False)
    print('{}: made {}'.format(datetime.now().isoformat(), filepath))
