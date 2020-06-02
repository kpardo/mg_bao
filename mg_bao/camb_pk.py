## import packages
import numpy as np
import sys, platform, os
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck15 as cosmo
import camb

from mg_bao.convenience import *

## define functions
def get_camb_spectra():
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06) ## CHECKKKKKKKK
    pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
    pars.WantTransfer = True
    pars.set_for_lmax(2500, lens_potential_accuracy=0);
    pars.set_matter_power(redshifts=[0.38, 800, 1100], kmax=10)
    pars.set_accuracy(AccuracyBoost=7.) ## need this to get good resolution
    results= camb.get_results(pars)
    kh,zs,PK = results.get_linear_matter_power_spectrum('delta_baryon', 'delta_baryon',hubble_units=False)
    transfer = results.get_matter_transfer_data()
    kh_tr = transfer.transfer_data[0,:,0]
    camb_tk = transfer.transfer_data[model.Transfer_b-1,:,1]
    primordial_pk = results.Params.scalar_power(kh*h)
    return kh*h, PK[0], PK[2], kh_tr*h, camb_tk, primordial_pk, PK[1]


## first, get CAMB spectra
cambk, cambpkz0, cambpkz1100, cambkt, camb_tk, primordial_pk, pkz750 = get_camb_spectra()

results = np.array([cambk, cambpkz0, cambpkz1100, cambkt, camb_tk, primordial_pk, pkz750]).T
np.savetxt('../data/camb_pk.dat', results)
