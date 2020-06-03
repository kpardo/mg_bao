'''
makes pk_plot
'''

import numpy as np
import pandas as pd

from mg_bao.get_ee import make_pkee
from mg_bao.camb_pk import get_camb_spectra
from mg_bao.find_pbb import get_pbb
from mg_bao.plotting import pk_plot
from mg_bao.constants import boss_h, RERUN_ANALYSIS


def rerun_analysis():
    print('re-doing analysis')
    make_pkee()
    get_camb_spectra()
    get_pbb()

def make_plot():
    ## first load data
    print('making plot')
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

    pk1100_lerr = pk_z1100 - pk_z1100_l
    pk1100_uerr = pk_z1100_u - pk_z1100
    pk1100yerr = [pk1100_lerr, pk1100_uerr]

    ## make figure
    filepath='/Users/kpardo/Dropbox/Apps/Overleaf/bao/power_spectra.png'
    pk_plot(planckk, pk_z1100, pk1100yerr, sdssk, sdsspk, sdsspk_err, cambk,
            cambpkz0, cambpkz1100, figpath=filepath)



def main():
    if RERUN_ANALYSIS:
        rerun_analysis()

    make_plot()


if __name__ == '__main__':
    main()
