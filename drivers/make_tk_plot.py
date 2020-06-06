
'''
makes tk_plot
'''


import numpy as np
import pandas as pd

from mg_bao.get_ee import make_pkee
from mg_bao.camb_pk import get_camb_spectra
from mg_bao.find_pbb import get_pbb
from mg_bao.calc_mg_form import make_tk
from mg_bao.plotting import pk_plot
from mg_bao.constants import boss_h, RERUN_ANALYSIS
from mg_bao.plotting import tk_plot



def rerun_analysis():
    print('re-doing analysis')
    make_pkee()
    get_camb_spectra()
    get_pbb()
    make_tk()

def make_plot():
    ## first load data
    print('making plot')
    camb = pd.read_csv('../results/data_products/camb_pk.dat')
    cambk = camb['k'].to_numpy()
    cambpkz0 = camb['pkz0'].to_numpy()
    cambpkz1100 = camb['pkz1100'].to_numpy()
    cambpkdiv = cambpkz0/cambpkz1100
    tkdat = pd.read_csv('../results/data_products/transfer.dat')
    k = tkdat['k']
    tk = tkdat['Tk']
    tk_l = tkdat['Tk_l']
    tk_u = tkdat['Tk_u']

    ## make figure
    filepath = '/Users/kpardo/Dropbox/Apps/Overleaf/bao/transfer.png'
    tk_plot(k, tk, tk_l, tk_u, cambk, cambpkdiv, filepath)


def main():
    if RERUN_ANALYSIS:
        rerun_analysis()

    make_plot()


if __name__ == '__main__':
    main()
