'''
makes greens plot
'''


import numpy as np
import pandas as pd

from mg_bao.get_ee import make_pkee
from mg_bao.camb_pk import get_camb_spectra
from mg_bao.find_pbb import get_pbb
from mg_bao.calc_mg_form import make_greens
from mg_bao.plotting import pk_plot
from mg_bao.constants import boss_h, RERUN_ANALYSIS
from mg_bao.plotting import greens_plot
from mg_bao.paths import RESULTSDIR



def rerun_analysis():
    print('re-doing analysis')
    make_pkee()
    get_camb_spectra()
    get_pbb()
    make_tk()
    make_greens(ext='zeros')
    make_greens(ext='const')

def make_plot():
    ## first load data
    print('making plot')
    greens= pd.read_csv(RESULTSDIR+'data_products/greens_zeros.dat')
    r = greens['r']
    Gr = greens['Gr']
    Gr_l = greens['Gr_l']
    Gr_u = greens['Gr_u']
    greens2 = pd.read_csv(RESULTSDIR+'data_products/greens_const.dat')
    r2 = greens2['r']
    Gr2 = greens2['Gr']
    Gr2_l = greens2['Gr_l']
    Gr2_u = greens2['Gr_u']

    cambgreens = pd.read_csv(RESULTSDIR+'data_products/cambgreens.dat')
    cambr = cambgreens['r']
    cambgr = cambgreens['Gr']
    ## make figure
    filepath = RESULTSDIR+'greens.png'
    greens_plot(r, Gr, Gr_l, Gr_u, r2, Gr2, Gr2_l, Gr2_u, cambr, cambgr, filepath)


def main():
    if RERUN_ANALYSIS:
        rerun_analysis()

    make_plot()


if __name__ == '__main__':
    main()
