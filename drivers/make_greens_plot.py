'''
makes greens plot
'''


import numpy as np
import pandas as pd

from mg_bao.get_ee import make_pkee
from mg_bao.camb_pk import get_camb_spectra
from mg_bao.find_pbb import get_pbb
from mg_bao.calc_mg_form import make_greens, make_tk
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
    make_greens()

def make_plot():
    ## first load data
    print('making plot')
    greens= pd.read_csv(RESULTSDIR+'data_products/greens.dat')
    r = greens['r']
    Gr = greens['Gr']
    Gr_l = greens['Gr_l']
    Gr_u = greens['Gr_u']

    ## make figure
    filepath = RESULTSDIR+'greens.png'
    greens_plot(r, Gr, Gr_l, Gr_u, filepath)


def main():
    if RERUN_ANALYSIS:
        rerun_analysis()

    make_greens()
    make_plot()


if __name__ == '__main__':
    main()
