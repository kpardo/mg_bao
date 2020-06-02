'''
plotting functions
'''


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from datetime import datetime

def paper_plot():
    sns.set_context("paper")
    sns.set_style('ticks')
    sns.set_palette('colorblind')
    figparams = {
            'text.latex.preamble': [r'\usepackage{amsmath}'],
            'text.usetex':True,
            'axes.labelsize':20.,
            'xtick.labelsize':16,
            'ytick.labelsize':16,
            'figure.figsize':[14., 12.],
            'font.family':'DejaVu Sans',
            'legend.fontsize':12}
    plt.rcParams.update(figparams)
    cs = plt.rcParams['axes.prop_cycle'].by_key()['color']

def savefig(fig, figpath, writepdf=True, dpi=450):
    ## stolen from luke
    fig.savefig(figpath, dpi=dpi, bbox_inches='tight')
    print('{}: made {}'.format(datetime.now().isoformat(), figpath))

    if writepdf:
        pdffigpath = figpath.replace('.png','.pdf')
        fig.savefig(pdffigpath, bbox_inches='tight', rasterized=True, dpi=dpi)
        print('{}: made {}'.format(datetime.now().isoformat(), pdffigpath))

    plt.close('all')

def pk_plot():
    '''
    makes power spectrum plot
    '''

