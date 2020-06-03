'''
makes some figures to check pkee
this is based on old code and I have not really cleaned it up.
'''
import matplotlib.pyplot as plt
import seaborn as sns
from mg_bao.plotting import paper_plot
import pandas as pd
paper_plot()
## check by making figures
# plt.figure()
# plt.errorbar(planckl[::5], planckcl[::5], yerr=[lclerr[::5], hclerr[::5]],fmt='o', c=cs[0], label=r'Planck Cls')
# plt.errorbar(actl, actcl, yerr=acterr, c=cs[1],fmt='o',label=r'ACT Cls')
# plt.errorbar(lbinmid, cls, yerr=[cls-cl_down, cl_upp-cls],fmt='o', c='black', label=r'Combined', zorder=3)
# plt.yscale('log')
# plt.legend()
# plt.xlabel(r'$l$')
# plt.ylabel(r'$C_l$')
# plt.tight_layout()
# plt.savefig('../results/planck_cls.png')
t = pd.read_csv('../results/data_products/pkee.dat')
yerr = [t['pkee']-t['pkee_l'], t['pkee_u'] - t['pkee']]
plt.figure()
plt.errorbar(t['k'], t['pkee'], yerr=yerr, fmt='o', color='black')
plt.yscale('log')
plt.xlabel(r'$k~[\rm{Mpc}^{-1}]$')
plt.ylabel(r'$P_{EE}(k)~[\rm{Mpc}^3]$');
plt.savefig('../results/planck_pkee.png')
