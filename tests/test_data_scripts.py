'''
test all scripts up to calc_mg_form
'''

from mg_bao.get_ee import *
from mg_bao.camb_pk import *
from mg_bao.find_pbb import *
from mg_bao.calc_mg_form import *

print('First get the EE power spectrum')
make_pkee()
print('Now get the CAMB spectra...this might take a minute')
#get_camb_spectra()
#get_camb_dmonly_spectra()
print('Now make the baryon PS at z=1100')
get_pbb()
print('Now get the transfer function')
make_tk()
print('Make the greens function')
make_greens(ext='zeros')
make_greens(ext='const')
