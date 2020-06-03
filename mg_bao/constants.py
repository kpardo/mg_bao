'''
constants used throughout project
'''

import numpy as np
from astropy.cosmology import FlatLambdaCDM

RERUN_ANALYSIS = False

## set cosmology to Planck 2018 Paper I Table 6
cosmo = FlatLambdaCDM(H0=67.32, Om0=0.3158, Ob0=0.03324)

boss_h = 0.676 ## h that BOSS uses. 
h = 0.6732 ## planck 2018 h
eta_star = cosmo.comoving_distance(1059.94).value ## z_drag from Planck 2018 cosmology paper Table 2, all Planck alone
rs = 147.09 ## try rs=r_drag from Planck 2018 same table as z_drag above
lstar = np.pi*eta_star/rs
dklss = np.pi/19. ##width of last scattering -- see  Bo & David's paper.
