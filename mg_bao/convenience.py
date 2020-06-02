'''
convenient functions

prange
pshape
psum

'''

import numpy as np

prange = lambda _: print('min:', np.min(_), 'max:', np.max(_))
pshape = lambda _: print(np.shape(_))
psum = lambda _: print(np.sum(_))
