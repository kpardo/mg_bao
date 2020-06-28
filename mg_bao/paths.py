'''
Defines paths for all the things.
'''

import os
from mg_bao import __path__

DATADIR = os.path.join(os.path.dirname(__path__._path[0]), 'data/')
RESULTSDIR = os.path.join(os.path.dirname(__path__._path[0]), 'results/')

