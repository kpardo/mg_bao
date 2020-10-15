'''
Defines paths for all the things.
'''

import os
from mg_bao import __path__

DATADIR = os.path.join(os.path.dirname(__path__._path[0]), 'data/')
RESULTSDIR = os.path.join(os.path.dirname(__path__._path[0]), 'results/')

alldirdict = {key:value for (key,value) in locals().items() if 'DIR' in key}

for i in alldirdict.values():
    if not os.path.exists(i):
        print("Creating needed directory:", i)
        os.makedirs(i)
