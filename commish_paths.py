import os

if 'NERSC_HOST' in os.environ:
    from commish_paths_nersc import *
else:
    from commish_paths_mzls import *
