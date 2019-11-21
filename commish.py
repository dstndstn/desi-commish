from glob import glob
import os
import numpy as np

def gfa_filename(expnum):
    pat = '/global/project/projectdirs/desi/spectro/data/*/%08i/gfa-%08i.fits.fz' % (expnum, expnum)
    fns = glob(pat)
    if len(fns) != 1:
        raise RuntimeError('Expected to find one file matching pattern', pat, '; got', len(fns))
    return fns[0]

def req_filename(expnum):
    pat = '/global/project/projectdirs/desi/spectro/data/*/%08i/request-%08i.json' % (expnum, expnum)
    fns = glob(pat)
    if len(fns) != 1:
        raise RuntimeError('Expected to find one file matching pattern', pat, '; got', len(fns))
    return fns[0]

def sub_guide_image(g):
    gclean = np.zeros((1032, 2048), np.float32)
    sub = g[:516, 50:1074]
    gclean [:516, :1024] = sub - np.median(sub)
    sub = g[516:, 50:1074]
    gclean [516:, :1024] = sub - np.median(sub)
    sub = g[:516, 1174:2198]
    gclean [:516, 1024:] = sub - np.median(sub)
    sub = g[516:, 1174:2198]
    gclean [516:, 1024:] = sub - np.median(sub)
    return gclean, 50, 50

def read_guide_image(expnum, ext='GUIDE0'):
    fn = gfa_filename(expnum)
    F = fitsio.FITS(fn)
    g = F[ext].read()
    gclean,x0,y0 = sub_guide_image(g)
    return gclean, x0, y0
