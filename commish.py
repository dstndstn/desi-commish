from glob import glob
import os
import numpy as np
import fitsio

def average_radec(ras, decs):
    xyz = np.zeros(3)
    for ra,dec in zip(ras, decs):
        xyz[0] += np.cos(np.deg2rad(ra)) * np.cos(np.deg2rad(dec))
        xyz[1] += np.sin(np.deg2rad(ra)) * np.cos(np.deg2rad(dec))
        xyz[2] += np.sin(np.deg2rad(dec))
    xyz /= np.sqrt(np.sum(xyz**2))
    dec = np.rad2deg(np.arcsin(xyz[2]))
    a = np.arctan2(xyz[1], xyz[0])
    a += (a<0) * (2.*np.pi)
    ra = np.rad2deg(a)
    return ra,dec

def read_good_pix_maps(guide='GUIDE0'):
    fn = '/global/cscratch1/sd/dstn/gfa-wcs/desi-gfa-goodpix-%s.fits' % guide
    return fitsio.read(fn)

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

def find_hot_pix(imgs):
    npeaks = np.zeros(imgs[0].shape)
    for g in imgs:
        peak = np.ones(g.shape, bool)
        peak[1:,:] *= (g[1:,:] > g[:-1,:])
        #print('After 1:', np.sum(peak), 'peaks')
        peak[:-1,:] *= (g[:-1,:] > g[1:,:])
        #print('After 2:', np.sum(peak), 'peaks')
        peak[:,1:] *= (g[:,1:] > g[:,:-1])
        #print('After 3:', np.sum(peak), 'peaks')
        peak[:,:-1] *= (g[:,:-1] > g[:,1:])
        #print('After 4:', np.sum(peak), 'peaks')
        peak[1:,1:] *= (g[1:,1:] > g[:-1,:-1])
        #print('After 5:', np.sum(peak), 'peaks')
        peak[1:,:-1] *= (g[1:,:-1] > g[:-1,1:])
        #print('After 6:', np.sum(peak), 'peaks')
        peak[:-1,1:] *= (g[:-1,1:] > g[1:,:-1])
        #print('After 7:', np.sum(peak), 'peaks')
        peak[:-1,:-1] *= (g[:-1,:-1] > g[1:,1:])
        #print('After 8:', np.sum(peak), 'peaks')
        npeaks += (peak*1)
    return (npeaks > int(len(imgs)*0.9))

