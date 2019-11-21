from glob import glob
import os
import numpy as np
import fitsio

import commish_paths

def write_astrometry_for_gfa_expnum(expnum, out_dir, **kwargs):
    #gnums=[0,2,3,5,7,8]):
    return write_astrometry_for_gfa_file(gfa_filename(expnum), out_dir, **kwargs) #gnums=gnums)

def write_astrometry_for_gfa_file(fn, out_dir, gnums=[0,2,3,5,7,8], left=False, config_fn=None):
    hdr = fitsio.read_header(fn, ext=1)
    skyra = hdr['SKYRA']
    skydec = hdr['SKYDEC']
    expnum = hdr['EXPID']
    imgfns = []

    crpix_arg = '--crpix-center'

    for gnum in gnums:
        gname = 'GUIDE%i' % gnum
        g,x0,y0 = read_guide_image_file(fn, ext=gname)
        good = read_good_pix_maps(gname)
        imgfn = os.path.join(out_dir, 'gfa-%i-%s.fits' % (expnum, gname))
        g = g * good

        gh,gw = g.shape

        if left:
            g = g[:, :gw//2]
            crpix_arg = '--crpix-x %.1f --crpix-y %.1f' % (gw+0.5, gh/2+0.5)
        fitsio.write(imgfn, g, clobber=True)
        imgfns.append(imgfn)

    if out_dir == '':
        out_dir = '.'

    if config_fn is None:
        config_fn = commish_paths.an_config_filename

    # NOTE -- definitely want --tweak-order 1, otherwise (with
    # --tweak-order 0 or --no-tweak) we get a SQUARE CD matrix!!
    cmd = ''

    an_dir = commish_paths.an_path
    if an_dir is not None:
        cmd = ('PATH=%s/bin:${PATH} ' % an_dir) + cmd

    an_py_path = commish_paths.an_py_path
    if an_py_path is not None:
        cmd = ('PYTHONPATH=%s:${PYTHONPATH} ' % an_py_path) + cmd
        
    cmd += (('solve-field --config %s --xscale 1.1' +
             ' --ra %f --dec %f --radius 2 --scale-low 0.18 --scale-high 0.24 --scale-units app --downsample 2' +
             ' --continue --tweak-order 1 --plot-scale 0.5 --objs 100' +
             ' --batch --dir %s %s ') %
            (config_fn, skyra, skydec, out_dir, crpix_arg))
    cmd += ' '.join(imgfns)
    print(cmd)
    os.system(cmd)

    wcsfns = [fn.replace('.fits', '.wcs') for fn in imgfns]
    wcsfns = [fn if os.path.exists(fn) else None for fn in wcsfns]
    return wcsfns, hdr

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
    dirnm = commish_paths.gfa_good_pix_map_dir
    fn = os.path.join(dirnm, 'desi-gfa-goodpix-%s.fits' % guide)
    return fitsio.read(fn)

def gfa_filename(expnum):
    gfa_dir = commish_paths.desi_dir
    pat = os.path.join(gfa_dir, '*', '%08i' % expnum, 'gfa-%08i.fits.fz' % expnum)
    fns = glob(pat)
    if len(fns) != 1:
        raise RuntimeError('Expected to find one file matching pattern', pat, '; got', len(fns))
    return fns[0]

def req_filename(expnum):
    dirnm = commish_paths.desi_dir
    pat = os.path.join(dirnm, '*', '%08i' % expnum, 'request-%08i.json' % expnum)
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
    return read_guide_image_file(fn, ext=ext)

def read_guide_image_file(fn, ext='GUIDE0'):
    F = fitsio.FITS(fn)
    g = F[ext].read()
    # Pull out the active silicon area
    return sub_guide_image(g)

def read_guide_image(expnum, ext='GUIDE0'):
    fn = gfa_filename(expnum)
    F = fitsio.FITS(fn)
    g = F[ext].read()
    return sub_guide_image(g)

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


