#! /usr/bin/env python3
'''
Assuming you're at NERSC, this script can use Astrometry.net to compute astrometric solutions
for the DESI GFA guide chips.
'''

import fitsio
import numpy as np
from glob import glob
import os

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
        #config_fn = '/global/project/projectdirs/cosmo/work/users/dstn/index-5000/cfg'
        config_fn = '/data1/astrometry-index/cfg'
        
    # NOTE -- definitely want --tweak-order 1, otherwise (with --tweak-order 0 or --no-tweak) we get a SQUARE CD matrix!!
    #an_dir = '/global/homes/d/dstn/astrometry'
    cmd = (('solve-field --config %s --xscale 1.1' +
            ' --ra %f --dec %f --radius 2 --scale-low 0.18 --scale-high 0.24 --scale-units app --downsample 2' +
            ' --continue --tweak-order 1 --plot-scale 0.5 --objs 100' +
            ' --batch --dir %s %s ') %
           (config_fn, skyra, skydec, out_dir, crpix_arg))
    # --no-remove-lines --uniformize 0
    cmd += ' '.join(imgfns)
    print(cmd)
    os.system(cmd)

    wcsfns = [fn.replace('.fits', '.wcs') for fn in imgfns]
    wcsfns = [fn if os.path.exists(fn) else None for fn in wcsfns]
    return wcsfns, hdr

def gfa_filename(expnum):
    pat = '/global/project/projectdirs/desi/spectro/data/*/%08i/gfa-%08i.fits.fz' % (expnum, expnum)
    fns = glob(pat)
    if len(fns) != 1:
        raise RuntimeError('Expected to find one file matching pattern', pat, '; got', len(fns))
    return fns[0]

def read_guide_image(expnum, ext='GUIDE0'):
    fn = gfa_filename(expnum)
    return read_guide_image_file(fn, ext=ext)

def read_guide_image_file(fn, ext='GUIDE0'):
    F = fitsio.FITS(fn)
    g = F[ext].read()
    # Pull out the active silicon area
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

def read_good_pix_maps(guide='GUIDE0'):
    #fn = '/global/cscratch1/sd/dstn/gfa-wcs/desi-gfa-goodpix-%s.fits' % guide
    fn = '/data1/desi-gfa-wcs/desi-gfa-goodpix-%s.fits' % guide
    return fitsio.read(fn)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('gfa', type=str, help='GFA exposure number or file name', nargs='+')
    parser.add_argument('--guides', type=str, help='GUIDE chips to use (comma-separated list of integers, eg 2,5 means GUIDE2 and GUIDE5', default=None)
    parser.add_argument('--out', type=str, help='Output directory (default: current directory', default='.')
    parser.add_argument('--left', default=False, action='store_true', help='Use left half of image only')
    parser.add_argument('--config', default=None, type=str, help='Use given astrometry.net config file')

    args = parser.parse_args()
    kwargs = {}
    if args.guides is not None:
        gnums = [int(x) for x in args.guides.split(',')]
        kwargs.update(gnums=gnums)

    if args.left:
        kwargs.update(left=True)

    if args.config:
        kwargs.update(config_fn=args.config)
        
    for g in args.gfa:
        expnum = None
        try:
            expnum = int(g)
        except:
            pass
        if expnum is not None:
            wcsfns,hdr = write_astrometry_for_gfa_expnum(expnum, args.out, **kwargs)
        else:
            wcsfns,hdr = write_astrometry_for_gfa_file(g, args.out, **kwargs)

        print('WCS files:', wcsfns)

        nok = 0
        xyz = np.zeros(3)
        for fn in wcsfns:
            if fn is None:
                continue
            nok += 1
            wcshdr = fitsio.read_header(fn)
            ra = wcshdr['CRVAL1']
            dec = wcshdr['CRVAL2']
            print('  ', fn, ': center (%.4f, %.4f)' % (ra, dec))
            xyz[0] += np.cos(np.deg2rad(ra)) * np.cos(np.deg2rad(dec))
            xyz[1] += np.sin(np.deg2rad(ra)) * np.cos(np.deg2rad(dec))
            xyz[2] += np.sin(np.deg2rad(dec))
        xyz /= nok
        if nok < len(wcsfns):
            print('Warning:', (len(wcsfns)-nok), 'chips did not get a WCS solution!')

        dec = np.rad2deg(np.arcsin(xyz[2]))
        a = np.arctan2(xyz[1], xyz[0])
        a += (a<0) * (2.*np.pi)
        ra = np.rad2deg(a)
        print('RA,Dec midpoint from', nok, 'chips: (%.4f, %.4f)' % (ra, dec))

        skyra  = hdr['SKYRA']
        skydec = hdr['SKYDEC']
        print('SKY RA,Dec: (%.4f, %.4f)' % (skyra, skydec))

        cosdec = np.cos(np.deg2rad(skydec))
        dr = (ra - skyra) * cosdec * 3600.
        dd = (dec - skydec) * 3600.
        print('Offset (GFA - SKY): %.1f, %.1f arcsec' % (dr, dd))

        if nok < len(wcsfns):
            print('Warning:', (len(wcsfns)-nok), 'chips did not get a WCS solution!')
            print('You might want to re-run with --guide 0,5 or --guide 0,5,2,7 or some other')
            print('combination of chips.')

if __name__ == '__main__':
    main()
    
