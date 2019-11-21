#! /usr/bin/env python3

'''
Assuming you're at NERSC or mayall-idl, this script can use
Astrometry.net to compute astrometric solutions for the DESI GFA guide
chips.
'''

import fitsio
import numpy as np
from glob import glob
import os

import commish_paths

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
    
