#!/usr/bin/env python3

import fitsio
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np
from glob import glob
from collections import Counter
from astrometry.util.util import Sip, Tan
import os
import json
from astrometry.util.fits import fits_table
from astrometry.util.starutil_numpy import hmsstring2ra

from commish import *

#desi_dir = '/global/project/projectdirs/desi/spectro/data/'
#J = json.load(open(desi_dir + '20191113/00027552/request-00027552.json'))

def dither_sequence(start, ending):
    dither = fits_table()
    dither.expnum = []
    dither.gfa_ra = []
    dither.gfa_dec = []
    dither.gfa05_ra = []
    dither.gfa05_dec = []
    dither.gfa27_ra = []
    dither.gfa27_dec = []
    dither.gfa38_ra = []
    dither.gfa38_dec = []
    dither.gfa_all_ra = []
    dither.gfa_all_dec = []
    dither.sky_ra = []
    dither.sky_dec = []
    headers = []

    for expnum in range(start, ending):
        fns = glob('/global/cscratch1/sd/dstn/gfa-wcs/gfa-%i-GUIDE?.wcs' % expnum)
        if len(fns) != 6:
            #print('Do not have', expnum)
            continue
        #print('Got', expnum)
        fns.sort()
        rr,dd = [],[]
        for fn in fns:
            print('Reading', fn)
            wcs = Tan(fn)
            rr.append(wcs.crval[0])
            dd.append(wcs.crval[1])
        fn = gfa_filename(expnum)
        hdr = fitsio.read_header(fn, ext=1)
        dither.expnum.append(expnum)
        headers.append(hdr)
        r,d = average_radec(rr, dd)
        dither.gfa_ra.append(r)
        dither.gfa_dec.append(d)
        r,d = average_radec([rr[0],rr[3]], [dd[0],dd[3]])
        dither.gfa05_ra.append(r)
        dither.gfa05_dec.append(d)
        r,d = average_radec([rr[1],rr[4]], [dd[1],dd[4]])
        dither.gfa27_ra.append(r)
        dither.gfa27_dec.append(d)
        r,d = average_radec([rr[2],rr[5]], [dd[2],dd[5]])
        dither.gfa38_ra.append(r)
        dither.gfa38_dec.append(d)
        dither.gfa_all_ra.append(rr)
        dither.gfa_all_dec.append(dd)
        dither.sky_ra.append(hdr['SKYRA'])
        dither.sky_dec.append(hdr['SKYDEC'])
    dither.to_np_arrays()
    dither.headers = headers
    return dither

def dither_plots(start, ending, name):
    dither = dither_sequence(start, ending)

    hdrs = dither.headers
    dither.delete_column('headers')
    #for h in hdrs:
    #    try:
    #        h['MJD-OBS']
    #    except:
    #        print('No MJD:', h)
    #mjdvals = [h['MJD-OBS'] for h in hdrs]
    #print('mjds:', mjdvals)
    dither.mjd = np.array([h.get('MJD-OBS',0) for h in hdrs])

    dither.program = np.array([h.get('PROGRAM','') for h in hdrs])

    #REQADC  = '24.13,27.29'        / [deg] requested ADC angles
    dither.reqadc1 = np.array([h.get('REQADC',(np.nan,np.nan))[0] for h in hdrs])
    dither.reqadc2 = np.array([h.get('REQADC',(np.nan,np.nan))[1] for h in hdrs])

    dither.adc1phi = np.array([h.get('ADC1PHI', np.nan) for h in hdrs])
    dither.adc2phi = np.array([h.get('ADC2PHI', np.nan) for h in hdrs])

    dither.adc1home = np.array([h.get('ADC1HOME', np.nan) for h in hdrs])
    dither.adc2home = np.array([h.get('ADC2HOME', np.nan) for h in hdrs])

    dither.adc1nrev = np.array([h.get('ADC1NREV', np.nan) for h in hdrs])
    dither.adc2nrev = np.array([h.get('ADC2NREV', np.nan) for h in hdrs])

    dither.adc1stat = np.array([h.get('ADC1STAT', np.nan) for h in hdrs])
    dither.adc2stat = np.array([h.get('ADC2STAT', np.nan) for h in hdrs])
    
    dither.writeto('%s.fits' % name)

    plt.figure(figsize=(8,8))
    plt.plot(dither.sky_ra, dither.sky_dec, 'rx-', label='SKY_RA, SKY_DEC')
    plt.plot(dither.gfa05_ra, dither.gfa05_dec, 'o-', label='Average of GUIDE0 + GUIDE5')
    plt.plot(dither.gfa27_ra, dither.gfa27_dec, 'o-', label='Average of GUIDE2 + GUIDE7')
    plt.plot(dither.gfa38_ra, dither.gfa38_dec, 'o-', label='Average of GUIDE3 + GUIDE8')
    sra,sdec = average_radec(dither.sky_ra, dither.sky_dec)
    plt.legend()
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.title('%s: SKY and GFA positions, near (%.1f, %.1f)' % (name, sra, sdec))
    plt.axis('equal');
    plt.savefig('%s-all.png' % name)
    
    plt.figure(figsize=(8,8))
    plt.plot(dither.gfa_ra, dither.gfa_dec, 'bo-', label='Mean GFA position (6 GUIDE chips)')
    plt.plot(dither.sky_ra, dither.sky_dec, 'rx-', label='SKY_RA,SKY_DEC header')
    plt.legend()
    plt.axis('equal')
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.title('Dither sequence %s' % name)
    plt.savefig('%s-avg.png' % name)
    
    plt.figure(figsize=(8,8))
    cosdec = np.cos(np.deg2rad(dither.sky_dec))
    dR = (dither.gfa_ra  - dither.sky_ra)*cosdec*3600.
    dD = (dither.gfa_dec - dither.sky_dec)*3600.
    plt.plot(dR, dD, 'ko-', label='GFA - SKY')
    for dr,dd,d in zip([dR[0],dR[-1]], [dD[0],dD[-1]], [dither[0], dither[-1]]):
        plt.text(dr, dd, '%i'%d.expnum, color='r', fontsize=16)
    plt.legend()
    plt.axis('equal')
    plt.title('Dither sequence %s: Mean GFA position MINUS SKY_RA,SKY_DEC' % name)
    plt.xlabel('delta-RA (arcsec)')
    plt.ylabel('delta-Dec (arcsec)')
    plt.savefig('%s-diff.png' % name)    
    
    return dither




def main():
    #dither_plots(27690, 27775+1, 'dither1');
    #dither_plots(27690, 27732+1, 'dither1a')
    #dither_plots(27733, 27775+1, 'dither1b')
    #dither_plots(27781, 27822+1, 'dither2');
    #dither_plots(27781, 27801+1, 'dither2a');
    #dither_plots(27802, 27822+1, 'dither2b');
    #dither_plots(27825, 27866+1, 'dither3');
    #dither_plots(27825, 27845+1, 'dither3a');
    #dither_plots(27846, 27866+1, 'dither3b');
    #dither_plots(27869, 27910+1, 'dither4');
    #dither_plots(27869, 27889+1, 'dither4a');
    #dither_plots(27890, 27910+1, 'dither4b');
    #dither_plots(27913, 27954+1, 'dither5');
    #dither_plots(27913, 27933+1, 'dither5a');
    #dither_plots(27934, 27954+1, 'dither5b');
    
    plt.figure(figsize=(8,8))
    seqs = ['1a','1b','2a','2b','3a','3b','4a','4b','5a','5b']
    for i,seq in enumerate(seqs):
        T = fits_table('dither%s.fits' % seq)
    
        if 'a' in seq:
            txt = 'ADC ON'
            c = 'r'
        else:
            txt = 'ADC OFF'
            c = 'b'
    
        if i%2 == 0:
            plt.clf()
    
        cosdec = np.cos(np.deg2rad(T.sky_dec))
        dR = (T.gfa_ra  - T.sky_ra)*cosdec*3600.
        dD = (T.gfa_dec - T.sky_dec)*3600.
        plt.plot(dR, dD, 'o-', color=c, label='%s (%s): Mean GFA - SKY' % (seq, txt))
        for dr,dd,d in zip([dR[0],dR[-1]], [dD[0],dD[-1]], [T[0], T[-1]]):
            plt.text(dr, dd, '%i'%d.expnum, color='k', fontsize=16)
    
        if i%2==1:
            plt.legend()
            plt.axis('equal')
            plt.title('Dither sequences with and without ADC')
            plt.xlabel('delta-RA (arcsec)')
            plt.ylabel('delta-Dec (arcsec)')
            plt.savefig('dithers-diff-%i.png' % (1+i//2))
    

'''
To produce all the WCS files:
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000278*/gfa* --guides 0 > /tmp/0.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000278*/gfa* --guides 2 > /tmp/2.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000278*/gfa* --guides 3 > /tmp/3.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000278*/gfa* --guides 5 > /tmp/5.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000278*/gfa* --guides 7 > /tmp/7.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000278*/gfa* --guides 8 > /tmp/8.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000279*/gfa* --guides 8 > /tmp/8b.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000279*/gfa* --guides 7 > /tmp/7b.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000279*/gfa* --guides 5 > /tmp/5b.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000279*/gfa* --guides 3 > /tmp/3b.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000279*/gfa* --guides 2 > /tmp/2b.log 2>&1 &
~/gfa-astrometry.py --out $CSCRATCH/gfa-wcs ~/desi/spectro/data/20191114/000279*/gfa* --guides 0 > /tmp/0b.log 2>&1 &
'''

main()
