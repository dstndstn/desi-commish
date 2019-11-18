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
#plt.rcParams['figure.figsize'] = [10, 5]

from commish import *

#expnums = [27552, 27617, 27619, 27621, 27623, 27625, 27627, 27629]
expnums = [27552]
gnums = [0, 5] #, 2, 7] #, 3, 8]

gfa_avgra = {}
gfa_avgdec = {}
headers = {}

desi_dir = '/global/project/projectdirs/desi/spectro/data/'

goodpix = {}
for gnum in gnums:
    goodpix[gnum] = read_good_pix_maps(guide='GUIDE%i'%gnum)

ii = 0
for expnum in expnums:
    F = fitsio.FITS(desi_dir + '20191113/%08i/guide-%08i.fits.fz' % (expnum, expnum))
    for gnum in gnums:
        gname = 'GUIDE%i' % gnum

        hdr = F[gname].read_header()
        print(hdr['DEVICE'])
        headers[(expnum,gnum)] = hdr

#         ########
        continue
#         ########

        gfa_stack = F[gname].read()
        good = goodpix[gnum]
        skyra = hdr['SKYRA']
        skydec = hdr['SKYDEC']

        gfasum = np.sum(gfa_stack, axis=0)
        print('Sum:', gfasum.shape)
        g,x0,y0 = sub_guide_image(gfasum)
        imgfn = '/global/cscratch1/sd/dstn/gfa-wcs/%i-sum-%s.fits' % (expnum, gname)
        fitsio.write(imgfn, g*good, clobber=True)
        cmd = (('solve-field --config /global/project/projectdirs/cosmo/work/users/dstn/index-5000/cfg --xscale 1.1' +
                ' --ra %f --dec %f --radius 2 --scale-low 0.18 --scale-high 0.24 --scale-units app --downsample 2' +
                ' --continue --crpix-center --dir /global/cscratch1/sd/dstn/gfa-wcs --tweak-order 1' +
                ' --no-remove-lines --uniformize 0 --plot-scale 0.5 --batch') % (skyra, skydec))
        cmd += ' ' + imgfn
        print(cmd, '> /tmp/%i.log 2>&1 &' % ii)
        
        continue

        #print(gfa_stack.shape)
        nims,h,w = gfa_stack.shape
        imgfns = []
        for i in range(nims):
            g,x0,y0 = sub_guide_image(gfa_stack[i,:,:])
            imgfn = '/global/cscratch1/sd/dstn/gfa-wcs/%i-%i-%s.fits' % (expnum, i, gname)
            fitsio.write(imgfn, g*good, clobber=True)
            imgfns.append(imgfn)
        #cmd = (('~/gfa-astrometry.py --guides %i --out /global/cscratch1/sd/dstn/gfa-wcs') %
        #       (gnum))
        cmd = (('solve-field --config /global/project/projectdirs/cosmo/work/users/dstn/index-5000/cfg --xscale 1.1' +
                ' --ra %f --dec %f --radius 2 --scale-low 0.18 --scale-high 0.24 --scale-units app --downsample 2' +
                ' --continue --crpix-center --dir /global/cscratch1/sd/dstn/gfa-wcs --tweak-order 1' +
                ' --no-remove-lines --uniformize 0 --plot-scale 0.5 --batch') % (skyra, skydec))
        cmd += ' ' + ' '.join(imgfns)
        print(cmd, '> /tmp/%i.log 2>&1 &' % ii)
        ii += 1
        
G = fits_table()
G.skyra = []
G.skydec = []
G.reqra = []
G.reqdec = []
G.gfara = []
G.gfadec = []
G.expnum = []
for gnum in gnums:
    G.set('gfa%ira'%gnum, [])
    G.set('gfa%idec'%gnum, [])

    G.set('gfa%ira_all'%gnum, [])
    G.set('gfa%idec_all'%gnum, [])

    G.set('gfa%ira_sum'%gnum, [])
    G.set('gfa%idec_sum'%gnum, [])
    
for expnum in expnums:

    hdr = headers[(expnum,gnum)]
    G.expnum.append(expnum)
    G.skyra.append(hdr['SKYRA'])
    G.skydec.append(hdr['SKYDEC'])
    G.reqra.append(hdr['REQRA'])
    G.reqdec.append(hdr['REQDEC'])

    for gnum in gnums:
        gname = 'GUIDE%i' % gnum

        ras,decs = [],[]
        for i in range(100):
            wcsfn = '/global/cscratch1/sd/dstn/gfa-wcs/%i-%i-%s.wcs' % (expnum, i, gname)
            if not os.path.exists(wcsfn):
                print('Does not exist:', wcsfn)
                break
            wcs = Tan(wcsfn)
            ras.append(wcs.crval[0])
            decs.append(wcs.crval[1])
            #print('CRPIX', wcs.crpix)

        G.get('gfa%ira_all'%gnum).append(ras)
        G.get('gfa%idec_all'%gnum).append(decs)
        G.get('gfa%ira'  % gnum).append(np.median(ras))
        G.get('gfa%idec' % gnum).append(np.median(decs))

        wcsfn = '/global/cscratch1/sd/dstn/gfa-wcs/%i-sum-%s.wcs' % (expnum, gname)
        if not os.path.exists(wcsfn):
            print('Does not exist:', wcsfn)
            break
        wcs = Tan(wcsfn)
        G.get('gfa%ira_sum' % gnum).append(wcs.crval[0])
        G.get('gfa%idec_sum' % gnum).append(wcs.crval[1])
        
G.to_np_arrays()

G.gfara = np.zeros(len(G))
G.gfadec = np.zeros(len(G))
for i,g in enumerate(G):
    r,d = average_radec([g.gfa0ra, g.gfa5ra], [g.gfa0dec, g.gfa5dec])
    G.gfara[i] = r
    G.gfadec[i] = d

G.writeto('guide.fits')



for i,g in enumerate(G):
    print('REQ RA,Dec:       (%.5f, %.5f)' % (g.reqra, g.reqdec))

    r,d = average_radec([g.gfa0ra_sum, g.gfa5ra_sum], [g.gfa0dec_sum, g.gfa5dec_sum])
    print('Summed GFA (0-5): (%.5f, %.5f)' % (r, d))

    cosdec = np.cos(np.deg2rad(g.reqdec))
    dr = (r - g.reqra ) * cosdec * 3600.
    dd = (d - g.reqdec) * 3600.
    print('Summed GFA - REQ: dRA,dDec %.1f, %.1f arcsec' % (dr,dd)) 

    dr = (g.gfara  - g.reqra ) * cosdec * 3600.
    dd = (g.gfadec - g.reqdec) * 3600.
    print('GFA - REQ: dRA,dDec %.1f, %.1f arcsec' % (dr,dd)) 
    for j,(r0,r5,d0,d5) in enumerate(zip(g.gfa0ra_all, g.gfa5ra_all, g.gfa0dec_all, g.gfa5dec_all)):
        r,d = average_radec([r0, r5], [d0, d5])
        dr = (r - g.reqra ) * cosdec * 3600.
        dd = (d - g.reqdec) * 3600.
        print('GFA(%i) - REQ: dRA,dDec %.1f, %.1f arcsec' % (j, dr,dd))



for i,g in enumerate(G):
    cosdec = np.cos(np.deg2rad(g.skydec))
    print('SKY RA,Dec:       (%.5f, %.5f)' % (g.skyra, g.skydec))

    r,d = average_radec([g.gfa0ra_sum, g.gfa5ra_sum], [g.gfa0dec_sum, g.gfa5dec_sum])
    print('Summed GFA (0-5): (%.5f, %.5f)' % (r, d))
    dr = (r - g.skyra ) * cosdec * 3600.
    dd = (d - g.skydec) * 3600.
    print('Summed GFA - SKY: dRA,dDec %.1f, %.1f arcsec' % (dr,dd)) 

    print('Medianed GFA (0-5): (%.5f, %.5f)' % (r, d))
    dr = (g.gfara  - g.skyra ) * cosdec * 3600.
    dd = (g.gfadec - g.skydec) * 3600.
    print('GFA - SKY: dRA,dDec %.1f, %.1f arcsec' % (dr,dd)) 

    for j,(r0,r5,d0,d5) in enumerate(zip(g.gfa0ra_all, g.gfa5ra_all, g.gfa0dec_all, g.gfa5dec_all)):
        dr = (r - g.skyra ) * cosdec * 3600.
        dd = (d - g.skydec) * 3600.
        print('GFA(%i) - SKY: dRA,dDec %.1f, %.1f arcsec' % (j, dr,dd))


    for j,(r0,r5,d0,d5) in enumerate(zip(g.gfa0ra_all, g.gfa5ra_all, g.gfa0dec_all, g.gfa5dec_all)):
        r,d = average_radec([r0, r5], [d0, d5])
        dr = (r - g.skyra ) * cosdec * 3600.
        dd = (d - g.skydec) * 3600.
        print('GFA(%i) - SKY: dRA,dDec %.1f, %.1f arcsec' % (j, dr,dd))


    

import sys
sys.exit(0)


#good0 = np.logical_not(find_hot_pix([g[0] for g in allg0]))
#allg5 = [read_guide_image(e, ext='GUIDE5') for e in exps]
#good5 = np.logical_not(find_hot_pix([g[0] for g in allg5]))
#for e,(g0,x0,y0) in zip(exps, allg0):
#    fitsio.write('/global/cscratch1/sd/dstn/gfa-dither/%i-g0.fits' % e, g0 * goodpix)


plt.figure(figsize=(8,8))
g0r = np.mean([wcs.crval[0] for wcs in g0tan])
g0d = np.mean([wcs.crval[1] for wcs in g0tan])
r0 = np.round(g0r, decimals=3)
d0 = np.round(g0d, decimals=3)
cosdec = np.cos(np.deg2rad(d0))
dR = [(wcs.crval[0]-r0)*cosdec*3600. for wcs in g0tan]
dD = [(wcs.crval[1]-d0)*3600. for wcs in g0tan]
plt.plot(dR, dD, 'bo-', label='GUIDE0 - (%.3f, %.3f)' % (r0,d0));
for r,d,e in zip(dR,dD, exps):
    plt.text(r, d, '%i'%e);
plt.title('GUIDE chip centers')# relative to (%.3f,%.3f)' % (r0,d0));
xl,xh = plt.xlim()
plt.xlim(xh,xl)
plt.xlabel('RA (arcsec)')
plt.ylabel('Dec (arcsec)');

g5r = np.mean([wcs.crval[0] for wcs in g5tan])
g5d = np.mean([wcs.crval[1] for wcs in g5tan])
r5 = np.round(g5r, decimals=3)
d5 = np.round(g5d, decimals=3)
dR5 = [(wcs.crval[0]-r5)*cosdec*3600. for wcs in g5tan]
dD5 = [(wcs.crval[1]-d5)*3600. for wcs in g5tan]
plt.plot(dR5, dD5, 'ro-', label='GUIDE5 - (%.3f, %.3f)' % (r5,d5));
for r,d,e in zip(dR5,dD5, exps):
    plt.text(r, d, '%i'%e);
plt.axis('equal')
plt.legend();
plt.savefig('guide-dradec.png')


# In[ ]:


plt.figure(figsize=(8,8))
dR0 = np.array(dR)
dD0 = np.array(dD)
dR5 = np.array(dR5)
dD5 = np.array(dD5)

dR = (dR0+dR5)/2.
dD = (dD0+dD5)/2.
plt.plot(dR, dD, 'ko-', label='(GUIDE0 + GUIDE5)/2 relative to (%.4f,%.4f)' % ((r0+r5)/2., (d0+d5)/2.));
xl,xh = plt.xlim()
plt.xlim(xh,xl)
plt.xlabel('RA (arcsec)')
plt.ylabel('Dec (arcsec)');
for r,d,e in zip(dR,dD, exps):
    plt.text(r, d, '%i'%e);
plt.title('Midpoint of GUIDE0 and GUIDE5')
plt.axis('equal')
plt.legend();
plt.savefig('guide-dradec-center.png')


# In[ ]:


allg2 = [read_guide_image(e, ext='GUIDE2') for e in exps]
good2 = np.logical_not(find_hot_pix([g[0] for g in allg2]))
allg7 = [read_guide_image(e, ext='GUIDE7') for e in exps]
good7 = np.logical_not(find_hot_pix([g[0] for g in allg7]))


# In[ ]:


for e,(g,x0,y0) in zip(exps, allg2):
    fitsio.write('/global/cscratch1/sd/dstn/gfa-dither/%i-g2.fits' % e, g * good2)
for e,(g,x0,y0) in zip(exps, allg7):
    fitsio.write('/global/cscratch1/sd/dstn/gfa-dither/%i-g7.fits' % e, g * good7)




g2tan = [Tan(fn) for fn in sorted(glob('/global/cscratch1/sd/dstn/gfa-dither/tanwcs/*-g2.wcs'))]
g7tan = [Tan(fn) for fn in sorted(glob('/global/cscratch1/sd/dstn/gfa-dither/tanwcs/*-g7.wcs'))]


# In[ ]:


plt.figure(figsize=(8,8))
g2r = np.mean([wcs.crval[0] for wcs in g2tan])
g2d = np.mean([wcs.crval[1] for wcs in g2tan])
r2 = np.round(g2r, decimals=3)
d2 = np.round(g2d, decimals=3)
cosdec = np.cos(np.deg2rad(d2))
dR2 = np.array([(wcs.crval[0]-r2)*cosdec*3600. for wcs in g2tan])
dD2 = np.array([(wcs.crval[1]-d2)*3600. for wcs in g2tan])
plt.plot(dR2, dD2, 'bo-', label='GUIDE2 - (%.3f, %.3f)' % (r2,d2));
for r,d,e in zip(dR2,dD2, exps):
    plt.text(r, d, '%i'%e);
plt.title('GUIDE chip centers')
xl,xh = plt.xlim()
plt.xlim(xh,xl)
plt.xlabel('RA (arcsec)')
plt.ylabel('Dec (arcsec)');

g7r = np.mean([wcs.crval[0] for wcs in g7tan])
g7d = np.mean([wcs.crval[1] for wcs in g7tan])
r7 = np.round(g7r, decimals=3)
d7 = np.round(g7d, decimals=3)
dR7 = np.array([(wcs.crval[0]-r7)*cosdec*3600. for wcs in g7tan])
dD7 = np.array([(wcs.crval[1]-d7)*3600. for wcs in g7tan])
plt.plot(dR7, dD7, 'ro-', label='GUIDE7 - (%.3f, %.3f)' % (r7,d7));
for r,d,e in zip(dR7,dD7, exps):
    plt.text(r, d, '%i'%e);
plt.axis('equal')
plt.legend();
plt.savefig('guide-dradec-27.png')


# In[ ]:


plt.figure(figsize=(8,8))
dR = (dR2+dR7)/2.
dD = (dD2+dD7)/2.
plt.plot(dR, dD, 'ko-', label='(GUIDE2 + GUIDE7)/2 relative to (%.4f,%.4f)' % ((r2+r7)/2., (d2+d7)/2.));
xl,xh = plt.xlim()
plt.xlim(xh,xl)
plt.xlabel('RA (arcsec)')
plt.ylabel('Dec (arcsec)');
for r,d,e in zip(dR,dD, exps):
    plt.text(r, d, '%i'%e);
plt.title('Midpoint of GUIDE2 and GUIDE7')
plt.axis('equal')
plt.legend();
plt.savefig('guide-dradec-center-27.png')


# In[ ]:


allg3 = [read_guide_image(e, ext='GUIDE3') for e in exps]
good3 = np.logical_not(find_hot_pix([g[0] for g in allg3]))
allg8 = [read_guide_image(e, ext='GUIDE8') for e in exps]
good8 = np.logical_not(find_hot_pix([g[0] for g in allg8]))


# In[ ]:


for e,(g,x0,y0) in zip(exps, allg3):
    fitsio.write('/global/cscratch1/sd/dstn/gfa-dither/%i-g3.fits' % e, g * good3)
for e,(g,x0,y0) in zip(exps, allg8):
    fitsio.write('/global/cscratch1/sd/dstn/gfa-dither/%i-g8.fits' % e, g * good8)


# In[ ]:


# NOTE -- definitely want --tweak-order 1, otherwise (with --tweak-order 0 or --no-tweak) we get a SQUARE CD matrix!!
# solve-field --config cfg --xscale 1.1 --ra 20.7 --dec 30.6 --radius 2 --scale-low 0.18 --scale-high 0.24 --scale-units app --downsample 2 -v --continue --crpix-center --dir $CSCRATCH/gfa-dither/tanwcs --batch --tweak-order 1 --no-remove-lines --uniformize 0 --plot-scale 0.5 $CSCRATCH/gfa-dither/*-g3.fits > g3.log 2>&1 &


# In[ ]:


g3tan = [Tan(fn) for fn in sorted(glob('/global/cscratch1/sd/dstn/gfa-dither/tanwcs/*-g3.wcs'))]
g8tan = [Tan(fn) for fn in sorted(glob('/global/cscratch1/sd/dstn/gfa-dither/tanwcs/*-g8.wcs'))]


# In[ ]:


plt.figure(figsize=(8,8))
g3r = np.mean([wcs.crval[0] for wcs in g3tan])
g3d = np.mean([wcs.crval[1] for wcs in g3tan])
r3 = np.round(g3r, decimals=3)
d3 = np.round(g3d, decimals=3)
cosdec = np.cos(np.deg2rad(d3))
dR3 = np.array([(wcs.crval[0]-r3)*cosdec*3600. for wcs in g3tan])
dD3 = np.array([(wcs.crval[1]-d3)*3600. for wcs in g3tan])
plt.plot(dR3, dD3, 'bo-', label='GUIDE3 - (%.3f, %.3f)' % (r3,d3));
for r,d,e in zip(dR3,dD3, exps):
    plt.text(r, d, '%i'%e);
plt.title('GUIDE chip centers')
xl,xh = plt.xlim()
plt.xlim(xh,xl)
plt.xlabel('RA (arcsec)')
plt.ylabel('Dec (arcsec)');

g8r = np.mean([wcs.crval[0] for wcs in g8tan])
g8d = np.mean([wcs.crval[1] for wcs in g8tan])
r8 = np.round(g8r, decimals=3)
d8 = np.round(g8d, decimals=3)
dR8 = np.array([(wcs.crval[0]-r8)*cosdec*3600. for wcs in g8tan])
dD8 = np.array([(wcs.crval[1]-d8)*3600. for wcs in g8tan])
plt.plot(dR8, dD8, 'ro-', label='GUIDE7 - (%.3f, %.3f)' % (r8,d8));
for r,d,e in zip(dR8,dD8, exps):
    plt.text(r, d, '%i'%e);
plt.axis('equal')
plt.legend();
plt.savefig('guide-dradec-38.png')


# In[ ]:


plt.figure(figsize=(8,8))
dR = (dR3+dR8)/2.
dD = (dD3+dD8)/2.
plt.plot(dR, dD, 'ko-', label='(GUIDE3 + GUIDE8)/2 relative to (%.4f,%.4f)' % ((r3+r8)/2., (d3+d8)/2.));
xl,xh = plt.xlim()
plt.xlim(xh,xl)
plt.xlabel('RA (arcsec)')
plt.ylabel('Dec (arcsec)');
for r,d,e in zip(dR,dD, exps):
    plt.text(r, d, '%i'%e);
plt.title('Midpoint of GUIDE3 and GUIDE8')
plt.axis('equal')
plt.legend();
plt.savefig('guide-dradec-center-38.png')


# In[ ]:


adc1phi = []
adc2phi = []
mid_ra = []
mid_dec = []
skyra = []
skydec = []
for e in exps:
    tans = glob('/global/cscratch1/sd/dstn/gfa-dither/tanwcs/%i-g*.wcs' % e)
    assert(len(tans) == 6)
    tans = [Tan(fn) for fn in tans]

    r = np.mean([wcs.crval[0] for wcs in tans])
    d = np.mean([wcs.crval[1] for wcs in tans])
    mid_ra.append(r)
    mid_dec.append(d)

    fn = gfa_filename(e)
    print('Reading', fn)
    hdr = fitsio.read_header(fn, ext=1)
    adc1phi.append(hdr['ADC1PHI'])
    adc2phi.append(hdr['ADC2PHI'])
    
    skyra.append(hdr['SKYRA'])
    skydec.append(hdr['SKYDEC'])

    
adc1phi = np.array(adc1phi)
adc2phi = np.array(adc2phi)
mid_ra = np.array(mid_ra)
mid_dec = np.array(mid_dec)
skyra = np.array(skyra)
skydec = np.array(skydec)


# In[ ]:


plt.plot(mid_ra, mid_dec, 'bo-');


# In[ ]:


plt.plot(exps, adc1phi - 50, 'b-', label='ADC1PHI - 50')
plt.plot(exps, adc2phi - 100, 'r-', label='ADC2PHI - 100')
plt.xlabel('Exposure Number')
plt.ylabel('relative PHI (deg)')
plt.legend()
plt.savefig('dither-adc.png')


# In[ ]:


cosdec = np.cos(np.deg2rad(np.mean(skydec)))
dr = (mid_ra - skyra)*3600.*cosdec
dd = (mid_dec - skydec)*3600
plt.plot(dr, dd, 'bo-');
plt.axis('equal')
plt.xlabel('GUIDE array midpoints - SKY_RA (arcsec)')
plt.ylabel('GUIDE array midpoints - SKY_DEC (arcsec)')
plt.text(dr[0], dd[0], '%i'%exps[0])
plt.text(dr[-1], dd[-1], '%i'%exps[-1])
plt.savefig('dither-dradec.png')


# In[ ]:


plt.figure(figsize=(8,8))
rd0 = np.vstack([wcs.crval for wcs in g0tan])
rd5 = np.vstack([wcs.crval for wcs in g5tan])
rd2 = np.vstack([wcs.crval for wcs in g2tan])
rd7 = np.vstack([wcs.crval for wcs in g7tan])
rd3 = np.vstack([wcs.crval for wcs in g3tan])
rd8 = np.vstack([wcs.crval for wcs in g8tan])

#rd = (rd0 + rd2 + rd3 + rd5 + rd7 + rd8) / 6.
#rr0 = np.mean(rd[:,0])
#dd0 = np.mean(rd[:,1])

rd = (rd0 + rd5)/2.
rr0 = np.mean(rd[:,0])
dd0 = np.mean(rd[:,1])
plt.plot(rd[:,0]-rr0, rd[:,1]-dd0, 'o-', label='(GUIDE0 + GUIDE5)/2')

rd = (rd2 + rd7)/2.
rr0 = np.mean(rd[:,0])
dd0 = np.mean(rd[:,1])
plt.plot(rd[:,0]-rr0, rd[:,1]-dd0, 'o-', label='(GUIDE2 + GUIDE7)/2')

rd = (rd3 + rd8)/2.
rr0 = np.mean(rd[:,0])
dd0 = np.mean(rd[:,1])
plt.plot(rd[:,0]-rr0, rd[:,1]-dd0, 'o-', label='(GUIDE3 + GUIDE8)/2')

xl,xh = plt.xlim()
plt.xlim(xh,xl)
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)');
#for r,d,e in zip(dR,dD, exps):
#    plt.text(r, d, '%i'%e);
plt.title('Midpoint of opposing GUIDE arrays')
plt.axis('equal')
plt.legend();
plt.savefig('guide-radec-all.png')


# In[ ]:


dR = (rd[:,0] - rr0) * 3600. * cosdec
dD = (rd[:,1] - dd0) * 3600.
dadc = (adc1phi - 50)
plt.plot(dR + dadc*0.8, dD + dadc*0.6, 'o-')
plt.xlabel('Corrected RA (arcsec)')
plt.ylabel('Corrected Dec (arcsec)')
plt.savefig('corr-adc.png')


# In[ ]:


for gnum,good in [(0, good0), (5, good5), (3, good3), (8, good8),
                  (2, good2), (7, good7)]:
    print(good.dtype)
    fitsio.write('/global/cscratch1/sd/dstn/gfa-wcs/desi-gfa-goodpix-GUIDE%i.fits' % gnum, good.astype(np.uint8), clobber=True)
    


# In[ ]:


for e in [27541, 27546, 27557, 27631, 27643, 27647, 27648, 27649, 27653, 27655, 27656, 27658, 27659,
         27660, 27661, 27662, 27663, 27664, 27665, 27666, 27667]:
    imgfns = []
    for gnum,good in [(0, good0), (5, good5), (3, good3), (8, good8),
                      (2, good2), (7, good7)]:
        gname = 'GUIDE%i' % gnum
        g,x0,y0 = read_guide_image(e, ext=gname)
        imgfn = '/global/cscratch1/sd/dstn/gfa-wcs/%i-g%i.fits' % (e, gnum)
        fitsio.write(imgfn, g * good)
        hdr = fitsio.read_header(gfa_filename(e), ext=2)
        skyra = hdr['SKYRA']
        skydec = hdr['SKYDEC']
        imgfns.append(imgfn)
    cmd = (('solve-field --config /global/project/projectdirs/cosmo/work/users/dstn/index-5000/cfg --xscale 1.1' +
            ' --ra %f --dec %f --radius 2 --scale-low 0.18 --scale-high 0.24 --scale-units app --downsample 2' +
            ' -v --continue --crpix-center --dir /global/cscratch1/sd/dstn/gfa-wcs --tweak-order 1' +
            ' --no-remove-lines --uniformize 0 --plot-scale 0.5 --batch') % (skyra, skydec))
    cmd += ' ' + ' '.join(imgfns)
    print(cmd)
        #os.system(cmd)


# In[ ]:





# In[ ]:





# In[ ]:


T = fits_table()

T.sky_ra = []
T.sky_dec = []
T.gfa_ra = []
T.gfa_dec = []
T.adc1phi = []
T.req_ha = []
T.req_dec = []
T.lst = []
T.expnum = []
T.exptime = []

for e in [27541, 27546, 27557, 27631, 27643, 27647, 27648, 27649, 27653, 27655, 27656, 27658, 27659,
         27660, 27661, 27662, 27663, 27664, 27665, 27666, 27667]:
    wcs = []
    missing = False
    for gnum,good in [(0, good0), (5, good5), (3, good3), (8, good8),
                      (2, good2), (7, good7)]:
        gname = 'GUIDE%i' % gnum
        wcsfn = '/global/cscratch1/sd/dstn/gfa-wcs/%i-g%i.wcs' % (e, gnum)
        if not os.path.exists(wcsfn):
            print('Missing:', wcsfn)
            missing = True
            break
        wcs.append(Tan(wcsfn))
    if missing:
        continue
        
    J = json.load(open(req_filename(e), 'r'))
    #print('Request:', J)
    T.req_ha.append(J.get('REQHA', np.nan))
    T.req_dec.append(J.get('REQDEC', np.nan))
    hdr = fitsio.read_header(gfa_filename(e), ext=1)
    sr = hdr['SKYRA']
    sd = hdr['SKYDEC']
    adc1 = hdr.get('ADC1PHI', np.nan)
    T.lst.append(hmsstring2ra(hdr['ST']))
    
    rd = np.vstack([w.crval for w in wcs])
    midrd = np.mean(rd, axis=0)
    print('GFA middle RA,Dec', midrd)
    
    T.sky_ra.append(sr)
    T.sky_dec.append(sd)
    T.gfa_ra.append(midrd[0])
    T.gfa_dec.append(midrd[1])
    T.adc1phi.append(adc1)
    T.expnum.append(e)
    T.exptime.append(hdr['EXPTIME'])

#sky_ra = np.array(sky_ra)
#sky_dec = np.array(sky_dec)
#gfa_ra = np.array(gfa_ra)
#gfa_dec = np.array(gfa_dec)
#adc1phi = np.array(adc1phi)
#req_ha = np.array(req_ha)
#req_dec = np.array(req_dec)
#lst = np.array(lst)
T.to_np_arrays()


# In[ ]:


T.req_ra = T.lst - T.req_ha


# In[ ]:


plt.plot(T.sky_ra, T.sky_dec, 'bx')
plt.plot(T.gfa_ra, T.gfa_dec, 'r.');
I = np.flatnonzero(np.isfinite(T.req_ra))
plt.plot(T.req_ra, T.req_dec, 'go');


# In[ ]:


adc1phi


# In[ ]:


cosdec = np.cos(np.deg2rad(sky_dec))
plt.plot(adc1phi, (gfa_ra - sky_ra)*cosdec * 3600., 'bo');
plt.xlabel('ADC1PHI (deg)')
plt.ylabel('GFA RA - SKY_RA (arcsec)')
plt.savefig('gfa-adc-ra.png');


# In[ ]:


plt.plot(adc1phi, (gfa_dec - sky_dec)*3600., 'ro')
plt.xlabel('ADC1PHI (deg)')
plt.ylabel('GFA Dec - SKY_DEC (arcsec)')
plt.savefig('gfa-adc-dec.png')


# In[ ]:


I = np.flatnonzero(np.isfinite(T.req_ra))
G = T[I]
cosdec = np.cos(np.deg2rad(G.gfa_dec))
plt.figure(figsize=(6,6))
dR = (G.gfa_ra - G.req_ra)*cosdec * 3600.
dD = (G.gfa_dec - G.req_dec)*3600
plt.plot(dR, dD, 'bo');
#plt.scatter(dR, dD, c=G.exptime);
for r,d,e in zip(dR,dD,G.expnum):
    plt.text(r, d, '%i'%e)
plt.xlabel('GFA_RA - REQ_RA (arcsec)')
plt.ylabel('GFA_DEC - REQ_DEC (arcsec)');

plt.savefig('gfa-req.png')


# In[38]:


desi_dir = '/global/project/projectdirs/desi/spectro/data/'
J = json.load(open(desi_dir + '20191113/00027552/request-00027552.json'))


# In[ ]:


req_ra = J['REQRA']
req_dec = J['REQDEC']


# In[ ]:


t = T[T.expnum == 27552]


# In[ ]:


t


# In[ ]:


F = fitsio.FITS(desi_dir + '20191113/00027552/guide-00027552.fits.fz')


# In[ ]:


len(F)


# In[ ]:


hdr = F[2].read_header()


# In[ ]:


hdr['DEVICE']


# In[ ]:


gfa_stack = F[2].read()


# In[ ]:


gfa_stack.shape


# In[ ]:


img = gfa_stack[2,:,:]
mn,mx = np.percentile(img.ravel(), [25,98])
plt.imshow(img, vmin=mn, vmax=mx);


# In[33]:




# In[ ]:


expnum


# In[147]:


expnum = 27552
#expnum = 27617
nims = 58

gfa_ras = []
gfa_decs = []
for i in range(nims):
    wcsfn = '/global/cscratch1/sd/dstn/gfa-wcs/%i-%i-%s.wcs' % (expnum, i, gname)
    wcs = Tan(wcsfn)
    gfa_ras.append(wcs.crval[0])
    gfa_decs.append(wcs.crval[1])

    #break

gfa_avgra[expnum] = np.median(gfa_ras)
gfa_avgdec[expnum] = np.median(gfa_decs)


# In[39]:


for expnum in expnums:
    F = fitsio.FITS(desi_dir + '20191113/%08i/guide-%08i.fits.fz' % (expnum, expnum))
    for gnum in [0, 5]:
        gname = 'GUIDE%i' % gnum
        hdr = F[gname].read_header()
        print(hdr['DEVICE'])
        headers[(expnum,gnum)] = hdr


# In[171]:


#for expnum in [27552, 27617]:
expnums = [27552, 27617, 27619, 27621, 27623, 27625, 27627, 27629]
for expnum in expnums:
    ras,decs = [],[]
    for gnum in [0, 5]:
        gname = 'GUIDE%i' % gnum
        wcspat = '/global/cscratch1/sd/dstn/gfa-wcs/%i-*-%s.wcs' % (expnum, gname)
        fns = glob(wcspat)
        gra,gdec = [],[]
        for fn in fns:
            wcs = Tan(fn)
            gra.append(wcs.crval[0])
            gdec.append(wcs.crval[1])
        #print('Expnum', expnum, 'GUIDE', gnum, 'RA', gra, 'Dec', gdec)
        ras.append(np.median(gra))
        decs.append(np.median(gdec))

        gfa_avgra[(expnum, gnum)] = np.median(gra)
        gfa_avgdec[(expnum, gnum)] = np.median(gdec)
        
    if len(ras):
        r,d = average_radec(ras, decs)
        gfa_avgra [expnum] = r
        gfa_avgdec[expnum] = d


# In[168]:


gfa_avgra, gfa_avgdec


# In[172]:
# for expnum in expnums:
#     F = fitsio.FITS(desi_dir + '20191113/%08i/guide-%08i.fits.fz' % (expnum, expnum))
#     for gnum in [0, 5]:
#         gname = 'GUIDE%i' % gnum
#         hdr = F[gname].read_header()
#         print(hdr['DEVICE'])
#         headers[(expnum,gnum)] = hdr
    


G = fits_table()
G.skyra = []
G.skydec = []
G.reqra = []
G.reqdec = []
G.gfara = []
G.gfadec = []
G.expnum = []

G.gfa0ra = []
G.gfa0dec = []
G.gfa5ra = []
G.gfa5dec = []

#expnums = [27552] #, 27617]
for e in expnums:
    gnum=0


# In[173]:


plt.plot(G.skyra, G.skydec, 'bx')
plt.plot(G.reqra, G.reqdec, 'ro')
plt.plot(G.gfara, G.gfadec, 'g.');
#plt.plot(G.gfa0ra, G.gfa0dec, 'm.');
#plt.plot(G.gfa5ra, G.gfa5dec, 'y.');
#plt.axis('equal')


# In[174]:


GG = G[np.isfinite(G.gfara)]
cosdec = np.cos(np.deg2rad(GG.reqdec))
dR1 = (GG.gfara - GG.reqra) * cosdec * 3600.
dD1 = (GG.gfadec - GG.reqdec) * 3600.
plt.figure(figsize=(8,8))
dR2 = (GG.skyra - GG.reqra) * cosdec * 3600.
dD2 = (GG.skydec - GG.reqdec) * 3600.
plt.plot(dR1, dD1, 'bx', label='GFA - REQ')
plt.plot(dR2, dD2, 'ro', label='SKY - REQ')
plt.plot([dR1, dR2], [dD1, dD2], 'k-')
#for r,d,e in zip((dR1+dR2)/2, (dD1+dD2)/2, expnums):
for r,d,e in zip((dR1*0.25+dR2*0.75), (dD1*0.25+dD2*0.75), expnums):
    plt.text(r,d,'%i'%e)
plt.legend();
plt.xlabel('RA difference (arcsec)')
plt.ylabel('Dec difference (arcsec)')
plt.axis('equal')
plt.savefig('gfa-sky-req.png');


# In[157]:


dD2


# In[175]:


G[np.isfinite(G.gfara)].writeto('gfa-sky-req.fits')


# In[ ]:


plt.plot(G.gfara - G.skyra, G.gfadec - G.skydec, 'bx', label='GFA - SKY')
plt.plot(G.gfara - G.reqra, G.gfadec - G.reqdec, 'ro', label='GFA - REQ')
plt.plot(np.vstack([G.gfara - G.reqra, G.gfara - G.skyra]),
         np.vstack([G.gfadec - G.reqdec, G.gfadec - G.skydec]), 'k-')
plt.legend();
plt.xlabel('RA difference, no cos(Dec) term, (deg)')
plt.ylabel('Dec difference (deg)')
plt.savefig('gfa-sky-req.png');


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# 2019-11-14 dither test 27690-27775


# In[ ]:





# In[ ]:





# In[ ]:


dither1 = fits_table()
dither1.expnum = []
dither1.gfa_ra = []
dither1.gfa_dec = []
dither1.gfa05_ra = []
dither1.gfa05_dec = []
dither1.gfa27_ra = []
dither1.gfa27_dec = []
dither1.gfa38_ra = []
dither1.gfa38_dec = []
dither1.gfa_all_ra = []
dither1.gfa_all_dec = []

dither1.sky_ra = []
dither1.sky_dec = []


wcsfns = []
for expnum in range(27690, 27775):
    fns = glob('/global/cscratch1/sd/dstn/gfa-wcs/gfa-%i-GUIDE?.wcs' % expnum)
    if len(fns) != 6:
        print('Do not have', expnum)
        continue
    fns.sort()
    #
    #fn2 = '/global/cscratch1/sd/dstn/gfa-wcs/gfa-%i-GUIDE2.wcs' % expnum
    #fn7 = '/global/cscratch1/sd/dstn/gfa-wcs/gfa-%i-GUIDE7.wcs' % expnum
    #if not (os.path.exists(fn2) and os.path.exists(fn7)):
    #    print('Do not have', expnum)
    #    continue
    rr,dd = [],[]
    #for fn in [fn2,fn7]:
    for fn in fns:
        wcs = Tan(fn)
        rr.append(wcs.crval[0])
        dd.append(wcs.crval[1])
    fn = gfa_filename(expnum)
    #print(fn)
    hdr = fitsio.read_header(fn, ext=1)
    
    dither1.expnum.append(expnum)
    dither1.gfa_ra.append(np.mean(rr))
    dither1.gfa_dec.append(np.mean(dd))

    dither1.gfa05_ra.append((rr[0] + rr[3])/2.)
    dither1.gfa05_dec.append((dd[0] + dd[3])/2.)
    dither1.gfa27_ra.append((rr[1] + rr[4])/2.)
    dither1.gfa27_dec.append((dd[1] + dd[4])/2.)
    dither1.gfa38_ra.append((rr[2] + rr[5])/2.)
    dither1.gfa38_dec.append((dd[2] + dd[5])/2.)

    dither1.gfa_all_ra.append(rr)
    dither1.gfa_all_dec.append(dd)

    dither1.sky_ra.append(hdr['SKYRA'])
    dither1.sky_dec.append(hdr['SKYDEC'])
    
dither1.to_np_arrays()


# In[ ]:


_,ng = dither1.gfa_all_ra.shape
for i in range(ng):
    plt.plot(dither1.gfa_all_ra[:,i], dither1.gfa_all_dec[:,i], 'o-', label='%i'%i)
plt.legend()


# In[ ]:


#plt.scatter(dither1.gfa_ra - dither1.sky_ra, dither1.gfa_dec - dither1.sky_dec,
#           c=dither1.expnum);

#plt.plot(dither1.gfa_ra - dither1.sky_ra, dither1.gfa_dec - dither1.sky_dec, 'bo-');

plt.plot(dither1.gfa05_ra - dither1.sky_ra, dither1.gfa05_dec - dither1.sky_dec, 'bo-');
plt.plot(dither1.gfa27_ra - dither1.sky_ra, dither1.gfa27_dec - dither1.sky_dec, 'ro-');
plt.plot(dither1.gfa38_ra - dither1.sky_ra, dither1.gfa38_dec - dither1.sky_dec, 'go-');

for d in dither1:
    plt.text(d.gfa38_ra - d.sky_ra, d.gfa38_dec - d.sky_dec, '%i'%d.expnum)


# In[ ]:


plt.figure(figsize=(8,8))
#plt.plot(dither1.gfa_ra, dither1.gfa_dec, 'bo-')
plt.plot(dither1.gfa05_ra, dither1.gfa05_dec, 'o-')
plt.plot(dither1.gfa27_ra, dither1.gfa27_dec, 'o-')
plt.plot(dither1.gfa38_ra, dither1.gfa38_dec, 'o-')
plt.plot(dither1.sky_ra, dither1.sky_dec, 'rx-')
plt.axis('equal')


# In[ ]:


plt.figure(figsize=(8,8))
plt.plot(dither1.gfa_ra, dither1.gfa_dec, 'bo-', label='Mean GFA position')
plt.plot(dither1.sky_ra, dither1.sky_dec, 'rx-', label='SKY_RA,SKY_DEC header')
plt.legend()
plt.axis('equal')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.savefig('dither1-20191114.png')


# In[ ]:


plt.subplot(2,1,1)
plt.plot(dither1.expnum, dither1.sky_ra, 'bo-')
plt.xlim(27720, 27740)
plt.subplot(2,1,2)
plt.plot(dither1.expnum, dither1.sky_dec, 'ro-')
plt.xlim(27720, 27740)


# In[ ]:


plt.figure(figsize=(8,8))
cosdec = np.cos(np.deg2rad(dither1.sky_dec))
dR = (dither1.gfa_ra - dither1.sky_ra)*cosdec*3600.
dD = (dither1.gfa_dec - dither1.sky_dec)*3600.
plt.plot(dR, dD, 'ko-', label='GFA - SKY')
for dr,dd,d in zip([dR[0],dR[-1]], [dD[0],dD[-1]], [dither1[0], dither1[-1]]):
    plt.text(dr, dd, '%i'%d.expnum, color='r', fontsize=16)
plt.legend()
plt.axis('equal')
plt.title('Mean GFA position MINUS SKY_RA,SKY_DEC')
plt.xlabel('delta-RA (arcsec)')
plt.ylabel('delta-Dec (arcsec)')
plt.savefig('dither1-20191114-diff.png')


# In[ ]:


dither2 = fits_table()
dither2.expnum = []
dither2.gfa_ra = []
dither2.gfa_dec = []
dither2.gfa05_ra = []
dither2.gfa05_dec = []
dither2.gfa27_ra = []
dither2.gfa27_dec = []
dither2.gfa38_ra = []
dither2.gfa38_dec = []
dither2.gfa_all_ra = []
dither2.gfa_all_dec = []
dither2.sky_ra = []
dither2.sky_dec = []


wcsfns = []
for expnum in range(27776, 27822):
                    #27955):
    fns = glob('/global/cscratch1/sd/dstn/gfa-wcs/gfa-%i-GUIDE?.wcs' % expnum)
    if len(fns) != 6:
        print('Do not have', expnum)
        continue
    print('Got', expnum)
    fns.sort()
    rr,dd = [],[]
    for fn in fns:
        wcs = Tan(fn)
        rr.append(wcs.crval[0])
        dd.append(wcs.crval[1])
    fn = gfa_filename(expnum)
    hdr = fitsio.read_header(fn, ext=1)
    dither2.expnum.append(expnum)
    dither2.gfa_ra.append(np.mean(rr))
    dither2.gfa_dec.append(np.mean(dd))
    dither2.gfa05_ra.append((rr[0] + rr[3])/2.)
    dither2.gfa05_dec.append((dd[0] + dd[3])/2.)
    dither2.gfa27_ra.append((rr[1] + rr[4])/2.)
    dither2.gfa27_dec.append((dd[1] + dd[4])/2.)
    dither2.gfa38_ra.append((rr[2] + rr[5])/2.)
    dither2.gfa38_dec.append((dd[2] + dd[5])/2.)
    dither2.gfa_all_ra.append(rr)
    dither2.gfa_all_dec.append(dd)
    dither2.sky_ra.append(hdr['SKYRA'])
    dither2.sky_dec.append(hdr['SKYDEC'])
dither2.to_np_arrays()


# In[ ]:


plt.figure(figsize=(8,8))
#plt.plot(dither1.gfa_ra, dither1.gfa_dec, 'bo-')
plt.plot(dither2.gfa05_ra, dither2.gfa05_dec, 'o-')
plt.plot(dither2.gfa27_ra, dither2.gfa27_dec, 'o-')
plt.plot(dither2.gfa38_ra, dither2.gfa38_dec, 'o-')
plt.plot(dither2.sky_ra, dither2.sky_dec, 'rx-')
plt.axis('equal');


# In[ ]:


plt.figure(figsize=(8,8))
plt.plot(dither2.gfa_ra, dither2.gfa_dec, 'bo-', label='Mean GFA position')
plt.plot(dither2.sky_ra, dither2.sky_dec, 'rx-', label='SKY_RA,SKY_DEC header')
plt.legend()
plt.axis('equal')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.savefig('dither2-20191114.png')


# In[ ]:


plt.figure(figsize=(8,8))
cosdec = np.cos(np.deg2rad(dither2.sky_dec))
dR = (dither2.gfa_ra  - dither2.sky_ra)*cosdec*3600.
dD = (dither2.gfa_dec - dither2.sky_dec)*3600.
plt.plot(dR, dD, 'ko-', label='GFA - SKY')
for dr,dd,d in zip([dR[0],dR[-1]], [dD[0],dD[-1]], [dither1[0], dither1[-1]]):
    plt.text(dr, dd, '%i'%d.expnum, color='r', fontsize=16)
plt.legend()
plt.axis('equal')
plt.title('Mean GFA position MINUS SKY_RA,SKY_DEC')
plt.xlabel('delta-RA (arcsec)')
plt.ylabel('delta-Dec (arcsec)')
plt.savefig('dither2-20191114-diff.png')


# In[142]:




# In[143]:


def dither_sequence(start, ending):
    dither2 = fits_table()
    dither2.expnum = []
    dither2.gfa_ra = []
    dither2.gfa_dec = []
    dither2.gfa05_ra = []
    dither2.gfa05_dec = []
    dither2.gfa27_ra = []
    dither2.gfa27_dec = []
    dither2.gfa38_ra = []
    dither2.gfa38_dec = []
    dither2.gfa_all_ra = []
    dither2.gfa_all_dec = []
    dither2.sky_ra = []
    dither2.sky_dec = []
    dither2.header = []

    for expnum in range(start, ending):
        fns = glob('/global/cscratch1/sd/dstn/gfa-wcs/gfa-%i-GUIDE?.wcs' % expnum)
        if len(fns) != 6:
            #print('Do not have', expnum)
            continue
        #print('Got', expnum)
        fns.sort()
        rr,dd = [],[]
        for fn in fns:
            wcs = Tan(fn)
            rr.append(wcs.crval[0])
            dd.append(wcs.crval[1])
        fn = gfa_filename(expnum)
        hdr = fitsio.read_header(fn, ext=1)
        dither2.expnum.append(expnum)
        dither2.header.append(hdr)
        r,d = average_radec(rr, dd)
        dither2.gfa_ra.append(r)
        dither2.gfa_dec.append(d)
        r,d = average_radec([rr[0],rr[3]], [dd[0],dd[3]])
        dither2.gfa05_ra.append(r)
        dither2.gfa05_dec.append(d)
        r,d = average_radec([rr[1],rr[4]], [dd[1],dd[4]])
        dither2.gfa27_ra.append(r)
        dither2.gfa27_dec.append(d)
        r,d = average_radec([rr[2],rr[5]], [dd[2],dd[5]])
        dither2.gfa38_ra.append(r)
        dither2.gfa38_dec.append(d)
        dither2.gfa_all_ra.append(rr)
        dither2.gfa_all_dec.append(dd)
        dither2.sky_ra.append(hdr['SKYRA'])
        dither2.sky_dec.append(hdr['SKYDEC'])
    dither2.to_np_arrays()
    return dither2


# In[144]:


def dither_plots(start, ending, name):
    dither = dither_sequence(start, ending)

    hdrs = dither.header
    dither.delete_column('header')
    dither.mjd = np.array([h['MJD-OBS'] for h in hdrs])
    #REQADC  = '24.13,27.29'        / [deg] requested ADC angles
    #for h in hdrs:
    #    print('REQADC:', h.get('REQADC', None))
    dither.reqadc1 = np.array([h.get('REQADC',(np.nan,np.nan))[0] for h in hdrs])
    dither.reqadc2 = np.array([h.get('REQADC',(np.nan,np.nan))[1] for h in hdrs])
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
    plt.title('SKY and GFA positions, near (%.1f, %.1f)' % (sra, sdec))
    plt.axis('equal');
    plt.savefig('%s-all.png' % name)
    
    plt.figure(figsize=(8,8))
    plt.plot(dither.gfa_ra, dither.gfa_dec, 'bo-', label='Mean GFA position')
    plt.plot(dither.sky_ra, dither.sky_dec, 'rx-', label='SKY_RA,SKY_DEC header')
    plt.legend()
    plt.axis('equal')
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
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
    plt.title('Mean GFA position MINUS SKY_RA,SKY_DEC')
    plt.xlabel('delta-RA (arcsec)')
    plt.ylabel('delta-Dec (arcsec)')
    plt.savefig('%s-diff.png' % name)    
    
    return dither


# In[145]:


dither_plots(27690, 27775+1, 'dither1')
#dither_plots(27781, 27822+1, 'dither2')
#dither_plots(27825, 27866+1, 'dither3')
#dither_plots(27869, 27910+1, 'dither4')
#dither_plots(27913, 27954+1, 'dither5')


# In[177]:


#dither1 = dither_plots(27690, 27775+1, 'dither1')
dither = dither_plots(27825, 27866+1, 'dither3')


# In[112]:


dither.about()


# In[178]:


from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz
from astropy.coordinates import EarthLocation
from astropy import units as u

from astropy.utils import iers
iers.conf.auto_download = False  

site = EarthLocation.of_site('kpno')
dither.sky_alt = np.zeros(len(dither))
dither.sky_az  = np.zeros(len(dither))
dither.gfa_alt = np.zeros(len(dither))
dither.gfa_az  = np.zeros(len(dither))
for i,d in enumerate(dither):
    time = Time(d.mjd, format='mjd', scale='tai')
    coords = SkyCoord(d.sky_ra, d.sky_dec, unit='deg')
    altaz = coords.transform_to(AltAz(obstime=time, location=site))
    #print('alt,az', altaz)
    dither.sky_alt[i] = altaz.alt.to_value(unit=u.deg)
    dither.sky_az [i] = altaz.az.to_value(unit=u.deg)

    coords = SkyCoord(d.gfa_ra, d.gfa_dec, unit='deg')
    altaz = coords.transform_to(AltAz(obstime=time, location=site))
    dither.gfa_alt[i] = altaz.alt.to_value(unit=u.deg)
    dither.gfa_az [i] = altaz.az.to_value(unit=u.deg)
    


# In[179]:


plt.figure(figsize=(8,8))
plt.plot(dither.sky_az, dither.sky_alt, 'bo', label='SKY')
plt.plot(dither.gfa_az, dither.gfa_alt, 'ro', label='GFA');
plt.legend()
plt.xlabel('Az (deg)')
plt.ylabel('Alt (deg)')
plt.savefig('dither3-altaz.png');


# In[180]:


plt.figure(figsize=(8,8))
plt.plot(dither.gfa_az - dither.sky_az, dither.gfa_alt - dither.sky_alt, 'ko');
for i in [0,-1]:
    plt.text((dither.gfa_az - dither.sky_az)[i], (dither.gfa_alt - dither.sky_alt)[i], '%i'%dither.expnum[i])
plt.xlabel('GFA az - SKY az (deg)')
plt.ylabel('GFA alt - SKY alt (deg)')
plt.savefig('dither3-daltaz.png')


# In[181]:


plt.figure(figsize=(8,8))
cosdec = np.cos(np.deg2rad(dither.sky_alt))
x = (dither.gfa_az - dither.sky_az)*cosdec * 3600.
y = (dither.gfa_alt - dither.sky_alt) * 3600.
plt.plot(x, y, 'ko')
for i in [0,-1]:
    plt.text(x[i], y[i], '%i'%dither.expnum[i])
plt.xlabel('GFA az - SKY az (arcsec)')
plt.ylabel('GFA alt - SKY alt (arcsec)')
plt.axis('equal');
plt.savefig('dither3-altaz-arcsec.png')


# In[ ]:





# In[ ]:





# In[94]:


plt.plot(dither.reqadc1)
plt.plot(dither.reqadc2);


# In[96]:


plt.plot(dither.reqadc1 - dither.reqadc2);


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[65]:


plt.plot(dither1.gfa_all_ra[:,0], dither1.gfa_all_dec[:,0])
plt.plot(dither1.gfa_all_ra[:,3], dither1.gfa_all_dec[:,3]);


# In[ ]:


dither_plots(27781, 27822+1, 'dither2')


# In[ ]:


dither_plots(27825, 27866+1, 'dither3')


# In[ ]:


dither_plots(27869, 27910+1, 'dither4')


# In[ ]:


dither_plots(27913, 27954+1, 'dither5')


# In[ ]:


dither3 = dither_sequence(27825, 27866+1)
dither = dither3
name = 'dither3'


# In[ ]:


plt.figure(figsize=(8,8))
#plt.plot(dither1.gfa_ra, dither1.gfa_dec, 'bo-')
plt.plot(dither.gfa05_ra, dither.gfa05_dec, 'o-')
plt.plot(dither.gfa27_ra, dither.gfa27_dec, 'o-')
plt.plot(dither.gfa38_ra, dither.gfa38_dec, 'o-')
plt.plot(dither.sky_ra, dither.sky_dec, 'rx-')
plt.axis('equal');


# In[ ]:


plt.figure(figsize=(8,8))
plt.plot(dither.gfa_ra, dither.gfa_dec, 'bo-', label='Mean GFA position')
plt.plot(dither.sky_ra, dither.sky_dec, 'rx-', label='SKY_RA,SKY_DEC header')
plt.legend()
plt.axis('equal')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.savefig('%s-20191114.png' % name)


# In[ ]:


plt.figure(figsize=(8,8))
cosdec = np.cos(np.deg2rad(dither.sky_dec))
dR = (dither.gfa_ra  - dither.sky_ra)*cosdec*3600.
dD = (dither.gfa_dec - dither.sky_dec)*3600.
plt.plot(dR, dD, 'ko-', label='GFA - SKY')
for dr,dd,d in zip([dR[0],dR[-1]], [dD[0],dD[-1]], [dither[0], dither[-1]]):
    plt.text(dr, dd, '%i'%d.expnum, color='r', fontsize=16)
plt.legend()
plt.axis('equal')
plt.title('Mean GFA position MINUS SKY_RA,SKY_DEC')
plt.xlabel('delta-RA (arcsec)')
plt.ylabel('delta-Dec (arcsec)')
plt.savefig('%s-20191114-diff.png' % name)


# In[ ]:




