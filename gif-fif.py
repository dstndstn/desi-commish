from astrometry.util.fits import fits_table
import json
import fitsio
from astrometry.util.util import Tan, Sip, fit_sip_wcs_py
from astrometry.util.starutil_numpy import radectoxyz, arcsec_between
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np

from platescale_ss import Tpsfunc, Rps, Mps, Sps, Tps
from petal_metrology import get_petal

for petal_id in [4, #5, 
                 6, 3, #8,
                 10, #11,
                 2, 7, #9
                 ]:
    print('Petal', petal_id)
    petal = get_petal(petal_id)

    g1x,g1y = petal.gfa_mm_to_focal_mm(petal.gfa.gif_1_mm_x, petal.gfa.gif_1_mm_y)
    g2x,g2y = petal.gfa_mm_to_focal_mm(petal.gfa.gif_2_mm_x, petal.gfa.gif_2_mm_y)

    print('Scatter in GIF1 fit positions: %.3f mm' % (np.mean(np.hypot(petal.gif1.x - g1x, petal.gif1.y - g1y))))
    print('Scatter in GIF2 fit positions: %.3f mm' % (np.mean(np.hypot(petal.gif2.x - g2x, petal.gif2.y - g2y))))

    # Now fit a SIP polynomial distortion model for how the echo22 optics projects onto GFA CCD pixels.
    # 
    # Evaluate a grid of points in CCD space and in RA,Dec, using the metrology transformations to go from GFA CCD pixels
    # to focal plane coordinates, and from there to Theta and RA,Dec.

    x0 = min([min(petal.gfa.gif_1_pix_x), min(petal.gfa.gif_2_pix_x), 0]) - 100
    y0 = min([min(petal.gfa.gif_1_pix_y), min(petal.gfa.gif_2_pix_y), 0]) - 100
    x1 = max([max(petal.gfa.gif_1_pix_x), max(petal.gfa.gif_2_pix_x), petal.ccdw]) + 100
    y1 = max([max(petal.gfa.gif_1_pix_y), max(petal.gfa.gif_2_pix_y), petal.ccdh]) + 100

    ccdgridpx, ccdgridpy = np.meshgrid(np.linspace(x0, x1, 20), np.linspace(y0, y1, 20))
    ccdgridpx = ccdgridpx.ravel()
    ccdgridpy = ccdgridpy.ravel()
    gridx, gridy = petal.gfa_pix_to_focal_mm(ccdgridpx, ccdgridpy)
    gridr = np.hypot(gridx, gridy)

    crpixx,crpixy = (petal.ccdw+1.)/2., (petal.ccdh+1.)/2.
    crx,cry = petal.gfa_pix_to_focal_mm(crpixx, crpixy)

    theta = Tpsfunc(gridr)
    gridu = theta * gridx / gridr
    gridv = theta * gridy / gridr
    crr = np.hypot(crx, cry)
    crd = Tpsfunc(crr)
    cru = crd * crx / crr
    crv = crd * cry / crr

    griddec = gridv
    gridra  = gridu / np.cos(np.deg2rad(griddec))
    starxyz = radectoxyz(gridra, griddec)
    fieldxy = np.vstack((ccdgridpx, ccdgridpy)).T
    weights = np.ones(len(gridra))
    crdec = crv[0]
    crra  = cru[0] / np.cos(np.deg2rad(crdec))
    ps = 0.2/3600.
    tan_in = Tan(crra, crdec, crpixx, crpixy, -ps, 0., 0., ps, float(petal.ccdw), float(petal.ccdh))
    sip_order = inv_order = 4
    sip = fit_sip_wcs_py(starxyz, fieldxy, weights, tan_in, sip_order, inv_order)

    hdr = fitsio.FITSHDR()
    sip.add_to_header(hdr)
    for i,(x,y) in enumerate(zip(petal.gfa.gif_1_pix_x, petal.gfa.gif_1_pix_y)):
        hdr.add_record(dict(name='GIF1X%i' % (i+1), value=x,
                            comment='Pinhole pixel pos'))
        hdr.add_record(dict(name='GIF1Y%i' % (i+1), value=y))
    for i,(x,y) in enumerate(zip(petal.gfa.gif_2_pix_x, petal.gfa.gif_2_pix_y)):
        hdr.add_record(dict(name='GIF2X%i' % (i+1), value=x))
        hdr.add_record(dict(name='GIF2Y%i' % (i+1), value=y))
    fn = 'sip-petal%i.fits' % petal_id
    fitsio.write(fn, None, header=hdr, clobber=True)
    print('Wrote', fn)
