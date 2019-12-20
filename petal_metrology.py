import os
import json
from pkg_resources import resource_filename
import numpy as np

from astrometry.util.fits import fits_table

from mappings import petal_id_to_gfa_num


class PetalMetrology(object):

    def __init__(self, fids, gfa_trans):
        self.fids = fids
        I = np.flatnonzero(fids.gif_num == 1)
        assert(len(I) == 1)
        self.gif1 = fids[I[0]]
        I = np.flatnonzero(fids.gif_num == 2)
        assert(len(I) == 1)
        self.gif2 = fids[I[0]]
        self.fifs = fids[fids.is_fif]

        G1,G2 = self.gif1, self.gif2
        Ti = gfa_trans
        v1 = np.array([np.mean(G1.x), np.mean(G1.y)])
        v2 = np.array([np.mean(G2.x), np.mean(G2.y)])
        vc = (v1 + v2) / 2.
        dv = v2 - v1
        p1 = np.array([np.mean(Ti.gif_1_mm_x), np.mean(Ti.gif_1_mm_y)])
        p2 = np.array([np.mean(Ti.gif_2_mm_x), np.mean(Ti.gif_2_mm_y)])
        pc = (p1 + p2) / 2.
        dp = p2 - p1
        th1 = np.arctan2(dv[1], dv[0])
        th2 = np.arctan2(dp[1], dp[0])
        dth = th2 - th1
        R = np.array([[np.cos(dth), np.sin(dth)],[-np.sin(dth), np.cos(dth)]])
        S = np.sqrt(np.sum(dv**2)) / np.sqrt(np.sum(dp**2))
        M = np.zeros((2,3), np.float32)
        M[:2,:2] = R * S
        M[:,2] = vc
        MI = np.zeros((2,3), np.float32)
        MI[:2,:2] = R.T / S
        MI[:,2] = pc
        self.pc = pc
        self.vc = vc
        self.M = M
        self.MI = MI
        for k in ['pix_x_coeffs', 'pix_y_coeffs', 'mm_x_coeffs', 'mm_y_coeffs']:
            setattr(self, k, gfa_trans.get(k))
        self.gfa = gfa_trans

        # GFA CCD bounds
        w,h = 2048, 1032
        self.ccdw, self.ccdh = w,h
        self.ccdbpx = np.array([0.5, 0.5, w+0.5, w+0.5, 0.5])
        self.ccdbpy = np.array([0.5, h+0.5, h+0.5, 0.5, 0.5])
        self.ccdbx,self.ccdby = self.gfa_pix_to_focal_mm(self.ccdbpx, self.ccdbpy)
        
    def gfa_mm_to_focal_mm(self, gfax, gfay):
        gfax = gfax.ravel()
        gfay = gfay.ravel()
        N = len(gfax)
        v = np.zeros((3,N))
        v[0,:] = gfax - self.pc[0]
        v[1,:] = gfay - self.pc[1]
        v[2,:] = 1.
        xy = np.matmul(self.M, v)
        return xy[0,:], xy[1,:]

    def focal_mm_to_gfa_mm(self, x, y):
        x = x.ravel()
        y = y.ravel()
        N = len(x)
        v = np.zeros((3,N))
        v[0,:] = x - self.vc[0]
        v[1,:] = y - self.vc[1]
        v[2,:] = 1.
        xy = np.matmul(self.MI, v)
        return xy[0,:], xy[1,:]

    def gfa_mm_to_gfa_pix(self, x, y):
        cox = self.pix_x_coeffs
        coy = self.pix_y_coeffs
        return (cox[0] + cox[1] * x + cox[2] * y,
                coy[0] + coy[1] * x + coy[2] * y)

    def gfa_pix_to_gfa_mm(self, x, y):
        cox = self.mm_x_coeffs
        coy = self.mm_y_coeffs
        return (cox[0] + cox[1] * x + cox[2] * y,
                coy[0] + coy[1] * x + coy[2] * y)

    def focal_mm_to_gfa_pix(self, x, y):
        gx,gy = self.focal_mm_to_gfa_mm(x, y)
        return self.gfa_mm_to_gfa_pix(gx, gy)

    def gfa_pix_to_focal_mm(self, x, y):
        gx,gy = self.gfa_pix_to_gfa_mm(x, y)
        return self.gfa_mm_to_focal_mm(gx, gy)


def get_petal(petal_id):
    gfa_num = petal_id_to_gfa_num[petal_id]

    # 
    datadir = resource_filename('desi_commish', 'data')
    fn = os.path.join(datadir, 'petal-metrology-json',
                      'petal%i.json' % petal_id)
    J = json.load(open(fn))

    Fids = fits_table()
    Fids.name = []
    Fids.petal_id = []
    Fids.device_loc = []
    Fids.xyz = np.zeros((len(J),4,3), np.float32)

    for i,(k,v) in enumerate(J.items()):
        Fids.name.append(k)
        Fids.petal_id.append(v['petal_id'])
        Fids.device_loc.append(v['device_loc'])
        for ipin in range(4):
            vv = v['pinhole%i' % (ipin+1)]
            Fids.xyz[i, ipin, 0] = vv['x']
            Fids.xyz[i, ipin, 1] = vv['y']
            Fids.xyz[i, ipin, 2] = vv['z']
    Fids.to_np_arrays()

    Fids.x = Fids.xyz[:,:,0]
    Fids.y = Fids.xyz[:,:,1]
    Fids.z = Fids.xyz[:,:,2]

    ## MAGIC numbers 541,542 are from DESI-0530 table "Positioners and Fiducial Locations"
    Fids.gif_num = np.array([{541:1, 542:2}.get(d,0) for d in Fids.device_loc])
    Fids.is_gif = np.array([d in [541, 542] for d in Fids.device_loc])
    Fids.is_fif = np.logical_not(Fids.is_gif)

    fn = os.path.join(datadir, 'gfa-metrology-transforms.fits')
    T = fits_table(fn)
    I = np.flatnonzero(T.gfa_num == gfa_num)
    assert(len(I) == 1)
    Ti = T[I[0]]
    
    petal = PetalMetrology(Fids, Ti)
    return petal
