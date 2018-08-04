import os
import numpy as np
import pylab as pl
from scipy.interpolate import interp1d
#--------------------------------------------------------------------
from astropy.cosmology import WMAP7 as p15

vc = 2.998e5 #km/s
G = 4.3011790220362e-09 # Mpc/h (Msun/h)^-1 (km/s)^2
apr =  206264.8        # 1/1^{''}

data_dir = os.path.join(os.environ['SIMS_GCRCATSIMINTERFACE_DIR'], 'data')

def Dc(z):
    res = p15.comoving_distance(z).value*p15.h
    return res

def Dc2(z1,z2):
    Dcz1 = (p15.comoving_distance(z1).value*p15.h)
    Dcz2 = (p15.comoving_distance(z2).value*p15.h)
    res = Dcz2-Dcz1+1e-8
    return res

def re_sv(sigmav, z1, z2):
    res = 4.0*np.pi*(sigmav/vc)**2.0*Dc2(z1, z2)/Dc(z2)*apr
    return res

# #--------------------------------------------------------------------
# # Calculate Einstein Radius according to Velocity Dispersion
# #
# def re_sv(sv,z1,z2):
    # Da_s = p15.angular_diameter_distance(z2).value
    # Da_ls = p15.angular_diameter_distance_z1z2(z1,z2).value
    # res = 4.0*np.pi*(sv**2.0/(const.c.value/1e3)**2.0)*Da_ls/Da_s*apr
    # return res
#--------------------------------------------------------------------

def e2le(e_in):

    e_tmp,lef_tmp = np.loadtxt(os.path.join(data_dir, "ell_lef.dat"),
                               comments='#',usecols=(0,1),unpack=True)
    f1 = interp1d(e_tmp, lef_tmp, kind='linear')

    return f1(e_in)

def make_r_coor(nc,dsx):
    bsz = nc*dsx
    x1 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2,x1 = np.meshgrid(x1,x2)
    return x1,x2


def make_c_coor(nc,dsx):
    bsz = nc*dsx
    x1,x2 = np.mgrid[0:(bsz-dsx):nc*1j,0:(bsz-dsx):nc*1j]-bsz/2.0+dsx/2.0
    return x1,x2


def kappa_sie(x0, y0, theta, ql, re, le, x, y):
    tr = np.pi * (theta / 180.0)  # + np.pi / 2.0

    cs = np.cos(tr)
    sn = np.sin(tr)

    sx = x - x0
    sy = y - y0

    sx_r = sx * cs + sy * sn
    sy_r = -sx * sn + sy * cs

    res = re*le/2.0*1.0/np.sqrt(sx_r*sx_r*ql + sy_r*sy_r/ql)
    return res


def alphas_sie(x0, y0, theta, ql, re, le, ext_shears, ext_angle, ext_kappa, x, y):  # SIE lens model
    tr = np.pi * (theta / 180.0)   + np.pi / 2.0

    sx = x - x0
    sy = y - y0

    cs = np.cos(tr)
    sn = np.sin(tr)

    sx_r = sx * cs + sy * sn
    sy_r = -sx * sn + sy * cs

    eql = np.sqrt(ql / (1.0 - ql**2.0))
    psi = np.sqrt(sx_r**2.0 * ql + sy_r**2.0 / ql)
    dx_tmp = (re * eql * np.arctan( sx_r / psi / eql))
    dy_tmp = (re * eql * np.arctanh(sy_r / psi / eql))

    dx = dx_tmp * cs - dy_tmp * sn
    dy = dx_tmp * sn + dy_tmp * cs

    # external shear
    tr2 = np.pi * (ext_angle / 180.0)
    cs2 = np.cos(2.0 * tr2)
    sn2 = np.sin(2.0 * tr2)
    dx2 = ext_shears * (cs2 * sx + sn2 * sy)
    dy2 = ext_shears * (sn2 * sx - cs2 * sy)

    # external kappa
    dx3 = ext_kappa * sx
    dy3 = ext_kappa * sy
    return dx*le + dx2 + dx3, dy*le + dy2 + dy3

#--------------------------------------------------------------------
# Rotate regular grids
#
def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x-xcen)*np.cos(phirad)+(y-ycen)*np.sin(phirad)
    ynew = (y-ycen)*np.cos(phirad)-(x-xcen)*np.sin(phirad)
    return (xnew,ynew)

#--------------------------------------------------------------------
# 2D Sersic Profile, Peak = 1.0
#
def sersic_2d(xi1,xi2,xc1,xc2,Reff_arc,ql,pha,ndex):
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    (xi1new,xi2new) = xy_rotate(xi1, xi2, xc1, xc2, pha)
    R_scale = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/Reff_arc
    res = np.exp(-bn*((R_scale)**(1.0/ndex)-1.0))
    return res

#--------------------------------------------------------------------
# Generate lensed Sersic Profile, peak = 1.0
#
def lensed_sersic_2d(xi1, xi2, yi1, yi2, source_cat, lens_cat):

    xlc1 = lens_cat['xl1']              # x position of the lens, arcseconds
    xlc2 = lens_cat['xl2']              # y position of the lens, arcseconds
    rlc  = lens_cat['rc']               # core size of Non-singular Isothermal Ellipsoid
    vd   = lens_cat['vd']               # velocity dispersion of the lens
    zl   = lens_cat['zl']               # redshift of the lens
    zs   = source_cat['zs']             # redshift of the source
    rle  = ole.re_sv(vd, zl, zs)        # Einstein radius of lens, arcseconds.
    ql   = lens_cat['ql']               # axis ratio b/a
    le   = ole.e2le(1.0 - ql)           # scale factor due to projection of ellpsoid
    phl  = lens_cat['phl']              # position angle of the lens, degree
    eshr = lens_cat['ext_shear']        # external shear
    eang = lens_cat['ext_angle']        # position angle of external shear
    ekpa = lens_cat['ext_kappa']        # external convergence

    #----------------------------------------------------------------------
    ai1, ai2 = ole.alphas_sie(xlc1, xlc2, phl, ql, rle, le,
                              esh, eang, ekpa, xi1, xi2)

    yi1 = xi1-ai1
    yi2 = xi2-ai2

    #----------------------------------------------------------------------
    ysc1     = source_cat['ys1']        # x position of the source, arcseconds
    ysc2     = source_cat['ys2']        # y position of the source, arcseconds
    mag_tot  = source_cat['mag_src']    # total magnitude of the source
    Reff_arc = source_cat['reff_src']   # Effective Radius of the source, arcseconds
    qs       = source_cat['qs']         # axis ratio of the source, b/a
    phs      = source_cat['phs']        # orientation of the source, degree
    ndex     = source_cat['index']      # index of the source

    #----------------------------------------------------------------------
    g_limage = sersic_2d(yi1,yi2,ysc1,ysc2,Reff_arc,qs,phs,ndex)
    g_source = sersic_2d(xi1,xi2,ysc1,ysc2,Reff_arc,qs,phs,ndex)

    mag_lensed = mag_tot - 2.5*np.log(np.sum(g_limage)/np.sum(g_source))

    return mag_lensed, g_limage
