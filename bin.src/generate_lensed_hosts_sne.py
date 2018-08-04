import numpy as np
import os
import pylab as pl
import argparse
import subprocess as sp
import astropy.io.fits as pyfits
import pandas as pd
import scipy.special as ss
import om10_lensing_equations as ole

data_dir = os.path.join(os.environ['SIMS_GCRCATSIMINTERFACE_DIR'], 'data')
twinkles_data_dir = os.path.join(os.environ['TWINKLES_DIR'], 'data')
outdefault = os.path.join(data_dir,'outputs')

parser = argparse.ArgumentParser(description='The location of the desired output directory')
parser.add_argument("--outdir", dest='outdir1', type=str, default = outdefault,
                    help='Output location for FITS stamps')
args = parser.parse_args()
outdir = args.outdir1

def random_location(Reff_src, qs, phs, ns):
    """Sample a random (x, y) location from the surface brightness
    profile of the galaxy. The input parameters are Sersic parameters for the host galaxy.
    Parameters: 
    -----------
    Reff_src: float 
        the effective radius in arcseconds, the radius within which half of the light is contained
    qs: float
        axis ratio of the source, b/a
    phs: float
        position angle of the galaxy in degrees
    ns: int
        Sersic index 
    """
    phs_rad = np.deg2rad(phs-90)

    bn = ss.gammaincinv(2. * ns, 0.5)
    z = np.random.random_sample()
    x = ss.gammaincinv(2. * ns, z)
    R = (x / bn)**ns * Reff_src
    theta = np.random.random() * 2 * np.pi
    
    xp, yp = R * np.cos(theta), R * np.sin(theta)
    xt = xp * np.sqrt(qs)
    yt = yp / np.sqrt(qs)
    dx, dy = np.linalg.solve([[np.cos(phs_rad), np.sin(phs_rad)],
                             [-np.sin(phs_rad), np.cos(phs_rad)]],
                             [xt, yt])
    return dx, dy


def check_random_locations():
    """Defines a random location to compare to"""
    
    npoints = 100000
    Reff_disk = 0.2
    qs_disk = 0.3
    phs_disk = 8.
    ns_disk = 1.0

    x_d = np.zeros(npoints)
    y_d = np.zeros(npoints)

    for i in xrange(npoints):
        x_d[i], y_d[i] = random_location(Reff_disk, qs_disk, phs_disk, ns_disk)


    bsz = 5.0
    nnn = 1000  # number of pixels per side
    dsx = bsz/nnn
    xi1, xi2 = ole.make_r_coor(nnn, dsx)

    src_disk = ole.sersic_2d(xi1,xi2,0.0,0.0,Reff_disk,qs_disk,phs_disk,ns_disk)
    src_disk_norm = src_disk/(np.sum(src_disk)*dsx*dsx)

    src_disk_px = np.sum(src_disk, axis=1)
    src_disk_norm_px = src_disk_px/(np.sum(src_disk_px)*dsx)

    src_disk_py = np.sum(src_disk, axis=0)
    src_disk_norm_py = src_disk_py/(np.sum(src_disk_py)*dsx)

    xsteps = xi1[:,0]
    #---------------------------------------------------------------

    from matplotlib.ticker import NullFormatter

    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    pl.figure(1, figsize=(8, 8))

    axScatter = pl.axes(rect_scatter)
    axHistx = pl.axes(rect_histx)
    axHisty = pl.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x_d, y_d)
    axScatter.contour(xi1, xi2, src_disk, colors=['k',])

    # now determine nice limits by hand:
    binwidth = 0.02
    xymax = max(np.max(np.abs(x_d)), np.max(np.abs(y_d)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    axScatter.set_xlim((-lim, lim))
    axScatter.set_ylim((-lim, lim))

    bins = np.arange(-lim, lim + binwidth, binwidth)
    axHistx.hist(x_d, bins=bins, density=1)
    axHistx.plot(xsteps, src_disk_norm_px, 'k-')

    axHisty.hist(y_d, bins=bins, density=1,orientation='horizontal')
    axHisty.plot(src_disk_norm_py, xsteps, 'k-')

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    
    return 0

def load_in_data_sne():
    """
    Reads in catalogs of host galaxy bulge and disk as well as om10 lenses
    """

    sne_host_bulge = pd.read_csv(os.path.join(data_dir,'sne_host_bulge.csv.gz'))
    sne_host_disk = pd.read_csv(os.path.join(data_dir,'sne_host_disk.csv.gz'))
    
    idx = sne_host_bulge['image_number'] == 0
    shb_purged = sne_host_bulge[:][idx]
    shd_purged = sne_host_disk[:][idx]

    sne_lens_cats = pd.read_csv(os.path.join(twinkles_data_dir,
                                             'dc2_sne_cat.csv'))
    idx = sne_lens_cats['imno'] == 0
    
    slc_purged = sne_lens_cats[:][idx]

    return slc_purged, shb_purged, shd_purged


def create_cats_sne(index, hdu_list, ahb_list, ahd_list):
    """
    Takes input catalogs and isolates lensing parameters as well as ra and dec of lens     
    Parameters: 
        -----------
    index: int
        Index for pandas data frame
    hdu_list:
        row of data frame that contains lens parameters
    ahb_list:
        row of data frame that contains lens galaxy parameters for the galactic bulge
    ahd_list:
        row of data frame that contains lens galaxy parameters for the galactic disk """
  
    twinkles_ID = ahb_list['twinkles_system'][index]
    UID_lens = ahb_list['uniqueId_lens'][index]
    Ra_lens = ahb_list['raPhoSim_lens'][index]
    Dec_lens = ahb_list['decPhoSim_lens'][index]
        
    idx = hdu_list['twinkles_sysno'] == twinkles_ID
    lid = hdu_list['sysno'][idx].values[0]
    xl1 = 0.0
    xl2 = 0.0
    vd = hdu_list['sigma'][idx].values[0]   # needed from OM10
    zd = hdu_list['zl'][idx].values[0]
    ql  = 1.0 - hdu_list['e'][idx].values[0]
    phi= hdu_list['theta_e'][idx].values[0]

    ys1 = hdu_list['snx'][idx].values[0]    # needed more discussion
    ys2 = hdu_list['sny'][idx].values[0]      # needed more discussion

    ext_shr = hdu_list['gamma'][idx].values[0]
    ext_phi = hdu_list['theta_gamma'][idx].values[0]

    #----------------------------------------------------------------------------
    lens_cat = {'xl1'        : xl1,
                'xl2'        : xl2,
                'ql'         : ql,
                'vd'         : vd,
                'phl'        : phi,
                'gamma'      : ext_shr,
                'phg'        : ext_phi,
                'zl'         : zd,
                'twinklesid' : twinkles_ID,
                'lensid'     : lid,
                'index'      : index,
                'UID_lens'   : UID_lens,
                'Ra_lens'    : Ra_lens,
                'Dec_lens'   : Dec_lens}

    #----------------------------------------------------------------------------
    mag_src_d = ahd_list['phosimMagNorm'][index]
    qs_d = ahd_list['minorAxis'][index]/ahd_list['majorAxis'][index]
    Reff_src_d = np.sqrt(ahd_list['minorAxis'][index]*ahd_list['majorAxis'][index])
    phs_d = ahd_list['positionAngle'][index]
    ns_d = ahd_list['sindex'][index]
    zs_d = ahd_list['redshift'][index]
    sed_src_d = ahd_list['sedFilepath'][index]
    
    dys2, dys1 = random_location(Reff_src_d, qs_d, phs_d, ns_d)
    ys1 = hdu_list['snx'][idx].values[0] - dys1    # needed more discussion
    ys2 = hdu_list['sny'][idx].values[0] - dys2    # needed more discussion

    srcsP_disk = {'ys1'          : ys1,
                  'ys2'          : ys2,
                  'mag_src'      : mag_src_d,
                  'Reff_src'     : Reff_src_d,
                  'qs'           : qs_d,
                  'phs'          : phs_d,
                  'ns'           : ns_d,
                  'zs'           : zs_d,
                  'sed_src'      : sed_src_d,
                  'components'   : 'disk'}
    
    #----------------------------------------------------------------------------

    mag_src_b = ahb_list['phosimMagNorm'][index]
    qs_b = ahb_list['minorAxis'][index]/ahb_list['majorAxis'][index]
    Reff_src_b = np.sqrt(ahb_list['minorAxis'][index]*ahb_list['majorAxis'][index])
    phs_b = ahb_list['positionAngle'][index]
    ns_b = ahb_list['sindex'][index]
    zs_b = ahb_list['redshift'][index]
    sed_src_b = ahb_list['sedFilepath'][index]
    
    srcsP_bulge = {'ys1'          : ys1,
                   'ys2'          : ys2,
                   'mag_src'      : mag_src_b,
                   'Reff_src'     : Reff_src_b,
                   'qs'           : qs_b,
                   'phs'          : phs_b,
                   'ns'           : ns_b,
                   'zs'           : zs_b,
                   'sed_src'      : sed_src_b,                         
                   'components'   : 'bulge'}
    

    return lens_cat, srcsP_bulge, srcsP_disk


def lensed_sersic_2d(xi1, xi2, yi1, yi2, source_cat, lens_cat):
    #Defines a magnitude of lensed host galaxy using 2d Sersic profile 
    #----------------------------------------------------------------------
    ysc1     = source_cat['ys1']        # x position of the source, arcseconds
    ysc2     = source_cat['ys2']        # y position of the source, arcseconds
    mag_tot  = source_cat['mag_src']    # total mitude of the source
    Reff_arc = source_cat['Reff_src']   # Effective Radius of the source, arcseconds
    qs       = source_cat['qs']         # axis ratio of the source, b/a
    phs      = source_cat['phs']        # orientation of the source, degree
    ndex     = source_cat['ns']         # index of the source

    #----------------------------------------------------------------------
    g_limage = ole.sersic_2d(yi1,yi2,ysc1,ysc2,Reff_arc,qs,phs,ndex)
    g_source = ole.sersic_2d(xi1,xi2,ysc1,ysc2,Reff_arc,qs,phs,ndex)

    mag_lensed = mag_tot - 2.5*np.log(np.sum(g_limage)/np.sum(g_source))

    return mag_lensed, g_limage


def generate_lensed_host(xi1, xi2, lens_P, srcP_b, srcP_d):
    """Does ray tracing of light from host galaxies using a non-singular isothermal ellipsoid profile.  
    Ultimately writes out a FITS image of the result of the ray tracing.      """
    dsx = 0.01
    xlc1 = lens_P['xl1']                # x position of the lens, arcseconds
    xlc2 = lens_P['xl2']                # y position of the lens, arcseconds
    rlc  = 0.0                          # core size of Non-singular Isothermal Ellipsoid
    vd   = lens_P['vd']                 # velocity dispersion of the lens
    zl   = lens_P['zl']                 # redshift of the lens
    zs   = srcP_b['zs']                 # redshift of the source
    rle  = ole.re_sv(vd, zl, zs)        # Einstein radius of lens, arcseconds.
    ql   = lens_P['ql']                 # axis ratio b/a
    le   = ole.e2le(1.0 - ql)           # scale factor due to projection of ellpsoid
    phl  = lens_P['phl']                # position angle of the lens, degree
    eshr = lens_P['gamma']              # external shear
    eang = lens_P['phg']                # position angle of external shear
    ekpa = 0.0                          # external convergence

    ximg, yimg = cross_check_with_lensed_sne(lens_P['twinklesid'])

    #----------------------------------------------------------------------
    ai1, ai2 = ole.alphas_sie(xlc1, xlc2, phl, ql, rle, le, eshr, eang, ekpa, xi1, xi2)

    yi1 = xi1 - ai1
    yi2 = xi2 - ai2
    #----------------------------------------------------------------------------

    lensed_mag_b, lensed_image_b = lensed_sersic_2d(xi1,xi2,yi1,yi2,srcP_b,lens_P)

    os.makedirs(os.path.join(outdir,'sne_lensed_bulges'), exist_ok=True)

    fits_limg_b = os.path.join(outdir,'sne_lensed_bulges/') + str(lens_P['UID_lens']) + "_" + str(lensed_mag_b) + "_bulge.fits"#\

    pyfits.writeto(fits_limg_b, lensed_image_b.astype("float32"), overwrite=True)

    # ----------------------------------------------------------------------------

    lensed_mag_d, lensed_image_d = lensed_sersic_2d(xi1,xi2,yi1,yi2,srcP_d,lens_P)

    os.makedirs(os.path.join(outdir,'sne_lensed_disks'), exist_ok=True)

    fits_limg_d = os.path.join(outdir,'sne_lensed_disks/') + str(lens_P['UID_lens']) + "_" + str(lensed_mag_b) + "_disk.fits" #\

    pyfits.writeto(fits_limg_d, lensed_image_d.astype("float32"), overwrite=True)

    return 0


def cross_check_with_lensed_sne(twinkles_ID):
    """ stack the lensed hosts and lensed points to verify the calculation
     make some plots."""
    sne_lens_cats = pd.read_csv(os.path.join(twinkles_data_dir,
                                             'dc2_sne_cat.csv'))
    ximgs = np.zeros((5))
    yimgs = np.zeros((5))
        
    idx = sne_lens_cats['twinkles_sysno'] == twinkles_ID
        
    imgnos = sne_lens_cats['imno'][idx]
    nimgs = len(sne_lens_cats['imno'][idx])
    
    ximgs[:nimgs] = sne_lens_cats['x'][idx].values
    yimgs[:nimgs] = sne_lens_cats['y'][idx].values
    
    return ximgs, yimgs

if __name__ == '__main__':

    dsx = 0.01  # pixel size per side, arcseconds
    nnn = 1000  # number of pixels per side
    xi1, xi2 = ole.make_r_coor(nnn, dsx)

    hdulist, ahb, ahd = load_in_data_sne()

    message_row = 0
    message_freq = 50
    for i, row in ahb.iterrows():
        if i >= message_row:
            print ("working on system ", i , "of", max(ahb.index))
            message_row += message_freq
        lensP, srcPb, srcPd = create_cats_sne(i, hdulist, ahb, ahd)
        generate_lensed_host(xi1, xi2, lensP, srcPb, srcPd)      
