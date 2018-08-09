import numpy as np
import os
import argparse
import pylab as pl
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

def load_in_data_agn():

    """
    Reads in catalogs of host galaxy bulge and disk as well as om10 lenses
    """
    agn_host_bulge = pd.read_csv(os.path.join(data_dir,'agn_host_bulge.csv.gz'))
    agn_host_disk = pd.read_csv(os.path.join(data_dir, 'agn_host_disk.csv.gz'))

    idx = agn_host_bulge['image_number'] == 0
    ahb_purged = agn_host_bulge[:][idx]
    ahd_purged = agn_host_disk[:][idx]

    lens_list = pyfits.open(os.path.join(twinkles_data_dir,
                                         'twinkles_lenses_v2.fits'))

    return lens_list, ahb_purged, ahd_purged


def create_cats_agns(index, hdu_list, ahb_list, ahd_list):
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

    twinkles_ID = ahd['twinkles_system'][index]
    UID_lens = ahd['uniqueId_lens'][index]
    Ra_lens = ahd['raPhoSim_lens'][index]
    Dec_lens = ahd['decPhoSim_lens'][index]

    idx = hdu_list[1].data['twinklesId'] == twinkles_ID
    lid = hdu_list[1].data['LENSID'][idx][0]
    xl1 = 0.0
    xl2 = 0.0
    vd = hdu_list[1].data['VELDISP'][idx][0]
    zd = hdu_list[1].data['ZLENS'][idx][0]
    ql  = 1.0 - hdu_list[1].data['ELLIP'][idx][0]
    phi= hdu_list[1].data['PHIE'][idx][0]

    ys1 = hdu_list[1].data['XSRC'][idx][0]
    ys2 = hdu_list[1].data['YSRC'][idx][0]

    ext_shr = hdu_list[1].data['GAMMA'][idx][0]
    ext_phi = hdu_list[1].data['PHIG'][idx][0]

    ximg = hdu_list[1].data['XIMG'][idx][0]
    yimg = hdu_list[1].data['YIMG'][idx][0]
    

    #----------------------------------------------------------------------------
    lens_cat = {'xl1'        : xl1,
                'xl2'        : xl2,
                'ql'         : ql,
                'vd'         : vd,
                'phl'        : phi,
                'gamma'      : ext_shr,
                'phg'        : ext_phi,
                'zl'         : zd,
                'ximg'       : ximg,
                'yimg'       : yimg,
                'twinklesid' : twinkles_ID,
                'lensid'     : lid,
                'index'      : index,
                'UID_lens'   : UID_lens,
                'Ra_lens'    : Ra_lens,
                'Dec_lens'   : Dec_lens}
    
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
    
    #----------------------------------------------------------------------------
    mag_src_d = ahd_list['phosimMagNorm'][index]
    qs_d = ahd_list['minorAxis'][index]/ahd_list['majorAxis'][index]
    Reff_src_d = np.sqrt(ahd_list['minorAxis'][index]*ahd_list['majorAxis'][index])
    phs_d = ahd_list['positionAngle'][index]
    ns_d = ahd_list['sindex'][index]
    zs_d = ahd_list['redshift'][index]
    sed_src_d = ahd_list['sedFilepath'][index]

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

    return lens_cat, srcsP_bulge, srcsP_disk


def lensed_sersic_2d(xi1, xi2, yi1, yi2, source_cat, lens_cat):
    #Defines a magnitude of lensed host galaxy using 2d Sersic profile	
    #----------------------------------------------------------------------
    ysc1     = source_cat['ys1']        # x position of the source, arcseconds
    ysc2     = source_cat['ys2']        # y position of the source, arcseconds
    mag_tot  = source_cat['mag_src']    # total magnitude of the source
    Reff_arc = source_cat['Reff_src']   # Effective Radius of the source, arcseconds
    qs       = source_cat['qs']         # axis ratio of the source, b/a
    phs      = source_cat['phs']        # orientation of the source, degree
    ns       = source_cat['ns']         # index of the source

    #----------------------------------------------------------------------
    g_limage = ole.sersic_2d(yi1,yi2,ysc1,ysc2,Reff_arc,qs,phs,ns)
    g_source = ole.sersic_2d(xi1,xi2,ysc1,ysc2,Reff_arc,qs,phs,ns)

    mag_lensed = mag_tot - 2.5*np.log(np.sum(g_limage)/np.sum(g_source))

    return mag_lensed, g_limage


def generate_lensed_host(xi1, xi2, lens_P, srcP_b, srcP_d):
    """Does ray tracing of light from host galaxies using
    a non-singular isothermal ellipsoid profile.  
    Ultimately writes out a FITS image of the result of the ray tracing.	"""
    dsx  = 0.01
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

    #----------------------------------------------------------------------
    ai1, ai2 = ole.alphas_sie(xlc1, xlc2, phl, ql, rle, le,
                              eshr, eang, ekpa, xi1, xi2)

    yi1 = xi1 - ai1
    yi2 = xi2 - ai2
    #----------------------------------------------------------------------------

    lensed_mag_b, lensed_image_b = lensed_sersic_2d(xi1,xi2,yi1,yi2,srcP_b,lens_P)

    os.makedirs(os.path.join(outdir,'agn_lensed_bulges'), exist_ok=True)

    fits_limg_b = os.path.join(outdir,'agn_lensed_bulges/') + str(lens_P['UID_lens']) + "_" + str(lensed_mag_b) + "_bulge.fits" 
 
    pyfits.writeto(fits_limg_b, lensed_image_b.astype("float32"), overwrite=True)

    #----------------------------------------------------------------------------

    lensed_mag_d, lensed_image_d = lensed_sersic_2d(xi1,xi2,yi1,yi2,srcP_d,lens_P)

    os.makedirs(os.path.join(outdir,'agn_lensed_disks'), exist_ok=True)

    fits_limg_d = os.path.join(outdir,'agn_lensed_disks/') + str(lens_P['UID_lens']) + "_" + str(lensed_mag_d) + "_disk.fits"
 
    pyfits.writeto(fits_limg_d, lensed_image_d.astype("float32"), overwrite=True)

    return 0


if __name__ == '__main__':

    dsx = 0.01  # pixel size per side, arcseconds
    nnn = 1000  # number of pixels per side
    xi1, xi2 = ole.make_r_coor(nnn, dsx)

    hdulist, ahb, ahd = load_in_data_agn()

    message_row = 0
    message_freq = 50
    for i, row in ahb.iterrows():
        if i >= message_row:
            print ("working on system ", i , "of", max(ahb.index))
            message_row += message_freq
        lensP, srcPb, srcPd = create_cats_agns(i, hdulist, ahb, ahd)
        generate_lensed_host(xi1, xi2, lensP, srcPb, srcPd)    
