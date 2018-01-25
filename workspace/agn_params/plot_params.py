"""
This script will generate the figures comparing distribution of
DC1 and DC2 AGN variability parameters in different m_i bins

Requires that get_dc1_agn_distributions.py has already been run,
and that the BH parameters from DC2 are in the data/ dir
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from desc.sims.GCRCatSimInterface import M_i_from_L_Mass
from desc.sims.GCRCatSimInterface import log_Eddington_ratio
from desc.sims.GCRCatSimInterface import k_correction
from desc.sims.GCRCatSimInterface import tau_from_params
from desc.sims.GCRCatSimInterface import SF_from_params

from lsst.utils import getPackageDir
from lsst.sims.photUtils import CosmologyObject
from lsst.sims.photUtils import Sed, BandpassDict

import numpy as np
import os
import time


def make_histogram(xx_in, dmag, cut_off=None, min_val = None, mode=None):
    if cut_off is not None:
        xx = xx_in[np.where(xx_in<=cut_off+dmag)]
    else:
        xx = xx_in
    #print xx.min(),xx.max()
    if min_val is None:
        min_val=xx.min()-dmag
    i_xx = np.round((xx-min_val)/dmag).astype(int)
    unique_ixx, ct = np.unique(i_xx, return_counts=True)

    if cut_off is not None:
        valid_out = np.where(unique_ixx*dmag+min_val<=cut_off)
        ct = ct[valid_out]
        unique_ixx =unique_ixx[valid_out]

    if mode == 'normalized':
        return unique_ixx*dmag+min_val, ct.astype(float)/float(len(xx))
    elif mode == 'cumulative':
        return unique_ixx*dmag+min_val, np.cumsum(ct.astype(float))

    return unique_ixx*dmag+min_val, ct.astype(int)


if __name__ == "__main__":

    fig_dir = 'param_figures'
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)
    else:
        if not os.path.isdir(fig_dir):
            raise RuntimeError('%s is not a dir' % fig_dir)

    dc2_dtype = np.dtype([('bhmass', float), ('accretion_rate', float),
                          ('redshift', float)])

    dc2_data = np.genfromtxt('data/proto_dc2_bh_params.txt', dtype=dc2_dtype)
    print('max z  %e' % dc2_data['redshift'].max())
    valid = np.where(np.logical_and(dc2_data['bhmass']!=0.0,
                                    dc2_data['accretion_rate']!=0.0))

    dc2_data = dc2_data[valid]

    mass_cut = np.where(dc2_data['bhmass']>=10.0**7)
    dc2_data = dc2_data[mass_cut]

    dc2_log_edd_rat = log_Eddington_ratio(dc2_data['bhmass'],
                                          dc2_data['accretion_rate'])

    dc2_abs_mag_i = M_i_from_L_Mass(dc2_log_edd_rat, np.log10(dc2_data['bhmass']))

    dc2_cosmo = CosmologyObject(H0=71.0, Om0=0.265)

    DM = dc2_cosmo.distanceModulus(redshift=dc2_data['redshift'])

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    bp_i = bp_dict['i']
    sed_dir = os.path.join(getPackageDir('sims_sed_library'),
                           'agnSED')
    sed_name = os.path.join(sed_dir, 'agn.spec.gz')
    if not os.path.exists(sed_name):
        raise RuntimeError('\n\n%s\n\nndoes not exist\n\n' % sed_name)
    base_sed = Sed()
    base_sed.readSED_flambda(sed_name)
    z_grid = np.arange(0.0, dc2_data['redshift'].max(), 0.01)
    k_grid = np.zeros(len(z_grid),dtype=float)

    for i_z, zz in enumerate(z_grid):
        ss = Sed(flambda=base_sed.flambda, wavelen=base_sed.wavelen)
        ss.redshiftSED(zz, dimming=True)
        k = k_correction(ss, bp_i, zz)
        k_grid[i_z] = k

    dc2_obs_mag_i = np.zeros(len(dc2_data), dtype=float)
    k_arr = np.interp(dc2_data['redshift'], z_grid, k_grid)
    print('calculating observed mag')
    for i_obj in range(len(dc2_data)):
       local_obs_mag = dc2_abs_mag_i[i_obj] + DM[i_obj] + k_arr[i_obj]
       dc2_obs_mag_i[i_obj] = local_obs_mag

    rng = np.random.RandomState(18231)

    dc2_tau = tau_from_params(dc2_data['redshift'],
                              dc2_abs_mag_i,
                              dc2_data['bhmass'],
                              rng=rng)

    dc2_sf = {}

    for bp in ('u', 'g', 'r', 'i', 'z', 'y'):
        dc2_sf[bp] = SF_from_params(dc2_data['redshift'], dc2_abs_mag_i,
                                    dc2_data['bhmass'],
                                    10.0*bp_dict[bp].calcEffWavelen()[0],
                                    rng=rng)

    dc1_dtype = np.dtype([('z', float), ('obs_m_i', float),
                          ('abs_m_i', float), ('tau', float),
                          ('sfu', float), ('sfg', float),
                          ('sfr', float), ('sfi', float),
                          ('sfz', float), ('sfy', float)])

    dc1_data = np.genfromtxt('data/dc1_agn_params.txt', dtype=dc1_dtype)

    z_cut = np.where(dc1_data['z']<=1.0)
    dc1_data = dc1_data[z_cut]

    plt.figsize=(30,30)
    xx1, yy1 = make_histogram(dc1_data['obs_m_i'], 0.05, mode='normalized')
    xx2, yy2 = make_histogram(dc2_obs_mag_i, 0.05, mode='normalized')
    fig_name = os.path.join(fig_dir,'obs_mag_dist.png')

    plt.subplot(2,1,1)
    l1, = plt.plot(xx1,yy1)
    l2, = plt.plot(xx2,yy2)
    plt.legend([l1,l2],['DC1','DC2'])
    plt.xlabel('m_i')
    plt.ylabel('dN/dm_i')
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)

    plt.subplot(2,1,2)
    l1, = plt.plot(xx1,yy1)
    l2, = plt.plot(xx2,yy2)
    plt.legend([l1,l2],['DC1','DC2'])
    #plt.xlabel('m_i')
    #plt.ylabel('dN/dm_i')
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.xlim(15.0,30.0)
    plt.title('zoom in', fontsize=10)
    plt.tight_layout()
    plt.savefig(fig_name)
    plt.close()

    """
    dc2_obs_mag_cut = np.where(dc2_obs_mag_i<=dc1_data['obs_m_i'].max())
    dc2_obs_mag_cut = np.where(dc2_obs_mag_i<=27.0)
    dc2_abs_mag_i = dc2_abs_mag_i[dc2_obs_mag_cut]
    dc2_tau = dc2_tau[dc2_obs_mag_cut]
    for bp in ('u', 'g', 'r', 'i', 'z', 'y'):
        dc2_sf[bp] = dc2_sf[bp][dc2_obs_mag_cut]

    dc1_obs_mag_cut = np.where(dc1_data['obs_m_i']<=27.0)
    dc1_data = dc1_data[dc1_obs_mag_cut]
    """

    M_i_max = (22.0, 23.0, 25.0, 27.0)
    M_i_min = (21.0, 22.0, 23.0, 25.0)
    plt.figsize = (30,30)
    fig_name = os.path.join(fig_dir,'agn_tau_dist.png')
    for i_bound in range(len(M_i_max)):
        mmin = M_i_min[i_bound]
        mmax = M_i_max[i_bound]
        dc1_valid = np.where(np.logical_and(dc1_data['obs_m_i']<=mmax,
                                            dc1_data['obs_m_i']>=mmin))

        dc2_valid = np.where(np.logical_and(dc2_obs_mag_i<=mmax,
                                            dc2_obs_mag_i>=mmin))



        plt.subplot(2,2,i_bound+1)
        xx1, yy1 = make_histogram(np.log10(dc1_data['tau'][dc1_valid]),0.05,
                                  mode='normalized')

        xx2, yy2 = make_histogram(np.log10(dc2_tau[dc2_valid]), 0.05,
                                  mode='normalized')

        l1, = plt.plot(xx1,yy1)
        l2, = plt.plot(xx2,yy2)
        plt.legend([l1,l2],['DC1','DC2'])
        plt.xlabel('log10(tau in days)')
        plt.ylabel('dN/dlog10(tau)')
        plt.title('%.2f <= m_i <= %.2f' % (mmin, mmax), fontsize=10)

    plt.tight_layout()
    plt.savefig(fig_name)
    plt.close

    for i_bound in range(len(M_i_max)):
        mmin = M_i_min[i_bound]
        mmax = M_i_max[i_bound]
        dc1_valid = np.where(np.logical_and(dc1_data['obs_m_i']<=mmax,
                                            dc1_data['obs_m_i']>=mmin))

        dc2_valid = np.where(np.logical_and(dc2_obs_mag_i<=mmax,
                                            dc2_obs_mag_i>=mmin))

        fig_name = os.path.join(fig_dir,
                                'agn_SF_dist_%.1f_%.1f.png' %
                                (np.abs(mmax), np.abs(mmin)))

        for i_bp, bp in enumerate(('u', 'g', 'r', 'i', 'z', 'y')):
            plt.subplot(3,2,i_bp+1)
            xx1, yy1 = make_histogram(dc1_data['sf%s' % bp][dc1_valid],0.05,
                                      mode='normalized')

            xx2, yy2 = make_histogram(dc2_sf[bp][dc2_valid], 0.05,
                                      mode='normalized')

            l1, = plt.plot(xx1,yy1)
            l2, = plt.plot(xx2,yy2)
            plt.legend([l1,l2],['DC1','DC2'])
            plt.xlabel('SF_%s (magnitudes)' % bp)
            plt.ylabel('dN/dSF_%s' % bp)
            if i_bp == 1:
                plt.title('%.2e DC1 objects; %.2e DC2 objects' %
                          (len(dc1_valid[0]), len(dc2_valid[0])),fontsize=10)
            elif i_bp == 0:
                plt.title('%.2f <= m_i <= %.2f' % (mmin, mmax), fontsize=10)

        plt.tight_layout()
        plt.savefig(fig_name)
        plt.close()

