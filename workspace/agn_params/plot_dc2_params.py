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

    valid = np.where(dc2_obs_mag_i<=24.0)

    tau_xx, tau_yy = make_histogram(dc2_tau[valid], 0.025)
    sf_xx, sf_yy = make_histogram(dc2_sf['i'][valid], 0.025)

    plt.figsize = (30,30)
    plt.subplot(1,2,1)
    plt.plot(tau_xx, tau_yy)
    plt.xlabel('$\\tau$ (days)')
    plt.ylabel('#')

    plt.subplot(1,2,2)
    plt.plot(sf_xx, sf_yy)
    plt.xlabel('SF_i')
    plt.ylabel('#')

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'param_distributions_180201.png'))
