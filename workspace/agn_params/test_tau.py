import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os
import time

from agn_param_module import M_i_from_L_Mass
from agn_param_module import log_Eddington_ratio
from agn_param_module import k_correction
from agn_param_module import tau_from_params
from lsst.sims.photUtils import BandpassDict, Sed
from lsst.utils import getPackageDir
from lsst.sims.photUtils import CosmologyObject

from test_m_i import make_histogram

if __name__ == "__main__":

    dtype = np.dtype([('bhmass', float), ('accretion_rate', float),
                      ('redshift', float)])

    data = np.genfromtxt('data/proto_dc2_bh_params.txt', dtype=dtype)
    print('max z  %e' % data['redshift'].max())
    valid = np.where(np.logical_and(data['bhmass']!=0.0,
                                    data['accretion_rate']!=0.0))

    data = data[valid]
    log_mbh = np.log10(data['bhmass'])

    # mass cut
    valid = np.where(log_mbh>=7.0)

    data = data[valid]
    log_mbh = log_mbh[valid]

    log_rat = log_Eddington_ratio(data['bhmass'], data['accretion_rate'])
    abs_mag = M_i_from_L_Mass(log_rat, log_mbh)
    cosmo = CosmologyObject(H0=71.0, Om0=0.265)
    DM = cosmo.distanceModulus(redshift=data['redshift'])

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    bp = bp_dict['i']
    sed_dir = os.path.join(getPackageDir('sims_sed_library'),
                           'agnSED')
    sed_name = os.path.join(sed_dir, 'agn.spec.gz')
    if not os.path.exists(sed_name):
        raise RuntimeError('\n\n%s\n\nndoes not exist\n\n' % sed_name)
    base_sed = Sed()
    base_sed.readSED_flambda(sed_name)
    z_grid = np.arange(0.0, data['redshift'].max(), 0.01)
    k_grid = np.zeros(len(z_grid),dtype=float)

    for i_z, zz in enumerate(z_grid):
        ss = Sed(flambda=base_sed.flambda, wavelen=base_sed.wavelen)
        ss.redshiftSED(zz, dimming=True)
        k = k_correction(ss, bp, zz)
        k_grid[i_z] = k

    obs_mag = np.zeros(len(data), dtype=float)
    k_arr = np.interp(data['redshift'], z_grid, k_grid)
    print('calculating observed mag')
    t_start = time.time()
    n_observable = 0
    for i_obj in range(len(data)):
       local_obs_mag = abs_mag[i_obj] + DM[i_obj] + k_arr[i_obj]
       obs_mag[i_obj] = local_obs_mag

    observable = np.where(obs_mag<=24.0)
    observable = np.where(np.logical_and(abs_mag<-26.0, abs_mag>-27.0))
    rng = np.random.RandomState(8812)
    tau = tau_from_params(data['redshift'], abs_mag, data['bhmass'], rng=rng)

    xx, yy = make_histogram(np.log10(tau[observable]), 0.01,
                            mode='normalized')

    plt.figsize=(30,30)
    plt.plot(xx,yy)
    plt.xlabel('log10(tau in days)')
    plt.ylabel('dN/dlog10(tau)')
    plt.savefig('tau_distribution.png')
    print(len(tau))
    print(len(tau[observable]))
