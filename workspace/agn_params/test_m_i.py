import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os
import time

from agn_param_module import M_i_from_L_Mass
from agn_param_module import log_Eddington_ratio
from agn_param_module import k_correction
from lsst.sims.photUtils import BandpassDict, Sed
from lsst.utils import getPackageDir
from lsst.sims.photUtils import CosmologyObject

def make_2d_histogram(xx, yy, dx, dy):
    """
    returns indices and counts of unique points on the map
    """
    i_color1 = np.round(xx/dx).astype(int)
    i_color2 = np.round(yy/dy).astype(int)
    dex_reverse = np.array([i_color1, i_color2])
    dex_arr = dex_reverse.transpose()
    # see http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    dex_raw = np.ascontiguousarray(dex_arr).view(np.dtype((np.void, dex_arr.dtype.itemsize*dex_arr.shape[1])))
    _, unique_rows, unique_counts = np.unique(dex_raw, return_index=True, return_counts=True)

    return unique_rows, unique_counts

def plot_color_mesh(xx, yy, dx, dy, vmin=None, vmax=None):
    i_x_arr = np.round((xx-xx.min())/dx).astype(int)
    i_y_arr = np.round((yy-yy.min())/dy).astype(int)
    new_x = i_x_arr*dx
    new_y = i_y_arr*dy
    dex_list, ct_list = make_2d_histogram(new_x, new_y, dx, dy)

    if i_x_arr.min()<0 or i_y_arr.min()<0:
        raise RuntimeError('negative dex')

    x_mesh=np.arange(xx.min(),xx.max()+0.1,dx)
    y_mesh=np.arange(yy.min(),yy.max()+0.1,dy)
    x_mesh,y_mesh = np.meshgrid(x_mesh,y_mesh,indexing='xy')
    print('x mesh shape ',x_mesh.shape)
    z_mesh = np.zeros(shape=x_mesh.shape, dtype=int)
    ct_1000b = 0

    for dex, ct in zip(dex_list, ct_list):
        ix = i_x_arr[dex]
        iy = i_y_arr[dex]
        z_mesh[iy][ix] += ct

    z_mesh = np.ma.masked_where(z_mesh==0,z_mesh)
    plt.pcolormesh(x_mesh,y_mesh,z_mesh, vmin=vmin, vmax=vmax)
                   #norm=matplotlib.colors.LogNorm(vmin=1.0,
                   #                               vmax=1.2e6))
    plt.colorbar(label='sources per pixel')

    return None

if __name__ == "__main__":

    mbh_grid = np.arange(7.5, 10.5, 0.1)
    ledd_grid = np.arange(-4.0, 0.5, 0.1)
    mbh = []
    ledd = []
    for mm in mbh_grid:
        for ll in ledd_grid:
            mbh.append(mm)
            ledd.append(ll)

    ledd = np.array(ledd)
    mbh = np.array(mbh)

    mi = M_i_from_L_Mass(ledd, mbh)
    valid = np.where(ledd>=-2.0)
    mbh = mbh[valid]
    ledd = ledd[valid]
    mi = mi[valid]

    plt.figsize = (30,30)
    plt.scatter(mbh, mi, c=ledd,
                cmap=plt.get_cmap('gist_rainbow_r'))
    plt.colorbar()
    plt.xlabel('Log(Mbh/Msun)')
    plt.ylabel('M_i')
    plt.xlim(7.3,11.0)
    plt.ylim(-29.4, -22.6)
    plt.gca().invert_yaxis()
    plt.savefig('macleod_15_analogy.png')
    plt.close()

    dtype = np.dtype([('bhmass', float), ('accretion_rate', float),
                      ('redshift', float)])

    data = np.genfromtxt('data/proto_dc2_bh_params.txt', dtype=dtype)
    valid = np.where(np.logical_and(data['bhmass']!=0.0,
                                    data['accretion_rate']!=0.0))

    data = data[valid]
    log_rat = log_Eddington_ratio(data['bhmass'], data['accretion_rate'])
    log_mbh = np.log10(data['bhmass'])
    m_i = M_i_from_L_Mass(log_rat, log_mbh)

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
       local_obs_mag = m_i[i_obj] + DM[i_obj] + k_arr[i_obj]
       obs_mag[i_obj] = local_obs_mag

    print('tot mags %d' % len(obs_mag))
    observable = np.where(obs_mag<=24.0)
    print('observable mags %d' % len(observable[0]))
    exit()

    valid =np.where(m_i<0.0)

    plt.figsize=(30,30)
    plt.subplot(2,1,1)
    plot_color_mesh(log_mbh[valid], m_i[valid], 0.1, 0.1)
    plt.xlabel('log(Mbh/Msun)')
    plt.ylabel('M_i')
    plt.subplot(2,1,2)
    plot_color_mesh(log_rat[valid], m_i[valid], 0.1, 0.1)
    plt.xlabel('log(L/L_Eddington)')
    plt.ylabel('M_i')
    plt.tight_layout()
    plt.savefig('m_i_distributions.png')

