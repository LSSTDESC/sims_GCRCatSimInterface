import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from agn_param_module import M_i_from_L_Mass
from agn_param_module import log_Eddington_ratio

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

    l_edd = [-0.5, -0.5, -0.5,
                 0.15, 0.15, 0.1,
                 -0.9, -0.9, -0.9, -0.9,
                 -1.4, -1.4, -1.4,
                 -1.8, -1.8, -1.8]
    mbh = [7.7, 9.0, 9.6,
               8.4, 8.2, 7.9,
               9.6, 8.9, 8.4, 10.1,
               10.1, 9.6, 9.1,
               8.9, 9.2, 9.4]
    m_i = [-23.2, -26.5, -27.8,
               -26.8, -26.4, -26.0,
               -26.8, -25.2, -23.6, -28.5,
               -26.9, -25.2, -24.0,
               -23.2, -23.6, -24.0]

    l_edd = np.array(l_edd)
    mbh = np.array(mbh)
    m_i_test = M_i_from_L_Mass(l_edd, mbh)
    err = 0.0
    for ii in range(len(l_edd)):
        # print(m_i[ii], m_i_test[ii], np.abs(m_i[ii]-m_i_test[ii]))
        err += (m_i[ii]-m_i_test[ii])**2
    print('rms err is %e' % (np.sqrt(err/len(m_i_test))))

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

    dtype = np.dtype([('bhmass', float), ('accretion_rate', float)])
    data = np.genfromtxt('data/proto_dc2_bh_params.txt', dtype=dtype)
    valid = np.where(np.logical_and(data['bhmass']!=0.0,
                                    data['accretion_rate']!=0.0))

    data = data[valid]
    log_rat = log_Eddington_ratio(data['bhmass'], data['accretion_rate'])
    log_mbh = np.log10(data['bhmass'])
    m_i = M_i_from_L_Mass(log_rat, log_mbh)
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

