import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from agn_param_module import M_i_from_L_Mass

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
        print(m_i[ii], m_i_test[ii], np.abs(m_i[ii]-m_i_test[ii]))
        err += (m_i[ii]-m_i_test[ii])**2
    print('err is %e' % (np.sqrt(err/len(m_i_test))))


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
    plt.savefig('m_i_v_mbh.png')

