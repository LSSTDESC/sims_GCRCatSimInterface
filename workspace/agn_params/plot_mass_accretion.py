import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from agn_param_module import log_Eddington_ratio

from test_m_i import make_histogram, plot_color_mesh

if __name__ == "__main__":

    dtype = np.dtype([('bhmass', float), ('accretion_rate', float),
                      ('redshift', float)])

    data = np.genfromtxt('data/proto_dc2_bh_params.txt', dtype=dtype)
    valid = np.where(np.logical_and(data['bhmass']>=np.power(10.0,7.5),
                                    data['accretion_rate']!=0.0))

    data = data[valid]

    l_edd_rat = log_Eddington_ratio(data['bhmass'], data['accretion_rate'])


    plt.figsize = (30,30)
    plt.subplot(2,1,1)
    plot_color_mesh(np.log10(data['bhmass']), l_edd_rat,
                    0.1, 0.1)
    plt.xlabel('Log10(Mbh/Msun)')
    plt.ylabel('Log10(L/L_Eddington)')

    xx, yy = make_histogram(l_edd_rat, 0.1)
    plt.subplot(2,1,2)
    plt.plot(xx,yy)
    plt.xlabel('Log10(Mbh/Msun)')
    plt.ylabel('Log10(L/L_Eddington')
    
    plt.tight_layout()
    plt.savefig('mass_v_accretion_rate.png')
    plt.close()
