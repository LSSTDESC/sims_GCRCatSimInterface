import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import os
import h5py

data_name = os.path.join(os.environ['SCRATCH'], 'sed_dust_grid',
                         'fit_mags_vs_cosmo_mags_10451.h5')
assert os.path.isfile(data_name)

with h5py.File(data_name, 'r') as data:
    plt.figure(figsize=(20,20))
    i_fig=0

    bp_list = 'ugrizy'
    for i_bp1 in range(len(bp_list)-2):
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        bp3 = bp_list[i_bp1+2]
        i_fig += 1
        plt.subplot(2,2,i_fig)
        c_c1 = (data['cosmo_%s'% bp1].value -
                data['cosmo_%s' % bp2].value)
        c_c2 = (data['cosmo_%s'% bp2].value -
                data['cosmo_%s' % bp3].value)
        a_c1 = (data['fit_%s' % bp1].value - data['fit_%s' % bp2].value)
        a_c2 = (data['fit_%s' % bp2].value - data['fit_%s' % bp3].value)
    
        counts, xbins, ybins = np.histogram2d(c_c1, c_c2, bins=50)
        c_contours = plt.contour(counts.transpose(),
                                 extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                                 colors='b')

        counts, xbins, ybins = np.histogram2d(a_c1, a_c2, bins=50)
        a_contours = plt.contour(counts.transpose(),
                                 extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                                 colors='r')

        plt.xlabel('%s-%s' % (bp1,bp2), fontsize=30)
        plt.ylabel('%s-%s' % (bp2,bp3), fontsize=30)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)

    plt.tight_layout()
    out_dir = os.path.join(os.environ['SCRATCH'], 'sed_dust_grid')
    plt.savefig(os.path.join(out_dir, 'color_color_comparison_10451.png'))
    plt.close()

    plt.figure(figsize=(20,20))
    i_fig=0

    bp_list = 'ugrizy'
    for i_bp1 in range(len(bp_list)-2):
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        bp3 = bp_list[i_bp1+2]
        i_fig += 1
        plt.subplot(2,2,i_fig)
        c_c1 = (data['cosmo_%s'% bp1].value -
                data['cosmo_%s' % bp2].value)
        c_c2 = (data['cosmo_%s'% bp2].value -
                data['cosmo_%s' % bp3].value)
        a_c1 = (data['fit_%s' % bp1].value - data['fit_%s' % bp2].value)
        a_c2 = (data['fit_%s' % bp2].value - data['fit_%s' % bp3].value)

        counts, xbins, ybins = np.histogram2d(c_c1-a_c1, c_c2-a_c2, bins=50)
        a_contours = plt.contour(counts.transpose(),
                                 extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                                 colors='r')

        plt.xlabel('$\Delta$ %s-%s  (cosmoDC2 - CatSim)' % (bp1,bp2), fontsize=30)
        plt.ylabel('$\Delta$ %s-%s (cosmoDC2 - CatSim)' % (bp2,bp3), fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(-0.2, 0.2)
        plt.ylim(-0.2, 0.2)

    plt.tight_layout()
    out_dir = os.path.join(os.environ['SCRATCH'], 'sed_dust_grid')
    plt.savefig(os.path.join(out_dir, 'delta_color_color_comparison_10451.png'))
    plt.close()
