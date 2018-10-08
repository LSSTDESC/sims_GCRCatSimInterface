import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

import os
import h5py

data_name = os.path.join(os.environ['SCRATCH'], 'extincted_galaxy_fit',
                         'galaxy_fitted_a0.01.h5')
assert os.path.isfile(data_name)

out_dir = os.path.join(os.environ['SCRATCH'], 'extincted_galaxy_fit')


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
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        if i_fig==3:
            plt.xlim(-0.2, 0.3)
            plt.ylim(-0.2, 0.3)
        if i_fig==4:
            plt.xlim(-0.2, 0.2)
            plt.ylim(-0.2, 0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'galacticus_dust_color_color_comparison_10451.png'))
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
    plt.savefig(os.path.join(out_dir, 'galacticus_dust_delta_color_color_comparison_10451.png'))
    plt.close()

    plt.figure(figsize=(20,20))
    i_fig=0
    for i_bp1 in range(len(bp_list)-1):
        i_fig += 1
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        c_c = data['cosmo_%s' % bp1].value - data['cosmo_%s' % bp2].value
        a_c = data['fit_%s' % bp1].value - data['fit_%s' % bp2].value
        plt.subplot(3,2,i_fig)
        plt.hist(c_c, color='b', bins=100, normed=True)
        plt.hist(a_c, color='r', bins=100, alpha=0.6, normed=True)
        if i_fig == 1:
            labels = ['cosmoDC2', 'CatSim']
            handles = [Rectangle((0,0),1,1,color='b',ec="k"),
                       Rectangle((0,0),1,1,color='r',ec="k")]
            plt.legend(handles, labels, loc=0, fontsize=20)
        plt.xlabel('%s-%s' % (bp1, bp2), fontsize=30)
        plt.ylabel('normed histogram', fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'galacticus_dust_color_1d_dist_10451.png'))
    plt.close()

    plt.figure(figsize=(20,20))
    i_fig=0
    for i_bp1 in range(len(bp_list)-1):
        i_fig += 1
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        c_c = data['cosmo_%s' % bp1].value - data['cosmo_%s' % bp2].value
        a_c = data['fit_%s' % bp1].value - data['fit_%s' % bp2].value
        plt.subplot(3,2,i_fig)
        plt.hist(c_c-a_c, color='r', bins=100, normed=True)

        plt.xlabel('$\Delta$ %s-%s (cosmoDC2 - CatSim)' % (bp1, bp2), fontsize=30)
        plt.ylabel('normed histogram', fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        if i_fig==1:
            plt.xlim(-0.5, 0.5)
        elif i_fig==2:
            plt.xlim(-0.3, 0.2)
        elif i_fig==3:
            plt.xlim(-0.2, 0.2)
        elif i_fig==4:
            plt.xlim(-0.2, 0.2)
        elif i_fig==5:
            plt.xlim(-0.15,0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'galacticus_dust_delta_color_1d_dist_10451.png'))
    plt.close()
