import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import os
import h5py

cosmoDC2_name = 'restframe_mags.h5'
cosmoDC2_name = os.path.join(os.environ['SCRATCH'], 'sed_dust_grid', 'restframe_10451.h5')
assert os.path.isfile(cosmoDC2_name)

catsim_name = 'sims_sed_grid.h5'
assert os.path.isfile(catsim_name)

cosmoDC2 = h5py.File(cosmoDC2_name, 'r')
catsim = h5py.File(catsim_name, 'r')

plt.figure(figsize=(20,20))
i_fig=0
bp_list = 'ugrizy'
for i_bp1 in range(len(bp_list)-2):
    bp1 = bp_list[i_bp1]
    bp2 = bp_list[i_bp1+1]
    bp3 = bp_list[i_bp1+2]
    i_fig += 1
    plt.subplot(2,2,i_fig)
    c_c1 = (cosmoDC2['Mag_true_%s_lsst_z0'% bp1].value -
            cosmoDC2['Mag_true_%s_lsst_z0' % bp2].value)
    c_c2 = (cosmoDC2['Mag_true_%s_lsst_z0'% bp2].value -
            cosmoDC2['Mag_true_%s_lsst_z0' % bp3].value)
    a_c1 = (catsim[bp1].value - catsim[bp2].value)
    a_c2 = (catsim[bp2].value - catsim[bp3].value)
    
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
plt.savefig('color_color_comparison.png')
plt.close()

stellar_mass = cosmoDC2['stellar_mass'].value
ln10 = np.log(10.0)
for log_mass in range(6,12):
    mass_min = np.exp(log_mass*ln10)
    mass_max = np.exp(log_mass*(ln10+1))
    valid = np.where(np.logical_and(stellar_mass>=mass_min,
                                    stellar_mass<mass_max))

    plt.figure(figsize=(20,20))
    i_fig=0
    for i_bp1 in range(len(bp_list)-2):
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        bp3 = bp_list[i_bp1+2]
        i_fig += 1
        plt.subplot(2,2,i_fig)
        c_c1 = (cosmoDC2['Mag_true_%s_lsst_z0'% bp1].value[valid] -
                cosmoDC2['Mag_true_%s_lsst_z0' % bp2].value[valid])
        c_c2 = (cosmoDC2['Mag_true_%s_lsst_z0'% bp2].value[valid] -
                cosmoDC2['Mag_true_%s_lsst_z0' % bp3].value[valid])
        a_c1 = (catsim[bp1].value - catsim[bp2].value)
        a_c2 = (catsim[bp2].value - catsim[bp3].value)
    
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
        if i_fig==1:
            plt.title('%d <= log10(M_star/M_sun) < %d' % (log_mass, log_mass+1), fontsize=30)

    plt.tight_layout()
    plt.savefig('color_color_comparison_%d.png' % log_mass)
    plt.close()


plt.figure(figsize=(20,20))
c = cosmoDC2['Mag_true_u_lsst_z0'].value-cosmoDC2['Mag_true_g_lsst_z0'].value
a = catsim['u'].value-catsim['g'].value
plt.hist(a,color='r',bins=200, normed=True)
plt.hist(c,color='b',bins=200,alpha=0.75,normed=True)
plt.xlabel('u-g',fontsize=30)
plt.xticks(fontsize=30)
plt.ylabel('normalized histogram', fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('ug_1d.png')
plt.close()

plt.figure(figsize=(20,20))
c = cosmoDC2['Mag_true_g_lsst_z0'].value-cosmoDC2['Mag_true_r_lsst_z0'].value
a = catsim['g'].value-catsim['r'].value
plt.hist(a,color='r',bins=200, normed=True)
plt.hist(c,color='b',bins=200,alpha=0.75,normed=True)
plt.xlabel('g-r',fontsize=30)
plt.xticks(fontsize=30)
plt.ylabel('normalized histogram', fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('gr_1d.png')
plt.close()
