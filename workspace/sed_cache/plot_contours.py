import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

import GCRCatalogs
from GCR import GCRQuery

import os
import h5py

data_name = os.path.join(os.environ['SCRATCH'], 'sed_dust_grid',
                         'fit_mags_vs_cosmo_mags_10451.h5')


data_name = os.path.join(os.environ['SCRATCH'], 'soln4',
                        'soln4_fitting_9556.h5')
assert os.path.isfile(data_name)

out_dir = os.path.join(os.environ['SCRATCH'], 'soln4')


with h5py.File(data_name, 'r') as data:
    plt.figure(figsize=(20,20))
    i_fig=0

    bp_list = 'ugrizy'
    for i_bp1 in range(len(bp_list)-2):
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        bp3 = bp_list[i_bp1+2]
        i_bp2 = i_bp1+1
        i_bp3 = i_bp1+2
        i_fig += 1
        plt.subplot(2,2,i_fig)
        c_c1 = (data['cosmo_%s'% bp1].value -
                data['cosmo_%s' % bp2].value)
        c_c2 = (data['cosmo_%s'% bp2].value -
                data['cosmo_%s' % bp3].value)
        a_c1 = (data['fit_lsst'].value[i_bp1] - data['fit_lsst'].value[i_bp2])
        a_c2 = (data['fit_lsst'].value[i_bp2] - data['fit_lsst'].value[i_bp3])

        counts, xbins, ybins = np.histogram2d(c_c1, c_c2, bins=50)
        c_contours = plt.contour(counts.transpose(),
                                 extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                                 colors='b')

        counts, xbins, ybins = np.histogram2d(a_c1, a_c2, bins=50)
        a_contours = plt.contour(counts.transpose(),
                                 extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                                 colors='r', linestyle='--')

        plt.xlabel('%s-%s (observed)' % (bp1,bp2), fontsize=30)
        plt.ylabel('%s-%s (observed)' % (bp2,bp3), fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        if i_fig==3:
            plt.xlim(-0.2, 0.3)
            plt.ylim(-0.2, 0.3)
        if i_fig==4:
            plt.xlim(-0.2, 0.2)
            plt.ylim(-0.2, 0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'soln4_color_color_comparison_9556.png'))
    plt.close()

    plt.figure(figsize=(20,20))
    i_fig=0

    bp_list = 'ugrizy'
    for i_bp1 in range(len(bp_list)-2):
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        bp3 = bp_list[i_bp1+2]
        i_bp2 = i_bp1+1
        i_bp3 = i_bp1 + 2
        i_fig += 1
        plt.subplot(2,2,i_fig)
        c_c1 = (data['cosmo_%s'% bp1].value -
                data['cosmo_%s' % bp2].value)
        c_c2 = (data['cosmo_%s'% bp2].value -
                data['cosmo_%s' % bp3].value)
        a_c1 = (data['fit_lsst'].value[i_bp1] - data['fit_lsst'].value[i_bp2])
        a_c2 = (data['fit_lsst'].value[i_bp2] - data['fit_lsst'].value[i_bp3])

        counts, xbins, ybins = np.histogram2d(c_c1-a_c1, c_c2-a_c2, bins=50)
        a_contours = plt.contour(counts.transpose(),
                                 extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                                 colors='r')

        plt.xlabel('$\Delta$ %s-%s  (cosmoDC2 - CatSim; observed)' % (bp1,bp2), fontsize=30)
        plt.ylabel('$\Delta$ %s-%s (cosmoDC2 - CatSim; observed)' % (bp2,bp3), fontsize=30)
        xmin = plt.xlim()[0]
        xmax = plt.xlim()[1]
        ymin = plt.ylim()[0]
        ymax = plt.ylim()[1]
        print(xmin,xmax,ymin,ymax)
        dx = 0.25*(xmax-xmin)
        xticks = np.arange(xmin,xmax+dx,dx)
        xformat = ['%.2e' % f for f in xticks]
        dy = 0.25*(ymax-ymin)
        yticks = np.arange(ymin,ymax+dy,dy)
        yformat = ['%.2e' % f for f in yticks]
        plt.xticks(xticks,xformat,fontsize=15)
        plt.yticks(yticks,yformat,fontsize=15)
        plt.axhline(0.0,linestyle='--',color='g',linewidth=3)
        plt.axvline(0.0,linestyle='--',color='g',linewidth=3)
        #plt.xlim(-0.2, 0.2)
        #plt.ylim(-0.2, 0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'soln4_delta_color_color_comparison_9556.png'))
    plt.close()

    plt.figure(figsize=(20,20))
    i_fig=0
    for i_bp1 in range(len(bp_list)-1):
        i_fig += 1
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        i_bp2 = i_bp1+1
        c_c = data['cosmo_%s' % bp1].value - data['cosmo_%s' % bp2].value
        a_c = data['fit_lsst'].value[i_bp1] - data['fit_lsst'].value[i_bp2]
        plt.subplot(3,2,i_fig)
        plt.hist(c_c, color='b', bins=100, normed=True)
        plt.hist(a_c, color='r', bins=100, alpha=0.6, normed=True)
        if i_fig == 1:
            labels = ['cosmoDC2', 'CatSim']
            handles = [Rectangle((0,0),1,1,color='b',ec="k"),
                       Rectangle((0,0),1,1,color='r',ec="k")]
            plt.legend(handles, labels, loc=0, fontsize=20)
        plt.xlabel('%s-%s (observed)' % (bp1, bp2), fontsize=30)
        plt.ylabel('normed histogram', fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'soln4_color_1d_dist_9556.png'))
    plt.close()

    plt.figure(figsize=(20,20))
    i_fig=0
    for i_bp1 in range(len(bp_list)-1):
        i_fig += 1
        bp1 = bp_list[i_bp1]
        bp2 = bp_list[i_bp1+1]
        i_bp2 = i_bp1 + 1
        c_c = data['cosmo_%s' % bp1].value - data['cosmo_%s' % bp2].value
        a_c = data['fit_lsst'].value[i_bp1] - data['fit_lsst'].value[i_bp2]
        plt.subplot(3,2,i_fig)
        plt.hist(c_c-a_c, color='r', bins=100, normed=True)

        plt.xlabel('$\Delta$ %s-%s (cosmoDC2 - CatSim; observed)' % (bp1, bp2), fontsize=30)
        plt.ylabel('normed histogram', fontsize=30)
        plt.yticks(fontsize=20)
        """
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
        """
        xmin = plt.xlim()[0]
        xmax = plt.xlim()[1]
        dx = 0.25*(xmax-xmin)
        xticks = np.arange(xmin,xmax,dx)
        xformat = ['%.2e' % f for f in xticks]
        plt.xticks(xticks,xformat,fontsize=15)
        plt.axvline(0.0,linestyle='--',color='g',linewidth=3)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'soln4_delta_color_1d_dist_9556.png'))
    plt.close()


    n_gal = data['fit_lsst'].value.shape[1]

    disk_mag = data['disk_magnorm'].value.transpose()
    bulge_mag = data['bulge_magnorm'].value.transpose()

    disk_mean = np.mean(disk_mag, axis=1)
    bulge_mean = np.mean(bulge_mag, axis=1)

    disk_mag_std = np.std(disk_mag, axis=1)
    bulge_mag_std = np.std(bulge_mag, axis=1)
    disk_mag_spread = np.amax(disk_mag, axis=1)-np.amin(disk_mag,axis=1)
    bulge_mag_spread = np.amax(bulge_mag, axis=1)-np.amin(bulge_mag,axis=1)

    plt.figure(figsize=(20,20))
    plt.subplot(1,2,1)

    valid = np.where(np.logical_and(disk_mean<1000.0,
                                    np.isfinite(disk_mag_std)))
    plt.hist(disk_mag_std[valid], bins=100, normed=True, color='b')
    valid = np.where(np.logical_and(bulge_mean<1000.0,
                                    np.isfinite(bulge_mag_std)))
    plt.hist(bulge_mag_std[valid], bins=100, normed=True, color='r',
             alpha=0.6)
    labels = ['disks', 'bulges']
    handles = [Rectangle((0,0),1,1,color='b',ec="k"),
               Rectangle((0,0),1,1,color='r',ec="k")]
    plt.legend(handles, labels, loc=0, fontsize=20)

    plt.xlabel('std dev of magNorm across bands', fontsize=20)
    plt.ylabel('normalized histogram', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.subplot(1,2,2)
    valid = np.where(np.logical_and(disk_mean<1000.0,
                                    np.isfinite(disk_mag_spread)))
    plt.hist(disk_mag_spread[valid], bins=100, normed=True, color='b')
    valid = np.where(np.logical_and(bulge_mean<1000.0,
                                    np.isfinite(bulge_mag_spread)))
    plt.hist(bulge_mag_spread[valid], bins=100, normed=True, color='r',
             alpha=0.6)

    plt.xlabel('max-min of magNorm across bands', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'soln4_magnorm_dist.png'))
    plt.close()


    query = GCRQuery('healpix_pixel==9556')
    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')
    qties = cat.get_quantities(['galaxy_id', 'redshift'],
                               native_filters=[query])

    np.testing.assert_array_equal(qties['galaxy_id'], data['galaxy_id'].value)

    rng = np.random.RandomState(812)
    with open(os.path.join(out_dir, 'example_sed_params.txt'), 'w') as out_file:
        out_file.write('# name z Av Rv magNorms\n')
        for range_val in np.arange(0.05,1.1, 0.25):
            valid = np.where(np.logical_and(disk_mean<1000.0,
                             np.logical_and(np.isfinite(disk_mag_spread),
                             np.logical_and(disk_mag_spread>range_val-0.1,
                                            disk_mag_spread<=range_val))))

            dexes = rng.choice(valid[0], size=4, replace=False)
            for dd in dexes:
                out_file.write('%s ' % (data['disk_sed'].value[dd].decode('utf-8')))
                out_file.write('%e ' % (qties['redshift'][dd]))
                out_file.write('%e %e ' % (data['disk_av'].value[dd],
                                           data['disk_rv'].value[dd]))
                for ii in range(6):
                    out_file.write('%e ' % (data['disk_magnorm'].value[ii][dd]))
                out_file.write('\n')
