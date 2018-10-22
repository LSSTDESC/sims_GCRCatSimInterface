import os
import h5py
import numpy as np
from lsst.sims.photUtils import Sed, BandpassDict
from lsst.sims.photUtils import getImsimFluxNorm
import GCRCatalogs
from GCR import GCRQuery

h_query = GCRQuery('healpix_pixel==9556')
cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')

in_dir = os.path.join(os.environ['SCRATCH'], 'sed_181019')
#in_dir = os.path.join(os.environ['SCRATCH'], 'sed_cache_181017')

assert os.path.isdir(in_dir)
in_file = os.path.join(in_dir, 'sed_cache_clever.h5')
#in_file = os.path.join(in_dir, 'test_10k.h5')
assert os.path.isfile(in_file)

qty_names = ['galaxy_id', 'redshift_true']
for bp in 'ugrizy':
    qty_names.append('mag_true_%s_lsst' % bp)
qties = cat.get_quantities(qty_names, native_filters=[h_query])

(tot_dict,
 hw_dict) = BandpassDict.loadBandpassesFromFiles()

f = h5py.File(in_file, 'r')

n_obj = len(f['galaxy_id'].value)
sorted_dexes = np.argsort(f['galaxy_id'].value)

galaxy_id = f['galaxy_id'].value[sorted_dexes]
disk_sedname = f['disk_sed'].value[sorted_dexes]
bulge_sedname = f['bulge_sed'].value[sorted_dexes]

bulge_av = f['bulge_av'].value[sorted_dexes]
bulge_rv = f['bulge_rv'].value[sorted_dexes]
disk_av = f['disk_av'].value[sorted_dexes]
disk_rv = f['disk_rv'].value[sorted_dexes]

disk_magnorm = {}
bulge_magnorm = {}
for ibp, bp in enumerate('ugrizy'):
    disk_magnorm[bp] = f['disk_magnorm'].value[ibp][sorted_dexes]
    bulge_magnorm[bp] = f['bulge_magnorm'].value[ibp][sorted_dexes]

in_the_set = np.in1d(qties['galaxy_id'], galaxy_id)
for k in qties.keys():
    qties[k]= qties[k][in_the_set]

idx = np.searchsorted(galaxy_id, qties['galaxy_id'])
np.testing.assert_array_equal(galaxy_id, qties['galaxy_id'][idx])

for k in qties.keys():
    qties[k] = qties[k][idx]

#for bp in 'ugrizy':
#    np.testing.assert_array_equal(f['cosmo_%s'%bp].value[sorted_dexes],
#                                  qties['mag_true_%s_lsst' % bp])


lsst_mags = np.zeros((6, n_obj), dtype=float)

sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']

ccm_w = None
d_maxes = np.zeros((6), dtype=float)
for i_obj in range(n_obj):
    for i_band, bp in enumerate('ugrizy'):
        disk_sed = Sed()
        disk_sed.readSED_flambda(os.path.join(sed_dir,
                                 disk_sedname[i_obj].decode('utf-8')))
        fnorm = getImsimFluxNorm(disk_sed,
                                 disk_magnorm[bp][i_obj])

        disk_sed.multiplyFluxNorm(fnorm)
        if ccm_w is None or not np.array_equal(disk_sed.wavelen, ccm_w):
            ccm_w = np.copy(disk_sed.wavelen)
            ax, bx = disk_sed.setupCCM_ab()
        disk_sed.addDust(ax, bx, A_v=disk_av[i_obj],
                         R_v=disk_rv[i_obj])

        disk_sed.redshiftSED(qties['redshift_true'][i_obj], dimming=True)

        bulge_sed = Sed()
        bulge_sed.readSED_flambda(os.path.join(sed_dir,
                                  bulge_sedname[i_obj].decode('utf-8')))
        fnorm = getImsimFluxNorm(bulge_sed,
                                 bulge_magnorm[bp][i_obj])

        bulge_sed.multiplyFluxNorm(fnorm)
        if ccm_w is None or not np.array_equal(bulge_sed.wavelen, ccm_w):
            ccm_w = np.copy(bulge_sed.wavelen)
            ax, bx = bulge_sed.setupCCM_ab()
        bulge_sed.addDust(ax, bx, A_v=bulge_av[i_obj],
                          R_v=bulge_rv[i_obj])

        bulge_sed.redshiftSED(qties['redshift_true'][i_obj], dimming=True)

        disk_flux = disk_sed.calcFlux(hw_dict[bp])
        bulge_flux = bulge_sed.calcFlux(hw_dict[bp])

        tot_mag = disk_sed.magFromFlux(disk_flux+bulge_flux)
        lsst_mags[i_band,i_obj] = tot_mag
        #d_mag = np.abs(tot_mag-qties['mag_true_%s_lsst' % bp][i_obj])
        #if d_mag>d_maxes[i_band]:
        #    d_maxes[i_band] = d_mag
        #print(d_maxes)

for i_band, bp in enumerate('ugrizy'):
    dmag = np.abs(lsst_mags[i_band,:]-qties['mag_true_%s_lsst' %bp])
    print(bp, dmag.max())
