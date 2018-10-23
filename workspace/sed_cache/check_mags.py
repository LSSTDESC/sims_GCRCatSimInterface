import os
import h5py
import numpy as np
from lsst.sims.photUtils import Sed, BandpassDict
from lsst.sims.photUtils import getImsimFluxNorm
import GCRCatalogs
from GCR import GCRQuery

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--healpix', type=int, default=10068)
    parser.add_argument('--in_file', type=str, default=None)
    parser.add_argument('--nsample', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=9166)
    parser.add_argument('--in_dir', type=str, default=None)

    args = parser.parse_args()
    rng = np.random.RandomState(args.seed)

    assert args.in_file is not None

    h_query = GCRQuery('healpix_pixel==%d' % args.healpix)
    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')

    assert os.path.isdir(args.in_dir)
    in_file = os.path.join(args.in_dir, args.in_file)
    assert os.path.isfile(in_file)

    qty_names = ['galaxy_id', 'redshift_true']
    for bp in 'ugrizy':
        qty_names.append('mag_true_%s_lsst' % bp)
    qties = cat.get_quantities(qty_names, native_filters=[h_query])

    (tot_dict,
     hw_dict) = BandpassDict.loadBandpassesFromFiles()

    with h5py.File(in_file, 'r') as f:

        sed_name_lookup = f['sed_names'].value.astype(str)

        subset = {}
        chosen_dexes = rng.choice(np.arange(len(f['galaxy_id'].value),dtype=int), size=args.nsample, replace=False)
        for k in f.keys():
            if k == 'sed_names':
                continue
            if 'magnorm' not in k:
                subset[k] = f[k].value[chosen_dexes]
            else:
                subset[k] = f[k].value[:,chosen_dexes]

    n_obj = len(subset['galaxy_id'])
    sorted_dexes = np.argsort(subset['galaxy_id'])

    galaxy_id = subset['galaxy_id'][sorted_dexes]
    disk_sedname = subset['disk_sed'][sorted_dexes]
    bulge_sedname = subset['bulge_sed'][sorted_dexes]

    bulge_av = subset['bulge_av'][sorted_dexes]
    bulge_rv = subset['bulge_rv'][sorted_dexes]
    disk_av = subset['disk_av'][sorted_dexes]
    disk_rv = subset['disk_rv'][sorted_dexes]

    disk_magnorm = {}
    bulge_magnorm = {}
    for ibp, bp in enumerate('ugrizy'):
        disk_magnorm[bp] = subset['disk_magnorm'][ibp][sorted_dexes]
        bulge_magnorm[bp] = subset['bulge_magnorm'][ibp][sorted_dexes]

    in_the_set = np.in1d(qties['galaxy_id'], galaxy_id)
    for k in qties.keys():
        qties[k]= qties[k][in_the_set]

    idx = np.searchsorted(galaxy_id, qties['galaxy_id'])
    np.testing.assert_array_equal(galaxy_id, qties['galaxy_id'][idx])

    for k in qties.keys():
        qties[k] = qties[k][idx]

    #for bp in 'ugrizy':
    #    np.testing.assert_array_equal(subset['cosmo_%s'%bp][sorted_dexes],
    #                                  qties['mag_true_%s_lsst' % bp])


    lsst_mags = np.zeros((6, n_obj), dtype=float)

    sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']

    ccm_w = None
    d_maxes = np.zeros((6), dtype=float)
    for i_obj in range(n_obj):
        for i_band, bp in enumerate('ugrizy'):
            disk_sed = Sed()
            disk_sed.readSED_flambda(os.path.join(sed_dir,
                                     sed_name_lookup[disk_sedname[i_obj]]))
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
                                      sed_name_lookup[bulge_sedname[i_obj]]))
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
