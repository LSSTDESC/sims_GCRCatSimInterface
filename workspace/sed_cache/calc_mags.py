import os
import numpy as np
import h5py
import copy
import multiprocessing
import GCRCatalogs
from GCR import GCRQuery
from lsst.utils import getPackageDir
from lsst.sims.photUtils import getImsimFluxNorm
from lsst.sims.photUtils import BandpassDict, Sed
from lsst.sims.photUtils import cache_LSST_seds
import time

def calc_fluxes(sed_name_arr, mag_norm, a_v, r_v, redshift):

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    sed_dir = getPackageDir('sims_sed_library')

    n_obj = len(sed_name_arr)

    fluxes_as_is = np.zeros((n_obj,6), dtype=float)
    fluxes_fixed = np.zeros((n_obj,6), dtype=float)

    ccm_w = None

    t_start = time.time()
    for ii in range(n_obj):
        if ii>0 and ii%1000 == 0:
            duration = (time.time()-t_start)/3600.0
            predicted = n_obj*duration/ii
            print('%d of %d in %.2e; pred %.2e' % (ii,n_obj,duration,predicted))

        ss = Sed()
        sed_name = os.path.join(sed_dir, sed_name_arr[ii])
        ss.readSED_flambda(sed_name)
        fnorm = getImsimFluxNorm(ss, mag_norm[ii])
        ss.multiplyFluxNorm(fnorm)
        if ccm_w is None or not np.array_equal(ccm_w, ss.wavelen):
            ccm_w = np.copy(ss.wavelen)
            a_x, b_x = ss.setupCCMab()

        clean_sed = copy.deepcopy(ss)
        ss.addCCMDust(a_x, b_x, A_v=a_v[ii],
                      R_v=r_v[ii])

        ss.redshiftSED(redshift[ii], dimming=True)
        local_flux = bp_dict.fluxListForSed(ss)
        fluxes_as_is[ii] = local_flux

        changed_dust = False
        new_rv = r_v[ii]
        new_av = a_v[ii]
        if new_rv<0.1:
            changed_dust = True
            new_rv = 0.1
            new_av = np.clip(new_av,0.0,1.0)
        elif new_rv<1.0:
            if new_av>1.0:
                changed_dust=True
                new_av=np.clip(new_av,0.0,1.0)
        elif new_av<0.0:
            changed_dust = True
            new_av=0.0

        if not changed_dust:
            fluxes_fixed[ii] = local_flux
        else:
            clean_sed.addCCMDust(a_x, b_x, A_v=new_av, R_v=new_rv)
            clean_sed.redshiftSED(redshift[ii],
                                  dimming=True)

            new_flux = bp_dict.fluxListForSed(clean_sed)
            fluxes_fixed[ii] = local_flux

    return fluxes_as_is, fluxes_fixed


def calc_fluxes_parallel(galaxy_id, sed_name, mag_norm,
                         a_v, r_v, zz, out_dict, lock):
    (fluxes,
     fluxes_fixed) = calc_fluxes(sed_name, mag_norm, a_v, r_v, zz)

    lock.acquire()
    ii = 0
    tag = '%d' % ii
    while 'fluxes_%s' % tag in out_dict:
        ii += 1
        tag = '%d' % ii
    lock.release()
    out_dict['fluxes_%s' % tag] = fluxes
    out_dict['fluxes_fixed_%s' % tag] = fluxes_fixed
    out_dict['galaxy_id_%s' % tag] = galaxy_id



if __name__ == "__main__":

    sed_dir = getPackageDir('sims_sed_library')
    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

    in_dir = os.path.join(os.environ['SCRATCH'], 'sed_cache')
    assert os.path.isdir(in_dir)

    disk_name = os.path.join(in_dir,'disk_offenders_10451.h5')
    assert os.path.isfile(disk_name)

    bulge_name = os.path.join(in_dir, 'bulge_offenders_10451.h5')
    assert os.path.isfile(bulge_name)

    cache_LSST_seds(wavelen_min=0.0, wavelen_max=1600.0,
                    cache_dir=in_dir)

    disk = h5py.File(disk_name, 'r')
    bulge = h5py.File(bulge_name, 'r')
    assert len(disk['galaxy_id'].value) == len(bulge['galaxy_id'].value)


    disk_gid = disk['galaxy_id'].value
    bulge_gid = bulge['galaxy_id'].value

    np.testing.assert_array_equal(disk_gid, bulge_gid)

    sed_name = disk['sed_name'].value.astype(str)
    mag_norm = disk['mag_norm'].value
    a_v = disk['A_v'].value
    r_v = disk['R_v'].value
    redshift = disk['redshift'].value
    gid = disk['galaxy_id'].value

    mgr = multiprocessing.Manager()
    disk_dict = mgr.dict()
    lock = multiprocessing.Lock()
    p_list = []

    n_max = len(gid)
    d_gal = n_max//60

    for i_start in range(0,n_max,d_gal):
        p = multiprocessing.Process(target=calc_fluxes_parallel,
                                    args=(gid[i_start:i_start+d_gal],
                                          sed_name[i_start:i_start+d_gal],
                                          mag_norm[i_start:i_start+d_gal],
                                          a_v[i_start:i_start+d_gal],
                                          r_v[i_start:i_start+d_gal],
                                          redshift[i_start:i_start+d_gal],
                                          disk_dict, lock))
        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()


    gid = bulge['galaxy_id'].value
    sed_name = bulge['sed_name'].value.astype(str)
    mag_norm = bulge['mag_norm'].value
    a_v = bulge['A_v'].value
    r_v = bulge['R_v'].value
    redshift = bulge['redshift'].value

    bulge_dict = mgr.dict()

    for i_start in range(0,n_max,d_gal):
        p = multiprocessing.Process(target=calc_fluxes_parallel,
                                    args=(gid[i_start:i_start+d_gal],
                                          sed_name[i_start:i_start+d_gal],
                                          mag_norm[i_start:i_start+d_gal],
                                          a_v[i_start:i_start+d_gal],
                                          r_v[i_start:i_start+d_gal],
                                          redshift[i_start:i_start+d_gal],
                                          bulge_dict, lock))
        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    disk.close()
    bulge.close()

    disk_gid = np.zeros(n_max, dtype=int)
    disk_fluxes = np.zeros((n_max, 6), dtype=float)
    disk_fluxes_fixed = np.zeros((n_max, 6), dtype=float)

    i_start = 0
    for ii in range(len(disk_dict)//3):
        local_gid = disk_dict['galaxy_id_%d' % ii]
        local_fluxes = disk_dict['fluxes_%d' % ii]
        local_fluxes_fixed = disk_dict['fluxes_fixed_%d' % ii]

        disk_gid[i_start:i_start+len(local_gid)] = local_gid
        disk_fluxes[i_start:i_start+len(local_gid)] = local_fluxes
        disk_fluxes_fixed[i_start:i_start+len(local_gid)] = local_fluxes_fixed
        i_start += len(local_gid)

    bulge_gid = np.zeros(n_max, dtype=int)
    bulge_fluxes = np.zeros((n_max,6), dtype=float)
    bulge_fluxes_fixed = np.zeros((n_max,6), dtype=float)

    i_start = 0
    for ii in range(len(bulge_dict)//3):
        local_gid = bulge_dict['galaxy_id_%d' % ii]
        local_fluxes = bulge_dict['fluxes_%d' % ii]
        local_fluxes_fixed = bulge_dict['fluxes_fixed_%d' % ii]

        bulge_gid[i_start:i_start+len(local_gid)] = local_gid
        bulge_fluxes[i_start:i_start+len(local_gid)] = local_fluxes
        bulge_fluxes_fixed[i_start:i_start+len(local_gid)] = local_fluxes_fixed
        i_start += len(local_gid)

    sorted_dex = np.argsort(disk_gid)
    disk_gid = disk_gid[sorted_dex]
    disk_fluxes = disk_fluxes[sorted_dex]
    disk_fluxes_fixed = disk_fluxes_fixed[sorted_dex]

    sorted_dex = np.argsort(bulge_gid)
    bulge_gid = bulge_gid[sorted_dex]
    bulge_fluxes = bulge_fluxes[sorted_dex]
    bulge_fluxes_fixed = bulge_fluxes_fixed[sorted_dex]

    np.testing.assert_array_equal(disk_gid, bulge_gid)

    out_name = os.path.join(in_dir, 'fluxes_10451.h5')
    out_file = h5py.File(out_name, 'w')
    out_file.create_dataset('disk_flux', data=disk_fluxes)
    out_file.create_dataset('disk_flux_fix', data=disk_fluxes_fixed)
    out_file.create_dataset('bulge_flux', data=bulge_fluxes)
    out_file.create_dataset('bluge_flux_fix', data=bulge_fluxes_fixed)
    out_file.create_dataset('galaxy_id', data=disk_gid)
    out_file.close()
    exit()

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')
    query = GCRQuery('healpix_pixel==10451')

    req = ['galaxy_id','stellar_mass_disk','stellar_mass_bulge']
    for bp in 'ugrizy':
        req.append('mag_true_%s_lsst' % bp)

    qties = cat.get_quantity(req, native_filters=[query])
