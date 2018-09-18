import os
import numpy as np
import h5py
import copy
import GCRCatalogs
from GCR import GCRQuery
from lsst.utils import getPackageDir
from lsst.sims.photUtils import getImsimFluxNorm
from lsst.sims.photUtils import BandpassDict, Sed
from lsst.sims.photUtils import cache_LSST_seds
import time

def calc_fluxes(galaxies):

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    sed_dir = getPackageDir('sims_sed_library')

    n_obj = len(galaxies['galaxy_id'].value)

    fluxes_as_is = np.zeros((n_obj,6), dtype=float)
    fluxes_fixed = np.zeros((n_obj,6), dtype=float)

    ccm_w = None

    sed_name_arr = galaxies['sed_name'].value.astype(str)

    t_start = time.time()
    for ii in range(n_obj):
        if ii>0 and ii%1000 == 0:
            duration = (time.time()-t_start)/3600.0
            predicted = n_obj*duration/ii
            print('%d in %.2e; pred %.2e' % (ii,duration,predicted))

        ss = Sed()
        sed_name = os.path.join(sed_dir, sed_name_arr[ii])
        ss.readSED_flambda(sed_name)
        fnorm = getImsimFluxNorm(ss, galaxies['mag_norm'].value[ii])
        ss.multiplyFluxNorm(fnorm)
        if ccm_w is None or not np.array_equal(ccm_w, ss.wavelen):
            ccm_w = np.copy(ss.wavelen)
            a_x, b_x = ss.setupCCMab()

        clean_sed = copy.deepcopy(ss)
        ss.addCCMDust(a_x, b_x, A_v=galaxies['A_v'].value[ii],
                      R_v=galaxies['R_v'].value[ii])

        ss.redshiftSED(galaxies['redshift'].value[ii], dimming=True)
        local_flux = bp_dict.fluxListForSed(ss)
        fluxes_as_is[ii] = local_flux

        changed_dust = False
        new_rv = galaxies['R_v'].value[ii]
        new_av = galaxies['A_v'].value[ii]
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
            clean_sed.redshiftSED(galaxies['redshift'].value[ii],
                                  dimming=True)

            new_flux = bp_dict.fluxListForSed(clean_sed)
            fluxes_fixed[ii] = local_flux

    return fluxes_as_is, fluxes_fixed


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

    (disk_fluxes,
     disk_fluxes_fixed) = calc_fluxes(disk)

    (bulge_fluxes,
     bulge_fluxes_fixed) = calc_fluxes(bulge)

    disk.close()
    bulge.close()

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
