import os
import sqlite3
import h5py
import numpy as np

from lsst.sims.utils import angularSeparation
from lsst.sims.utils import getRotSkyPos
from lsst.sims.utils import defaultSpecMap
import lsst.sims.photUtils as sims_photUtils
from lsst.sims.utils import ObservationMetaData
import lsst.sims.catUtils.mixins.VariabilityMixin as variability
from lsst.sims.coordUtils import chipNameFromRaDecLSST

import multiprocessing
import time

class VariabilityGenerator(variability.StellarVariabilityModels,
                           variability.MLTflaringMixin,
                           variability.ParametrizedLightCurveMixin):

    def __init__(self, chunk):
        """
        chunk = (simobjid, htmid_6, sedFilename, magNorm, ebv, varParamStr,
                 parallax, ra, dec)
        """
        self.photParams = sims_photUtils.PhotometricParameters(nexp=1,
                                                              exptime=30)
        self._actually_calculated_columns = []
        for bp in 'ugrizy':
            self._actually_calculated_columns.append('lsst_%s' % bp)
        sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']
        self.lsstBandpassDict = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()
        self.quiescent_mags = {}
        for bp in 'ugrizy':
            self.quiescent_mags[bp] = np.zeros(len(chunk), dtype=float)

        self.parallax = np.zeros(len(chunk), dtype=float)
        self.ebv = np.zeros(len(chunk), dtype=float)
        self.varParamStr = np.empty(len(chunk), dtype=(str,300))
        self.simobjid = np.empty(len(chunk), dtype=int)
        ccm_wav = None
        for i_star, star in enumerate(chunk):
            self.simobjid[i_star] = int(star[0])
            full_file_name = os.path.join(sed_dir,
                                          defaultSpecMap[star[2]])
            spec = sims_photUtils.Sed()
            spec.readSED_flambda(full_file_name)
            fnorm = sims_photUtils.getImsimFluxNorm(spec, float(star[3]))
            spec.multiplyFluxNorm(fnorm)
            if ccm_wav is None or not np.array_equal(spec.wavelen, ccm_wav):
                ccm_wav = np.copy(spec.wavelen)
                a_x, b_x = spec.setupCCM_ab()
            ebv_val = float(star[4])
            spec.addDust(a_x, b_x, ebv=ebv_val, R_v=3.1)
            mag_list = self.lsstBandpassDict.magListForSed(spec)
            for i_bp, bp in enumerate('ugrizy'):
                self.quiescent_mags[bp][i_star] = mag_list[i_bp]
            self.ebv[i_star] = ebv_val
            self.varParamStr[i_star] = star[5]
            self.parallax[i_star] = float(star[6])

        self.parallax *= (np.pi/648000000.0) # convert from mas to radians

    def column_by_name(self, name):
        if name.startswith('quiescent'):
            bp = name[-1]
            return np.copy(self.quiescent_mags[bp])
        elif name == 'ebv':
            return np.copy(self.ebv)
        elif name=='parallax':
            return np.copy(self.parallax)
        elif name=='simobjid':
            return np.copy(self.simobjid)
        raise RuntimeError("\n\nCannot get column %s\n\n" % name)


def do_photometry(chunk,
                  obs_lock, obs_metadata_dict,
                  star_lock, star_data_dict):
    """
    make sure that chunk is all from one htmid_6

    output_dict will accumulate metadata

    query = "SELECT simobjid, htmid_6, sedFilename, magNorm, ebv, "
    query += "varParamStr, parallax, ra, decl FROM stars "
    query += "WHERE htmid_6=%d" % htmid

    """
    dummy_sed = sims_photUtils.Sed()
    data_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth'
    assert os.path.isdir(data_dir)
    htmid_lookup_name = os.path.join(data_dir, 'htmid_6_to_obsHistID_lookup.h5')
    assert os.path.isfile(htmid_lookup_name)

    metadata_dict = {}
    metadata_keys = ['obsHistID', 'ra', 'dec', 'rotTelPos', 'mjd', 'filter']
    htmid = int(chunk[0][1])
    with h5py.File(htmid_lookup_name, 'r') as in_file:
        valid_obsid = in_file['%d' % htmid][()]
        for k in metadata_keys:
            metadata_dict[k] = in_file[k][()]

    valid_obsid = np.sort(valid_obsid)
    valid_dex = np.searchsorted(metadata_dict['obsHistID'], valid_obsid)
    np.testing.assert_array_equal(metadata_dict['obsHistID'][valid_dex],
                                  valid_obsid)

    for k in metadata_dict.keys():
        metadata_dict[k] = metadata_dict[k][valid_dex]

    has_dmag = False
    while not has_dmag:
        with obs_lock:
            if obs_metadata_dict['running_dmag']>=5:
                continue
            else:
                obs_metadata_dict['running_dmag'] += 1

        t_start = time.time()
        var_gen = VariabilityGenerator(chunk)
        dmag_raw = var_gen.applyVariability(var_gen.varParamStr,
                                            expmjd=metadata_dict['mjd'])

        dmag_raw = dmag_raw.transpose([1,2,0])
        # make the columns (star, mjd, filter)
        assert dmag_raw.shape == (len(chunk), len(metadata_dict['mjd']), 6)

        dmag = np.zeros((len(chunk), len(metadata_dict['mjd'])),
                        dtype=float)

        for i_mjd in range(len(metadata_dict['mjd'])):
            dmag[:,i_mjd] = dmag_raw[:,i_mjd,metadata_dict['filter'][i_mjd]]

        del dmag_raw
        has_dmag = True
        with obs_lock:
            obs_metadata_dict['running_dmag'] -= 1
            print('running dmag %d' % obs_metadata_dict['running_dmag'])

    quiescent_fluxes = np.zeros((len(chunk),6), dtype=float)
    for i_bp, bp in enumerate('ugrizy'):
        quiescent_fluxes[:, i_bp] = dummy_sed.fluxFromMag(var_gen.column_by_name('quiescent_%s' % bp))

    t_dmag = time.time()-t_start
    print('generated dmag in %e seconds' % t_dmag)

    print('dmag ',dmag.shape)
    print('expmjd ',metadata_dict['mjd'].shape)

    star_ra = np.array([c[7] for c in chunk])
    star_dec = np.array([c[8] for c in chunk])

    obs_mask = np.zeros((len(chunk), len(metadata_dict['mjd'])),
                        dtype=bool)
    print('made mask')

    fov_radius = 2.1 # in degrees
    t_start = time.time()
    for i_obs in range(len(metadata_dict['mjd'])):
        ra = np.degrees(metadata_dict['ra'][i_obs])
        dec = np.degrees(metadata_dict['dec'][i_obs])

        ## back-of-the-envelope implementation
        ## (does not use focal plane model with chip gaps, etc.)
        #ang_dist = angularSeparation(ra, dec, star_ra, star_dec)
        #valid = ang_dist<fov_radius

        rotTel = np.degrees(metadata_dict['rotTelPos'][i_obs])
        obs_md = ObservationMetaData(pointingRA=ra, pointingDec=dec,
                                     mjd=metadata_dict['mjd'][i_obs])
        rotSky = getRotSkyPos(ra, dec, obs_md, rotTel)
        obs_md.rotSkyPos = rotSky
        chip_name = chipNameFromRaDecLSST(star_ra, star_dec,
                                          obs_metadata=obs_md).astype(str)
        valid = np.char.find(chip_name, 'None')<0

        obs_mask[:,i_obs] = valid

    for i_star in range(len(chunk)):
        if obs_mask[i_star].sum()==0:
            print('did not observe %d' % int(chunk[i_star][0]))

    t_flux = time.time()
    print('calculating dflux')
    dflux_arr = []
    for i_star in range(len(chunk)):
        star_filter = metadata_dict['filter'][obs_mask[i_star]]
        q_flux = quiescent_fluxes[i_star][star_filter]
        assert len(q_flux) == obs_mask[i_star].sum()
        q_mag = dummy_sed.magFromFlux(q_flux)
        tot_mag = q_mag + dmag[i_star,obs_mask[i_star]]
        tot_flux = dummy_sed.fluxFromMag(tot_mag)
        dflux = tot_flux-q_flux
        dflux_arr.append(dflux)
        assert len(dflux) == obs_mask[i_star].sum()
    print('dflux took %e seconds' % (time.time()-t_flux))

    with obs_lock:
        obs_metadata_dict['mjd'].append(metadata_dict['mjd'])
        obs_metadata_dict['obsHistID'].append(metadata_dict['obsHistID'])
        obs_metadata_dict['filter'].append(metadata_dict['filter'])

    with star_lock:
        local_simobjid = var_gen.column_by_name('simobjid')
        star_data_dict['simobjid'].append(local_simobjid)
        star_data_dict['ra'].append(star_ra)
        star_data_dict['dec'].append(star_dec)
        for i_bp, bp in enumerate('ugrizy'):
            star_data_dict['quiescent_%s' % bp].append(quiescent_fluxes[:,i_bp])
        for i_star in range(len(local_simobjid)):
            if len(dflux_arr[i_star]) == 0:
                continue
            star_data_dict['lc_id'].append(local_simobjid[i_star])
            star_data_dict['lc_flux'].append(dflux_arr[i_star])
            star_data_dict['lc_obsHistID'].append(metadata_dict['obsHistID'][obs_mask[i_star]])


if __name__ == "__main__":

    cache_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth'
    assert os.path.isdir(cache_dir)
    sims_photUtils.cache_LSST_seds(wavelen_min=0.0,
                                   wavelen_max=1500.0,
                                   cache_dir=cache_dir)

    kplr_dummy = variability.ParametrizedLightCurveMixin()
    kplr_dummy.load_parametrized_light_curves()

    lookup_name = os.path.join(cache_dir, 'htmid_6_to_obsHistID_lookup.h5')
    assert os.path.isfile(lookup_name)

    htmid_to_ct = {}
    with h5py.File(lookup_name, 'r') as in_file:
        for k in in_file.keys():
            try:
                htmid = int(k)
            except ValueError:
                continue

            n_obs = len(in_file['%d' % htmid][()])
            htmid_to_ct[htmid] = n_obs

    stellar_db_name = os.path.join(cache_dir, 'dc2_stellar_db.db')
    assert os.path.isfile(stellar_db_name)

    mgr = multiprocessing.Manager()
    obs_lock = mgr.Lock()
    star_lock = mgr.Lock()
    obs_metadata_dict = mgr.dict()
    star_data_dict = mgr.dict()

    obs_metadata_dict['running_dmag'] = 0
    obs_metadata_dict['mjd'] = mgr.list()
    obs_metadata_dict['obsHistID'] = mgr.list()
    obs_metadata_dict['filter'] = mgr.list()

    star_data_dict['simobjid'] = mgr.list()
    star_data_dict['ra'] = mgr.list()
    star_data_dict['dec'] = mgr.list()
    for bp in 'ugrizy':
        star_data_dict['quiescent_%s' % bp] = mgr.list()
    star_data_dict['lc_flux'] = mgr.list()
    star_data_dict['lc_obsHistID'] = mgr.list()
    star_data_dict['lc_id'] = mgr.list()

    n_threads = 20

    with sqlite3.connect(stellar_db_name) as conn:
        cursor = conn.cursor()

        ct =0
        t_start = time.time()
        for htmid in [8978]:
            chunk_size = 50000000//htmid_to_ct[htmid]
            print('chunk_size %d' % chunk_size)

            query = "SELECT simobjid, htmid_6, sedFilename, magNorm, ebv, "
            query += "varParamStr, parallax, ra, decl FROM stars "
            query += "WHERE htmid_6=%d" % htmid

            data_iterator = cursor.execute(query)
            chunk = data_iterator.fetchmany(chunk_size)
            p_list = []
            while len(chunk)>0:
                p = multiprocessing.Process(target=do_photometry,
                                      args=(chunk,
                                            obs_lock, obs_metadata_dict,
                                            star_lock, star_data_dict))
                p.start()
                p_list.append(p)
                while len(p_list)>=n_threads:
                    e_list = list([p.exitcode for p in p_list])
                    was_popped = False
                    for ii in range(len(e_list)-1,-1,-1):
                        if e_list[ii] is not None:
                            assert p_list[ii].exitcode is not None
                            p_list.pop(ii)
                            was_popped = True
                            ct += chunk_size
                    if was_popped:
                        duration = (time.time()-t_start)/3600.0
                        per = duration/ct
                        print('ran %d in %e hrs (%e)' % (ct, duration, per))

                chunk = data_iterator.fetchmany(chunk_size)

            for p in p_list:
                p.join()

        print('that took %e hrs ' % ((time.time()-t_start)/3600.0))
