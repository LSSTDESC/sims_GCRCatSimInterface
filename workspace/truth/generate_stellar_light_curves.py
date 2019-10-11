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
        self.varParamStr = np.empty(len(chunk), dtype=(str,100))
        ccm_wav = None
        for i_star, star in enumerate(chunk):
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
        raise RuntimeError("\n\nCannot get column %s\n\n" % name)


def do_photometry(chunk, my_lock, output_dict):
    """
    make sure that chunk is all from one htmid_6
    """
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

    t_start = time.time()
    var_gen = VariabilityGenerator(chunk)
    dmag_raw = var_gen.applyVariability(var_gen.varParamStr,
                                        expmjd=metadata_dict['mjd'])

    dmag_raw = dmag_raw.transpose([1,2,0]) # make the columns (star, mjd, filter)
    assert dmag_raw.shape == (len(chunk), len(metadata_dict['mjd']), 6)

    dmag = np.zeros((len(chunk), len(metadata_dict['mjd'])),
                    dtype=float)

    for i_mjd in range(len(metadata_dict['mjd'])):
        dmag[:,i_mjd] = dmag_raw[:,i_mjd,metadata_dict['filter'][i_mjd]]

    del dmag_raw

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

    htmid = 8978
    chunk_size = 50000000//htmid_to_ct[htmid]
    print('chunk_size %d' % chunk_size)

    query = "SELECT simobjid, htmid_6, sedFilename, magNorm, ebv, "
    query += "varParamStr, parallax, ra, decl FROM stars "
    query += "WHERE htmid_6=%d" % htmid

    mgr = multiprocessing.Manager()
    my_lock = mgr.Lock()
    output_dict = mgr.dict()

    with sqlite3.connect(stellar_db_name) as conn:
        cursor = conn.cursor()
        data_iterator = cursor.execute(query)
        chunk = data_iterator.fetchmany(chunk_size)
        t_start = time.time()
        do_photometry(chunk, my_lock, output_dict)
        print('that took %e seconds' % (time.time()-t_start))
