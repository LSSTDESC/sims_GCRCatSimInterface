import os
import sqlite3
import h5py
import numpy as np

from lsst.sims.utils import defaultSpecMap
import lsst.sims.photUtils as sims_photUtils
from lsst.sims.utils import ObservationMetaData
import lsst.sims.catUtils.mixins.VariabilityMixin as variability

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
    metadata_keys = ['obsHistID', 'ra', 'dec', 'rotTelPos', 'mjd']
    htmid = int(chunk[0][1])
    with h5py.File(htmid_lookup_name, 'r') as in_file:
        valid_obsid = in_file['%d' % htmid][()]
        for k in metadata_keys:
            metadata_dict[k] = in_file[k][()]

    valid_obsid = np.sort(valid_obsid)
    valid_dex = np.searchsorted(metadata_dict['obsHistID'], valid_obsid)
    np.testing.assert_array_equal(metadata_dict['obsHistID'][valid_dex],
                                  valid_obsid)

    expmjd_arr = metadata_dict['mjd'][valid_dex]

    var_gen = VariabilityGenerator(chunk)
    dmag = var_gen.applyVariability(var_gen.varParamStr,
                                    expmjd=expmjd_arr)

    print('dmag ',dmag.shape)
    print('expmjd ',expmjd_arr.shape)    


if __name__ == "__main__":

    cache_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth'
    assert os.path.isdir(cache_dir)
    sims_photUtils.cache_LSST_seds(wavelen_min=0.0,
                                   wavelen_max=1500.0,
                                   cache_dir=cache_dir)

    kplr_dummy = variability.ParametrizedLightCurveMixin()
    kplr_dummy.load_parametrized_light_curves()

    stellar_db_name = os.path.join(cache_dir, 'dc2_stellar_db.db')
    assert os.path.isfile(stellar_db_name)


    query = "SELECT simobjid, htmid_6, sedFilename, magNorm, ebv, "
    query += "varParamStr, parallax, ra, decl FROM stars "
    query += "WHERE htmid_6=8978"

    mgr = multiprocessing.Manager()
    my_lock = mgr.Lock()
    output_dict = mgr.dict()

    with sqlite3.connect(stellar_db_name) as conn:
        cursor = conn.cursor()
        data_iterator = cursor.execute(query)
        chunk = data_iterator.fetchmany(10000)
        t_start = time.time()
        do_photometry(chunk, my_lock, output_dict)
        print('that took %e seconds' % (time.time()-t_start))
