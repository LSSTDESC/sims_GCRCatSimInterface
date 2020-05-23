import sys
from collections import namedtuple, defaultdict
import sqlite3
import numpy as np
import pandas as pd
import lsst.sims.photUtils as sims_photUtils
import desc.sims_truthcatalog as stc

Metadata = namedtuple('Metadata', ['expMJD', 'band'])

class LensedAgnTruthCat:
    def __init__(self):
        with sqlite3.connect('/global/cfs/cdirs/descssim/DC2/'
                             'minion_1016_desc_dithered_v4_trimmed.db') \
                             as conn:
            query = 'select obsHistID, expMJD, filter from summary'
            cursor = conn.execute(query)
            self.md = dict()
            for visit, expMJD, band in cursor:
                self.md[visit] = Metadata(expMJD, band)

        self.sed_file = stc.find_sed_file('agnSED/agn.spec.gz')
        self.truth_cat = sqlite3.connect('/global/cfs/cdirs/descssim/DC2/'
                                         'Run3.0i/truth_tables/'
                                         'updated_lensed_agn_truth.db')
        self.photParams = sims_photUtils.PhotometricParameters(nexp=1,
                                                               exptime=30,
                                                               gain=1)
        self.bp_dict \
            = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()

    def __call__(self, visit):
        mjd, band = self.md[visit]
        query = f'''select unique_id, ra, dec, redshift, t_delay, magnorm,
                    magnification, seed, agn_tau_{band}, agn_sf_{band},
                    av_mw, rv_mw, lens_cat_sys_id from lensed_agn'''
        cursor = self.truth_cat.execute(query)
        data = defaultdict(list)
        for row in cursor:
            (unique_id, ra, dec, z, t_delay, magnorm, magnification,
             seed, tau, sf, av, rv, lens_cat_sys_id) = row
            mjds = np.array([mjd - t_delay])
            magnorm += (stc.agn_mag_norms(mjds, z, tau, sf, seed, start_date=58350.)[0]
                        - 2.5*np.log10(np.abs(magnification)))
            synth_phot = stc.SyntheticPhotometry(self.sed_file, magnorm, z,
                                                 gAv=av, gRv=rv)
            data['id'].append(unique_id)
            data['ra'].append(ra)
            data['dec'].append(dec)
            data['flux'].append(synth_phot.calcFlux(band))
            data['adu'].append(synth_phot.sed.calcADU(self.bp_dict[band],
                                                      self.photParams))
            data['lens_sys_id'].append(lens_cat_sys_id)
        return pd.DataFrame(data=data)


if __name__ == '__main__':
    lensed_agn_truth_cat = LensedAgnTruthCat()
    visit = int(sys.argv[1])
    df = lensed_agn_truth_cat(visit)
    df.to_pickle(f'lensed_agn_fluxes_v{visit}.pkl')
