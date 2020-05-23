import sys
from collections import namedtuple, defaultdict
import sqlite3
import numpy as np
import pandas as pd
import desc.sims_truthcatalog as stc

Metadata = namedtuple('Metadata', ['expMJD', 'band'])

class LensedSneTruthCat:
    def __init__(self):
        with sqlite3.connect('/global/cfs/cdirs/descssim/DC2/'
                             'minion_1016_desc_dithered_v4_trimmed.db') \
                             as conn:
            query = 'select obsHistID, expMJD, filter from summary'
            cursor = conn.execute(query)
            self.md = dict()
            for visit, expMJD, band in cursor:
                self.md[visit] = Metadata(expMJD, band)

        self.truth_cat = sqlite3.connect('/global/cfs/cdirs/descssim/DC2/'
                                         'Run3.0i/truth_tables/'
                                         'updated_lensed_sne_truth.db')
    def __call__(self, visit):
        mjd, band = self.md[visit]
        query = f'''select unique_id, ra, dec, redshift, t_delay,
                    magnification, t0, x0, x1, c,
                    av_mw, rv_mw, lens_cat_sys_id from lensed_sne'''
        cursor = self.truth_cat.execute(query)
        data = defaultdict(list)
        for row in cursor:
            (unique_id, ra, dec, z, t_delay, magnification,
             t0, x0, x1, c, av, rv, lens_cat_sys_id) = row
            sp_factory = stc.SNSynthPhotFactory(z=z, t0=t0, x0=x0, x1=x1,
                                                c=c, snra=ra, sndec=dec)
            synth_phot = sp_factory.create(mjd - t_delay)
            data['id'].append(unique_id)
            data['ra'].append(ra)
            data['dec'].append(dec)
            data['flux'].append(synth_phot.calcFlux(band))
            data['lens_sys_id'].append(lens_cat_sys_id)
        return pd.DataFrame(data=data)


if __name__ == '__main__':
    lensed_sne_truth_cat = LensedSneTruthCat()
    visit = int(sys.argv[1])
    df = lensed_sne_truth_cat(visit)
    df.to_pickle(f'lensed_sne_fluxes_v{visit}.pkl')
