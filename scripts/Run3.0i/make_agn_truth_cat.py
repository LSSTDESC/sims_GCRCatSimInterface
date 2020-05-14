import sys
from collections import namedtuple, defaultdict
import json
import sqlite3
import numpy as np
import pandas as pd
import desc.sims_truthcatalog as stc

Metadata = namedtuple('Metadata', ['expMJD', 'band'])

class AgnTruthCat:
    def __init__(self):
        with sqlite3.connect('/global/projecta/projectdirs/lsst/groups/SSim/'
                             'DC2/minion_1016_desc_dithered_v4_trimmed.db') \
                             as conn:
            query = 'select obsHistID, expMJD, filter from summary'
            cursor = conn.execute(query)
            self.md = dict()
            for visit, expMJD, band in cursor:
                self.md[visit] = Metadata(expMJD, band)

        self.sed_file = stc.find_sed_file('agnSED/agn.spec.gz')
        self.truth_cat = sqlite3.connect('agn_cosmoDC2_v1.1.4_ddf.db')

    def __call__(self, visit):
        mjd, band = self.md[visit]
        mjds = np.array([mjd])
        query = f'''select galaxy_id, ra, dec, redshift, magNorm, varParamStr
                    from agn_params'''
        cursor = self.truth_cat.execute(query)
        data = defaultdict(list)
        for row in cursor:
            galaxy_id, ra, dec, z, magnorm, varParamStr = row
            unique_id = 1024*galaxy_id + 117
            pars = json.loads(varParamStr)['p']
            tau = pars[f'agn_tau_{band}']
            sf = pars[f'agn_sf_{band}']
            seed = pars['seed']
            magnorm += stc.agn_mag_norms(mjds, z, tau, sf, seed)[0]
            synth_phot = stc.SyntheticPhotometry(self.sed_file, magnorm, z)
            synth_phot.add_MW_dust(ra, dec)
            data['id'].append(unique_id)
            data['ra'].append(ra)
            data['dec'].append(dec)
            data['flux'].append(synth_phot.calcFlux(band))
        return pd.DataFrame(data=data)


if __name__ == '__main__':
    agn_truth_cat = AgnTruthCat()
    visit = 709680
    df = agn_truth_cat(visit)
    df.to_pickle(f'agn_fluxes_v{visit}.pkl')
