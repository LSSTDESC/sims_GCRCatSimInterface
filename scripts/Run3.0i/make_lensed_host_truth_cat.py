import os
import glob
from collections import namedtuple, defaultdict
import sqlite3
from astropy.io import fits
import pandas as pd
import lsst.sims.photUtils
import desc.sims_truthcatalog as stc

Metadata = namedtuple('Metadata', ['expMJD', 'band'])

class HostsTruthCat:
    def __init__(self, image_dir):
        with sqlite3.connect('/global/projecta/projectdirs/lsst/groups/SSim/'
                             'DC2/minion_1016_desc_dithered_v4_trimmed.db') \
                             as conn:
            query = 'select obsHistID, expMJD, filter from summary'
            cursor = conn.execute(query)
            self.md = dict()
            for visit, expMJD, band in cursor:
                self.md[visit] = Metadata(expMJD, band)

        self.truth_cat = sqlite3.connect('../truth_tables/host_truth.db')
        self.image_dir = image_dir

    def __call__(self, visit):
        _, band = self.md[visit]
        phot_params = lsst.sims.photUtils.PhotometricParameters(
            exptime=30, nexp=1, gain=1, readnoise=0, darkcurrent=0,
            bandpass=band)
        data = defaultdict(list)

        df_agns = pd.read_sql('select * from agn_hosts', self.truth_cat)
        df_sne = pd.read_sql('select * from sne_hosts', self.truth_cat)

        for component in ('bulge', 'disk'):
            data = self._process_df(df_agns, band, data, component, 'agn',
                                    phot_params=phot_params)
            data = self._process_df(df_sne, band, data, component, 'sne',
                                    phot_params=phot_params)
        return pd.DataFrame(data=data)

    def _get_mag_norms(self, component, host_type, band):
        files = glob.glob(os.path.join(self.image_dir,
                                       f'{host_type}_lensed_{component}s',
                                       '*.fits'))
        mag_norms = dict()
        filter_ = list('ugrizy').index(band)
        for item in files:
            with fits.open(item) as hdus:
                header = hdus[0].header
                mag_norms[header['LENS_ID']] = header[f'MAGNORM{band.upper()}']
        return mag_norms

    def _process_df(self, df, band, data, component, host_type,
                    phot_params=None):
        mag_norms = self._get_mag_norms(component, host_type, band)
        for iloc in range(len(df)):
            row = df.iloc[iloc]
            ra = row['ra_lens']
            dec = row['dec_lens']
            unique_sys_id = str(row['unique_id'])
            if unique_sys_id not in mag_norms:
                continue
            redshift = row['redshift']
            gAv = row['av_mw']
            gRv = row['rv_mw']
            unique_id = f'{unique_sys_id}'# + component[0]
            sed_file = stc.find_sed_file(
                row[f'sed_{component}_host'].lstrip('b').strip("'"))
            mag_norm = mag_norms[unique_sys_id]
            iAv = row[f'av_internal_{component}']
            iRv = row[f'rv_internal_{component}']
            synth_phot = stc.SyntheticPhotometry(sed_file, mag_norm,
                                                 redshift=redshift,
                                                 iAv=iAv, iRv=iRv,
                                                 gAv=gAv, gRv=gRv)
            data['id'].append(unique_id)
            data['ra'].append(ra)
            data['dec'].append(dec)
            data['flux'].append(synth_phot.calcFlux(band))
            data['adu'].append(synth_phot.sed.calcADU(synth_phot.bp_dict[band],
                                                      phot_params))
            data['host_type'].append(host_type)
            data['component'].append(component)
        return data

if __name__ == '__main__':
    hosts_truth_cat = HostsTruthCat('../FITS_stamps')
    visit = 709680
    df = hosts_truth_cat(visit)
    df.to_pickle(f'host_fluxes_v{visit}.pkl')
