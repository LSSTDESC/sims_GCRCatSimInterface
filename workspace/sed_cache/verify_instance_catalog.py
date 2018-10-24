import os
import numpy as np
import pandas as pd
import healpy

import GCRCatalogs
from GCR import GCRQuery

from lsst.sims.utils import angularSeparation

colnames = ['obj', 'uniqueID', 'ra', 'dec', 'magnorm', 'sed', 'redshift', 'g1', 'g2',
            'kappa', 'dra', 'ddec', 'src_type', 'major', 'minor',
            'positionAngle', 'sindex', 'dust_rest', 'rest_av', 'rest_rv',
            'dust_obs', 'obs_av', 'obs_rv']

col_types = {'magnorm': float, 'redshift': float,
             'rest_av': float, 'rest_rv': float,
             'sed': bytes, 'uniqueID': int}

data_dir = os.path.join(os.environ['SCRATCH'], 'instcat_181024_verify', '00277065')
assert os.path.isdir(data_dir)

disk_file = os.path.join(data_dir, 'disk_gal_cat_277065.txt.gz')
assert os.path.isfile(disk_file)

bulge_file = os.path.join(data_dir, 'bulge_gal_cat_277065.txt.gz')
assert os.path.isfile(bulge_file)

knots_file = os.path.join(data_dir, 'knots_cat_277065.txt.gz')
assert os.path.isfile(knots_file)

disk_df = pd.read_csv(disk_file, delimiter=' ',
                        compression='gzip', names=colnames, dtype=col_types)
disk_df['galaxy_id'] = pd.Series(disk_df['uniqueID']//1024, index=disk_df.index)
disk_df = disk_df.set_index('galaxy_id')

bulge_df = pd.read_csv(bulge_file, delimiter=' ',
                       compression='gzip', names=colnames, dtype=col_types)
bulge_df['galaxy_id'] = pd.Series(bulge_df['uniqueID']//1024, index=bulge_df.index)
bulge_df = bulge_df.set_index('galaxy_id')

for ii in range(len(colnames)):
    colnames[ii] = colnames[ii]+'_knots'

knots_df = pd.read_csv(knots_file, delimiter=' ',
                       compression='gzip', names=colnames, dtype=col_types)
knots_df['galaxy_id'] = pd.Series(knots_df['uniqueID_knots']//1024, index=knots_df.index)
knots_df = knots_df.set_index('galaxy_id')



#galaxy_df = disk_df.join(bulge_df, how='outer', lsuffix='_disk', rsuffix='_bulge')
wanted_col = ['sed', 'magnorm', 'redshift', 'rest_av', 'rest_rv', 'ra', 'dec']
galaxy_df = disk_df[wanted_col].join(bulge_df[wanted_col], how='outer', lsuffix='_disk', rsuffix='_bulge')
for ii in range(len(wanted_col)):
    wanted_col[ii] = wanted_col[ii]+'_knots'
galaxy_df = galaxy_df.join(knots_df[wanted_col], how='outer', rsuffix='_knots')

galaxy_df = galaxy_df.sort_index()

ra_center = np.nanmedian(galaxy_df['ra_disk'].values)
dec_center = np.nanmedian(galaxy_df['dec_disk'].values)

dd = angularSeparation(ra_center, dec_center, galaxy_df['ra_disk'].values, galaxy_df['dec_disk'].values)
radius_deg = np.nanmax(dd)
ra_rad = np.radians(ra_center)
dec_rad = np.radians(dec_center)
vv = np.array([np.cos(ra_rad)*np.cos(dec_rad),
               np.sin(ra_rad)*np.cos(dec_rad),
               np.sin(dec_rad)])

healpix_list = healpy.query_disc(32, vv, np.radians(radius_deg),
                                 nest=False, inclusive=True)

print('healpix list')
print(healpix_list)
print(ra_center, dec_center)

hp_query = GCRQuery()
for hp in healpix_list:
    hp_query |= GCRQuery('healpix_pixel==%d' % hp)

print('built final df')
cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')
gid = cat.get_quantities('galaxy_id', native_filters=[hp_query])['galaxy_id']
print('loaded galaxy_id')
valid_dexes = np.where(np.in1d(gid, galaxy_df.index.values))
print('got valid_dexes')
print(len(valid_dexes[0]))
print(len(galaxy_df))
