import os
import numpy as np
import pandas as pd

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
wanted_col = ['sed', 'magnorm', 'redshift', 'rest_av', 'rest_rv']
galaxy_df = disk_df[wanted_col].join(bulge_df[wanted_col], how='outer', lsuffix='_disk', rsuffix='_bulge')
for ii in range(len(wanted_col)):
    wanted_col[ii] = wanted_col[ii]+'_knots'
galaxy_df = galaxy_df.join(knots_df[wanted_col], how='outer', rsuffix='_knots')
print(galaxy_df)
print(len(galaxy_df))
