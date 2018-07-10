import os
import tempfile
import time
from desc.sims.GCRCatSimInterface import write_sprinkled_param_db
from desc.sims.GCRCatSimInterface import write_sprinkled_lc
from lsst.sims.utils import ObservationMetaData
from lsst.sims.photUtils import BandpassDict

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sql_dir', type=str, default=None)
parser.add_argument('--fov', type=float, default=None)
args = parser.parse_args()

if args.fov is None:
    raise RuntimeError("must specify fov")

bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

yaml_file = 'proto-dc2_v2.1.2'

agn_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
assert os.path.isdir(agn_dir)

agn_db = os.path.join(agn_dir, 'agn_db_mbh_7.0_m_i_30.0.sqlite')
assert os.path.isfile(agn_db)

obs = ObservationMetaData(pointingRA=55.064,
                          pointingDec=-29.783,
                          mjd=59580.0,
                          boundType='circle',
                          boundLength=args.fov)

opsim_db = '/global/projecta/projectdirs/lsst'
opsim_db += '/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'

pointing_dir = os.path.join(os.environ['HOME'], 'DC2_Repo', 'data', 'Run1.1')

h5_name = os.path.join(os.environ['SCRATCH'], 'h5_file_test', 'agn_lc.h5')
h5_name = os.path.join(os.environ['SCRATCH'], 'sql_lc_test_small', 'agn_lc.db')

if os.path.exists(h5_name):
    os.unlink(h5_name)

#sql_dir = tempfile.mkdtemp(dir=os.environ['SCRATCH'],
#                           prefix='sprinkled_sql_')


t_start = time.time()
(sql_file_name, table_names) = write_sprinkled_param_db(obs,
                                                 field_ra=55.064,
                                                 field_dec=-29.783,
                                                 agn_db=agn_db,
                                                 yaml_file=yaml_file,
                                                 out_dir=args.sql_dir,
                                                 bp_dict=bp_dict)

print('\nwrote\n%s\n' % sql_file_name)
print('fov %e in %e sec' % (args.fov, time.time()-t_start))
exit()

write_sprinkled_lc(h5_name, obs,
                   pointing_dir, opsim_db,
                   sql_file_name=sql_file_name,
                   bp_dict=bp_dict)

print('\nparam db\n%s\n' % sql_file_name)
