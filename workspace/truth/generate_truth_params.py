"""
This is the script that I used to extract parameters on all
of the extragalactic sources from protoDC2 and store them
in a sqlite database.

It took about 5 hours to run on protoDC2 running on Cori's
login node.

To scale it up for cosmoDC2, we will need to parallelize it
and submit it to the slurm queue.
"""

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
parser.add_argument('--sql_file', type=str, default='sprinkled_objects.sqlite')
parser.add_argument('--fov', type=float, default=None)
parser.add_argument('--yaml', type=str, default='proto-dc2_v2.1.2')
parser.add_argument('--agn', type=str,
                    default=os.path.join('/global/projecta/projectdirs',
                                         'lsst/groups/SSim/DC2',
                                         'agn_db_mbh_7.0_m_i_30.0.sqlite'))
args = parser.parse_args()

if args.fov is None:
    raise RuntimeError("must specify fov")

bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

assert os.path.isfile(args.agn)

obs = ObservationMetaData(pointingRA=55.064,
                          pointingDec=-29.783,
                          mjd=59580.0,
                          boundType='circle',
                          boundLength=args.fov)

opsim_db = '/global/projecta/projectdirs/lsst'
opsim_db += '/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'

pointing_dir = os.path.join(os.environ['HOME'], 'DC2_Repo', 'data', 'Run1.1')

t_start = time.time()
(sql_file_name, table_names) = write_sprinkled_param_db(obs,
                                                 field_ra=55.064,
                                                 field_dec=-29.783,
                                                 agn_db=args.agn,
                                                 yaml_file=args.yaml,
                                                 out_dir=args.sql_dir,
                                                 out_file=args.sql_file)

print('\nwrote\n%s\n' % sql_file_name)
print('fov %e in %e sec' % (args.fov, time.time()-t_start))
