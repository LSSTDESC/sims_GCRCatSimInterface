import os

param_dir = os.path.join('/astro', 'store', 'pogo3', 'danielsf')
param_dir = os.path.join(param_dir, 'truth_181002')

assert os.path.isdir(param_dir)

param_file = os.path.join(param_dir, 'truth_params_181002_protodc2_v3.db')

assert os.path.isfile(param_file)

ptng_dir = os.path.join('/local', 'lsst', 'danielsf', 'DC2-production', 'data')
ptng_dir = os.path.join(ptng_dir, 'Run1.1')

assert os.path.isdir(ptng_dir)

opsim_db = os.path.join('/local', 'lsst', 'danielsf', 'OpSimData')
opsim_db = os.path.join(opsim_db, 'minion_1016_desc_dithered_v4.db')

assert os.path.isfile(opsim_db)

from lsst.sims.utils import ObservationMetaData
obs_tot = ObservationMetaData(pointingRA=55.064, pointingDec=-29.783,
                              boundType='circle', boundLength=4.0)

out_name = os.path.join(param_dir, 'run_1.2_trial_lc.db')

assert not os.path.isfile(out_name)

from lsst.sims.photUtils import BandpassDict

bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

from desc.sims.GCRCatSimInterface import write_sprinkled_lc

write_sprinkled_lc(out_name, obs_tot, ptng_dir, opsim_db,
                   sql_file_name=param_file,
                   bp_dict=bp_dict)
