import pickle
import os
import time
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

opsim_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
opsim_file = os.path.join(opsim_dir, 'minion_1016_desc_dithered_v4.db')
assert os.path.isfile(opsim_file)

out_file = os.path.join(os.environ['SCRATCH'],
                        'minion_1016_desc_dithered_dict.p')

obs_gen = ObservationMetaDataGenerator(opsim_file)

t_start = time.time()
obs_md = obs_gen.getObservationMetaData(boundLength=2.1,
                                        boundType='circle',
                                        obsHistID=(-10,1000000000))

print('getting records took %e' % (time.time()-t_start))

out_dict = {}

t_start = time.time()
for obs in obs_md:
    out_dict[obs.OpsimMetaData['obsHistID']] = obs

with open(out_file, 'wb') as out_file:
    pickle.dump(out_dict, out_file)
print('output took %e' % (time.time()-t_start))
