import os
import h5py
import numpy as np

dtype = np.dtype([('gid', int), ('hpid', int), ('sed', str, 100),
                  ('magnorm', float)])

baseline_disks = np.genfromtxt('baseline/disks.txt', dtype=dtype)

h5_name = os.path.join(os.environ['SCRATCH'], 'sed_cache',
                  'disk_10451.h5')

assert os.path.isfile(h5_name)

disks = h5py.File(h5_name, 'r')


dexes = np.in1d(disks['galaxy_id'].value, baseline_disks['gid'])
gid = disks['galaxy_id'].value[dexes]
np.testing.assert_array_equal(gid, baseline_disks['gid'])

sed = disks['sed_name'].value[dexes]
magnorm = disks['mag_norm'].value[dexes]

print(gid[1],baseline_disks['gid'][1])
print(magnorm[1],baseline_disks['magnorm'][1])
np.testing.assert_array_equal(magnorm, baseline_disks['magnorm'])

np.testing.assert_array_equal(sed.astype(str), baseline_disks['sed'])
