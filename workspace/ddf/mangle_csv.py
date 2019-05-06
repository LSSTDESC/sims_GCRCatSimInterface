"""
This script will take the dc2_block*_v2.csv files with DDF observations
and mangle them into the same schema as def_cadence.csv, which is the
schema expected by make_cadence_files
"""
import numpy as np
import os

cadence_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/DDF_cadences'
assert os.path.isdir(cadence_dir)

file_name_list = ['dc2_blocks1_v2.csv', 'dc2_blocks2_v2.csv']

out_file_name = os.path.join(cadence_dir, 'dc2_block_cadences.csv')
if os.path.isfile(out_file_name):
    raise RuntimeError("\n\n%s\n\nalready exists\n" % out_file_name)

rng = np.random.RandomState(815)

id_used = set()
with open(out_file_name, 'w') as out_file:
    out_file.write('obsHistId,ra,dec,mjd,rawSeeing,filter,rotSkyPos\n')
    for file_name in file_name_list:
        full_name = os.path.join(cadence_dir, file_name)
        assert os.path.isfile(full_name)
        with open(full_name,'r') as in_file:
            for line in in_file:
                if line.startswith('exp'):
                    continue
                p = line.strip().split(',')
                ii = int(p[-1])
                assert ii not in id_used
                id_used.add(ii)
                rot = rng.random_sample()*2.0*np.pi
                out_file.write('%s,0.92721,-0.49044,%s,%s,%s,%.5f\n' %
                (p[4],p[0],p[3],p[1],rot))
