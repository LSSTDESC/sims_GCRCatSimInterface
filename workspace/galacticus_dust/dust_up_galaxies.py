import os
import numpy as np
import gzip
import time

from lsst.utils import getPackageDir
from extinction_model import extinction_curve

input_dir = os.path.join(getPackageDir('sims_sed_library'), 'galaxySED')
assert os.path.isdir(input_dir)

output_dir = os.path.join(os.environ['SCRATCH'], 'extincted_galaxy_seds')
assert os.path.isdir(output_dir)

a_grid = np.arange(0.01, 1.01, 0.01)
dtype = np.dtype([('wav', float), ('fl', float)])
list_of_sed_files = os.listdir(input_dir)

t_start = time.time()
for i_file, sed_name in enumerate(list_of_sed_files):
    if i_file>0 and i_file%10==0:
        duration = (time.time()-t_start)/3600.0
        pred = len(list_of_sed_files)*duration/i_file
        print('run %d of %d; %.2e hrs left' %
        (i_file, len(list_of_sed_files), pred))

    full_name = os.path.join(os.path.join(input_dir, sed_name))
    data = np.genfromtxt(full_name, dtype=dtype)
    wav_ang = 10.0*data['wav']
    for aa in a_grid:
        out_name = sed_name.replace('.spec','_a%.2f.spec' % aa)
        dust = extinction_curve(wav_ang, aa)
        with gzip.open(os.path.join(output_dir, out_name), 'w') as out_file:
            for ww, dd, ff in zip(data['wav'], dust, data['fl']):
                line='%.4f %.6e\n' % (ww, dd*ff)
                out_file.write(line.encode('utf-32'))
