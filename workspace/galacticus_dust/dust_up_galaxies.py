import os
import numpy as np
import gzip
import time

import multiprocessing

from lsst.utils import getPackageDir
from extinction_model import extinction_curve

def add_dust_to_galaxy(sed_name, input_dir, out_dir, a_grid):
    dtype = np.dtype([('wav', float), ('fl', float)])
    full_name = os.path.join(os.path.join(input_dir, sed_name))
    data = np.genfromtxt(full_name, dtype=dtype)
    wav_ang = 10.0*data['wav']
    for aa in a_grid:
        out_name = sed_name.replace('.spec','_a%.2f.spec' % aa)
        dust = a_grid[aa]
        with gzip.open(os.path.join(output_dir, out_name), 'w') as out_file:
            for ww, dd, ff in zip(data['wav'], dust, data['fl']):
                line='%.4f %.6e\n' % (ww, dd*ff)
                out_file.write(line.encode('utf-8'))


input_dir = os.path.join(getPackageDir('sims_sed_library'), 'galaxySED')
assert os.path.isdir(input_dir)

output_dir = os.path.join(os.environ['SCRATCH'], 'extincted_galaxy_seds')
assert os.path.isdir(output_dir)

a_grid = np.arange(0.01, 1.01, 0.01)
list_of_sed_files = os.listdir(input_dir)

a_grid_dict = {}

dtype = np.dtype([('wav', float), ('f', float)])
wav = None
for fname in list_of_sed_files:
    full_name = os.path.join(input_dir, fname)
    data = np.genfromtxt(full_name, dtype=dtype)
    if wav is not None:
        assert np.abs(wav-data['wav']).max()<1.0e-6
    else:
        wav = np.copy(data['wav'])
        wav_ang=10.0*wav
        for aa in a_grid:
            dust = extinction_curve(wav, aa)
            a_grid_dict[aa] = dust

print('starting fit')
p_list = []
t_start = time.time()
for i_file, sed_name in enumerate(list_of_sed_files):
    p = multiprocessing.Process(target=add_dust_to_galaxy,
                                args=(sed_name, input_dir,
                                      output_dir, a_grid_dict))
    p.start()
    p_list.append(p)
    if len(p_list)>=63 or i_file==(len(list_of_sed_files)-1):
        for p in p_list:
            p.join()
        p_list = []
        duration = (time.time()-t_start)/3600.0
        pred = len(list_of_sed_files)*duration/i_file
        print('run %d of %d; %.2e hrs left' %
        (i_file, len(list_of_sed_files), pred-duration))
