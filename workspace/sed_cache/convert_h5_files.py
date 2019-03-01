import json
import h5py
import os
import numpy as np
import time
import GCRCatalogs
from GCR import GCRQuery

data_dir = os.path.join(os.environ['SCRATCH'], 'sed_cache_181017')
assert os.path.isdir(data_dir)

galaxy_dir = os.path.join(os.environ['SIMS_SED_LIBRARY_DIR'], 'galaxySED')
assert os.path.isdir(galaxy_dir)
list_of_files = os.listdir(galaxy_dir)
list_of_files.sort()

sed_name_dict = {}
name_array = np.empty(len(list_of_files), dtype=(bytes, 100))
for i_nn, nn in enumerate(list_of_files):
    full_name = os.path.join('galaxySED', nn)
    sed_name_dict[full_name] = i_nn
    name_array[i_nn] = full_name


hp_list = [10069, 10196, 9940, 10068, 10195, 10197, 10323, 10324, 10325, 10447, 10448]

cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')

t_start = time.time()
for hp in hp_list:
    print('processing %d; %.1e minutes' % (hp, (time.time()-t_start)/60.0))
    in_file_name = os.path.join(data_dir, 'sed_fit_%d.h5' % hp)
    assert os.path.isfile(in_file_name)

    hp_query = GCRQuery('healpix_pixel==%d' % hp)
    cat_qties = cat.get_quantities(['galaxy_id', 'ra', 'dec'], native_filters=[hp_query])

    sorted_dex = np.argsort(cat_qties['galaxy_id'])
    for k in cat_qties.keys():
        cat_qties[k] = cat_qties[k][sorted_dex]

    with h5py.File(in_file_name, 'r') as in_file:
        disk_names = in_file['disk_sed'].value.astype(str)
        disk_name_idx = np.array([sed_name_dict[nn] for nn in disk_names])
        print('    got disks')
        bulge_names = in_file['bulge_sed'].value.astype(str)
        bulge_name_idx = np.array([sed_name_dict[nn] for nn in bulge_names])
        print('    got bulges')

        mapped_dex = np.searchsorted(cat_qties['galaxy_id'], in_file['galaxy_id'].value)
        np.testing.assert_array_equal(cat_qties['galaxy_id'][mapped_dex],
                                      in_file['galaxy_id'].value)

        out_file_name = os.path.join(data_dir, 'sed_fit_idx_%d.h5' % hp)
        with h5py.File(out_file_name, 'w') as out_file:
            for k in in_file.keys():
                if k == 'disk_sed' or k == 'bulge_sed':
                    continue
                out_file.create_dataset(k, data=in_file[k].value)
            d_sed = out_file.create_dataset('disk_sed', data=disk_name_idx)
            b_sed = out_file.create_dataset('bulge_sed', data=bulge_name_idx)
            out_file.create_dataset('sed_names', data=name_array)
            out_file.create_dataset('ra', data=cat_qties['ra'][mapped_dex])
            out_file.create_dataset('dec', data=cat_qties['dec'][mapped_dex])
