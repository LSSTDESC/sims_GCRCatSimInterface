import sys
import os
import h5py
import GCRCatalogs
import re
import numpy as np

catalog_filename_template='cosmoDC2_v1.0_{}'

zlos = [0, 1, 2] 
zhis = [z+1 for z in zlos]
healpixel = 9556

output_dir = '/global/homes/k/kovacs/CatSim_CosmoDC2_Match/CatSim_CosmoDC2_Data'
output_file_template = 'catsim_fit_solution_{}_z_{}_{}_healpix_{}.hdf5'

filedict = {'1':{'catsim_file_dir':'/global/homes/k/kovacs/CatSim_CosmoDC2_Match/CatSim_CosmoDC2_Data',
            'catsim_file_name':'extincted_fit_solution_1',
            'catsim_file_template':'{}_{}{}.hdf5',
            'frames':[''],
           },
            '2_lsst':{'catsim_file_dir':'/global/cscratch1/sd/danielsf/extincted_galaxy_lsst_fits',
                      'catsim_file_name':'galaxy_fit_lsst_a0.01',
                      'catsim_file_template':'{}_h{}{}.h5',
                      'frames':['', '_observer'],
           },
            '2_sed':{'catsim_file_dir':'/global/cscratch1/sd/danielsf/extincted_galaxy_fit/',
                     'catsim_file_name':'galaxy_fit',
                     'catsim_file_template':'{}_h{}_181008b{}.h5',
                     'frames':['_observer'],
           },
      }

def rewrite_file(solution='1'):
    gc=GCRCatalogs.load_catalog(catalog_filename_template.format(healpixel))
    catsim_files =[]
    location = filedict[solution]

    for f in location['frames']:
        fr = f if '2_sed' not in solution else ''  #fix up naming inconsistency
        catsim_file = os.path.join(location['catsim_file_dir'], 
                                   location['catsim_file_template'].format(location['catsim_file_name'], healpixel, fr))
        catsim_files.append(catsim_file)
        print('Reading CatSim-fit file {}'.format(catsim_file))

    data_block_start = 0
    for zlo, zhi in zip(zlos, zhis):
        data = gc.get_quantities(['redshift', 'galaxy_id', 'magnification'], native_filters=[(lambda z: z == zlo, 'redshift_block_lower')])
        datlen = len(data['galaxy_id'])
        out_file_name = output_file_template.format(solution, zlo, zhi, healpixel)
        print('Processing block {}-{} with {} elements and saving to {}'.format(zlo, zhi, datlen, out_file_name))
        #write data block to output file
        with h5py.File(os.path.join(output_dir, out_file_name), 'w') as out_file:
            cgroup = out_file.create_group('CatSim_fits')
            for f, catsim_file in zip(location['frames'], catsim_files):
                f_id = 'obs' if 'obs' in f else 'rest'
                with h5py.File(catsim_file, 'r') as catsim:
                    #check galaxy_id's match
                    assert np.array_equal(catsim['galaxy_id'].value[data_block_start:data_block_start+datlen], 
                                          data['galaxy_id'])
                    for k in list(catsim.keys()):
                        key = k.replace('cosmo_', 'cosmoDC2_'+f_id+'_')
                        if solution != '1':
                            key = key.replace('fit_', 'CatSim_fit_'+f_id+'_')
                        print(k, key)
                        if key not in cgroup.keys():
                            cgroup.create_dataset(key, data=catsim[k].value[data_block_start:data_block_start+datlen])

            cgroup.create_dataset('magnification', data=data['magnification'])
            data_block_start += datlen

            out_file.close()

    return    

                                        
