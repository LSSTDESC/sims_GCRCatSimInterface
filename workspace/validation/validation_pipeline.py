import os
import numpy as np

import desc.sims.GCRCatSimInterface.validation as validation

import time

import argparse
import subprocess
import shutil

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--file_name', type=str)
    parser.add_argument('--seed', type=int, default=4561)
    args = parser.parse_args()

    rng = np.random.RandomState(args.seed)

    project_dir = os.path.join('/global/projecta/projectdirs',
                               'lsst/groups/SSim/DC2/cosmoDC2_v1.1.4')

    agn_db = os.path.join(project_dir,
                          'agn_db_mbh7_mi30_sf4.db')

    already_run = set()

    scratch_dir = os.path.join(os.environ['SCRATCH'], 'test_validation')
    assert os.path.isdir(scratch_dir)

    orig_file = args.file_name
    p = subprocess.Popen(['cp',orig_file,scratch_dir])
    p.wait()
    instcat_tar_name = os.path.basename(orig_file)
    copied_instcat_name = os.path.join(scratch_dir, instcat_tar_name)
    untarred_instcat_name = os.path.join(scratch_dir,
                                         instcat_tar_name.replace('.tar.gz',''))
    print('tar name: ',instcat_tar_name)
    p = subprocess.Popen(['tar','-C',scratch_dir,
                          '-xf',os.path.join(scratch_dir, instcat_tar_name)])

    p.wait()
    print('done untarring')

    log_name = 'validated_catalogs.txt'

    while len(already_run)<1:
        obs_dir = instcat_tar_name.replace('.tar.gz','')
        obsid = int(obs_dir)
        if obsid in already_run:
            continue

        cat_dir = os.path.join(scratch_dir)

        t_start = time.time()
        #validate_sne(cat_dir, obsid, out_file=f_out)
        mag_seed = rng.randint(0,10000)
        validation.validate_instance_catalog_magnitudes(cat_dir, obsid,
                                                        seed=mag_seed,
                                                        nrows=10000)
        validation.validate_instance_catalog_positions(cat_dir, obsid, 2.1)
        validation.validate_agn_mags(cat_dir, obsid, agn_db)
        print('\n\nvalidated agn after %e seconds' % (time.time()-t_start))
        already_run.add(obsid)
        with open(log_name, 'a') as out_file:
            out_file.write('%d validated\n' % obsid)

        os.unlink(copied_instcat_name)
        shutil.rmtree(untarred_instcat_name)
