import os
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.baseCatalogModels import StarObj
from desc.sims.GCRCatSimInterface import get_obs_md
from desc.sims.GCRCatSimInterface import make_instcat_header

import multiprocessing

import time

def patch_dir(dir_name=None, opsim_db=None):

    print('running on %s' % dir_name)

    if not hasattr(patch_dir, 'star_db'):
        star_db = StarObj(database='LSSTCATSIM',
                          host='fatboy.phys.washington.edu',
                          port=1433, driver='mssql+pymssql')

        patch_dir.star_db = star_db

        obs_gen = ObservationMetaDataGenerator(opsim_db)
        patch_dir.obs_gen = obs_gen

    t_start = time.time()

    n_complete = 0
    list_of_files = os.listdir(dir_name)
    n_files = len(list_of_files)
    print('%d in %s' % (n_files, dir_name))
    for i_file, file_name in enumerate(list_of_files):

        duration = (time.time()-t_start)
        predicted = n_files*duration/(i_file+1)
        if i_file>0 and i_file%20==0:
            print('ran %d of %d in %s; pred %e' %
            (i_file, n_files,dir_name,predicted))

        if file_name == 'job_log.txt':
            continue
        if not file_name.startswith('job_log'):
            continue

        job_log = os.path.join(dir_name, file_name)
        with open(job_log, 'r') as in_file:
            job_lines = in_file.readlines()
        if 'all done' in job_lines[-1]:
            n_complete += 1
            continue

        completed_catalogs = []
        for line in job_lines:
            if 'wrote star' in line:
                completed_catalogs.append('star')
            if 'wrote knots' in line:
                completed_catalogs.append('knots')
            if 'wrote bulge' in line:
                completed_catalogs.append('bulge_gal')
            if 'wrote disk' in line:
                completed_catalogs.append('disk_gal')
            if 'wrote agn' in line:
                completed_catalogs.append('agn_gal')
            if 'wrote SNe' in line:
                completed_catalogs.append('sne')

        params = os.path.basename(job_log).strip().replace('.txt','').split('_')
        if len(params)<3:
            continue

        out_dir = os.path.join(dir_name, params[2])
        if not os.path.isdir(out_dir):
            raise RuntimeError("%s is not a dir\n %s" % (out_dir,job_log))

        obshistid = int(params[2])

        out_name = os.path.join(out_dir,'phosim_cat_%d.txt' % obshistid)
        if os.path.exists(out_name):
            raise RuntimeError("%s already exists" % out_name)

        obs_md = get_obs_md(patch_dir.obs_gen,
                            obshistid, fov=2.1,
                            dither=True)

        sub_catalogs = os.listdir(out_dir)
        object_catalogs = []
        for sub_name in completed_catalogs:
            file_name ='%s_cat_%d.txt' % (sub_name, obshistid)

            if not (os.path.isfile(os.path.join(out_dir, file_name)) or
            os.path.isfile(os.path.join(out_dir, file_name+'.gz'))):
                raise RuntimeError("%s doesn't exist" %
                                   os.path.join(out_dir, file_name))

            object_catalogs.append(file_name)

        assert len(object_catalogs) > 0

        make_instcat_header(star_db, obs_md, out_name,
                            object_catalogs=object_catalogs)

    print('n_complete %d' % n_complete)
    print('took %e hours' % ((time.time()-t_start)/3600.0))

if __name__ == "__main__":

    opsim_db = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
    assert os.path.isdir(opsim_db)
    opsim_db = os.path.join(opsim_db, 'minion_1016_desc_dithered_v4.db')

    assert os.path.isfile(opsim_db)

    project_scratch = '/global/cscratch1/sd/desc/DC2/Run2.0i/instCat/'
    assert os.path.isdir(project_scratch)

    subdir_list = ['180914',
                   'edison_packed_submit_idx400_size200',
                   'edison_packed_submit_idx600_size200',
                   'edison_packed_submit_idx800_size200']

    dir_list = []
    for subdir in subdir_list:
        full_dir = os.path.join(project_scratch, subdir)
        if not os.path.isdir(full_dir):
            raise RuntimeError("%s is not a dir" % full_dir)
        dir_list.append(full_dir)

    #patch_dir(dir_name, opsim_db)
    j_list = []
    for dir_name in dir_list:
        p = multiprocessing.Process(target=patch_dir,
                                    kwargs={'dir_name':dir_name,
                                            'opsim_db':opsim_db})
        p.start()
        j_list.append(p)

    for p in j_list:
        p.join()

    print("all done")
