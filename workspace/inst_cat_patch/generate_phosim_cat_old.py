import os
import multiprocessing
from generate_phosim_cat import patch_dir

if __name__ == "__main__":

    opsim_db = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
    assert os.path.isdir(opsim_db)
    opsim_db = os.path.join(opsim_db, 'minion_1016_desc_dithered_v4.db')

    project_scratch = '/global/cscratch1/sd/desc/DC2/Run2.0i/instCat/'
    assert os.path.isdir(project_scratch)

    parent_dir = os.path.join(project_scratch,
                              'cori_haswell_submit_dix2200_size200')

    dirs_to_patch = []
    for bp in 'ugrizy':
        bp_dir = os.path.join(parent_dir, '%s-WFD')
        obs_list = os.listdir(bp_dir)
        for obs in obs_list:
            full_dir = os.path.join(bp_dir, obs, 'instCat')
            assert os.path.isdir(full_dir)
            dirs_to_patch.append(full_dir)

    p_list = []
    for dir_name in dirs_to_patch:
        p = multiprocessing.Process(target=patch_dir,
                                    args=(dir_name, opsim_db))
        p.start()
        p_list.append(p)
        if len(p_list)>=24:
            for p in p_list:
                p.join()
            p_list = []
    for p in p_list:
        p.join()
