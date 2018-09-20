import os
import time
import multiprocessing
from generate_phosim_cat import patch_dir

if __name__ == "__main__":

    opsim_db = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
    assert os.path.isdir(opsim_db)
    opsim_db = os.path.join(opsim_db, 'minion_1016_desc_dithered_v4.db')

    project_scratch = '/global/cscratch1/sd/desc/DC2/Run2.0i/instCat/'
    assert os.path.isdir(project_scratch)

    parent_dir = os.path.join(project_scratch,
                              'cori_haswell_submit_idx2200_size200')

    dirs_to_patch = []
    for bp in 'ugrizy':
        bp_dir = os.path.join(parent_dir, '%s-WFD' % bp)
        obs_list = os.listdir(bp_dir)
        for obs in obs_list:
            full_dir = os.path.join(bp_dir, obs, 'instCat')
            assert os.path.isdir(full_dir)
            dirs_to_patch.append(full_dir)

    p_list = []
    t_start = time.time()
    d_dir = len(dirs_to_patch)//23
    for i_dir in range(0, len(dirs_to_patch), d_dir):
        p = multiprocessing.Process(target=patch_dir,
                                    args=(dirs_to_patch[i_dir:i_dir+d_dir], opsim_db))
        p.start()
        p_list.append(p)
        if len(p_list)>=24:
            for p in p_list:
                p.join()
            p_list = []
            duration = (time.time()-t_start)/3600.0
            predicted = len(dirs_to_patch)*duration/i_dir
            print('%d in %.2e; predict %.2e' % (i_dir, duration, predicted))

    for p in p_list:
        p.join()

    print('all done')
