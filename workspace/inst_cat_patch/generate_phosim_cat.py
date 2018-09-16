import os
from lsst.sims.catUtils.baseCatalogModels import StarObj
from desc.sims.GCRCatSimInterface import get_obs_md
from desc.sims.GCRCatSimInterface import make_instcat_header

def patch_dir(dir_name):

    if not hasattr(patch_dir, 'star_db'):
        star_db = StarObj(database='LSSTCATSIM',
                          host='fatboy.phys.washington.edu',
                          port=1433, driver='mssql+pymssql')

        patch_dir.star_db = star_db

    n_complete = 0
    list_of_files = os.listdir(dir_name)
    for file_name in list_of_files:
        if file_name == 'job_log.txt':
            continue
        if not file_name.startswith('job_log'):
            continue

        job_log = os.path.join(dir_name, file_name)
        is_complete = False
        with open(job_log, 'r') as in_file:
            job_lines = in_file.readlines()
        if 'all done' in job_lines[-1]:
            is_complete = True
            n_complete += 1

    print('n_complete %d' % n_complete)

if __name__ == "__main__":

    dir_name = '/global/cscratch1/sd/desc/DC2/Run2.0i/instCat/edison_packed_submit_idx200_size200'

    patch_dir(dir_name) 
