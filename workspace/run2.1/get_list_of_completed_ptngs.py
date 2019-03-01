import shutil
import os

generated_obshistid = []

super_dir = os.path.join('/global/cscratch1/sd/desc/DC2/Run2.1i/instCat')
assert os.path.isdir(super_dir)

sub_dir_list = os.listdir(super_dir)
for sub_dir in sub_dir_list:
    if not os.path.isdir(os.path.join(super_dir, sub_dir)):
        continue

    obs_dir = os.listdir(os.path.join(super_dir, sub_dir))
    for obs in obs_dir:
        try:
            obsid = int(obs)
            generated_obshistid.append(obsid)
        except ValueError:
            pass

other_dir = os.path.join('/global/projecta/projectdirs/lsst/groups',
                         'SSim/DC2/Run2.1i/instCat')

other_dir_list = [other_dir]
other_dir_list.append(os.path.join('/global/projecta/projectdirs/lsst/groups',
                                   'SSim/DC2/Run2.1i/instCat_0'))

wld_keep = 0
wld_delete = 0
for other_dir in other_dir_list:
    job_log_list = os.listdir(other_dir)
    for job_log in job_log_list:
        if not job_log.startswith('job'):
            continue

        p = job_log.replace('.txt','').split('_')[-1]
        obsid = int(p)
        full_name = os.path.join(other_dir, job_log)
        is_complete = False
        with open(full_name, 'r') as in_file:
            in_lines = in_file.readlines()
        if 'all done' in in_lines[-1]:
            is_complete = True

        if is_complete:
            wld_keep += 1
            generated_obshistid.append(obsid)
        else:
            wld_delete += 1
            target_dir = os.path.join(other_dir,'%.8d' % obsid)
            print('deleting %d' % obsid)
            shutil.rmtree(target_dir)
            os.unlink(full_name)

print('keep ',wld_keep)
print('delete ',wld_delete)

with open('list_of_completed_obshistid.txt', 'w') as out_file:
    for obs in generated_obshistid:
        out_file.write('%d\n' % obs)
