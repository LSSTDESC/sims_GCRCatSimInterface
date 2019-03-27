import os

out_file = 'list_of_completed_obshistid.txt'
instcat_dir = '/global/projecta/projectdirs/lsst/production/'
instcat_dir = os.path.join(instcat_dir, 'DC2_ImSim/Run2.1i/instCat')
assert os.path.isdir(instcat_dir)

super_contents = os.listdir(instcat_dir)
with open(out_file, 'w') as out_file:
    for sub_dir in super_contents:
        if 'to' not in sub_dir:
            continue
        full_name = os.path.join(instcat_dir, sub_dir)
        if not os.path.isdir(full_name):
            continue
        sub_contents = os.listdir(full_name)
        for name in sub_contents:
            if not name.endswith('.tar.gz'):
                continue
            obsid = int(name.replace('.tar.gz',''))
            out_file.write('%d\n' % obsid)
