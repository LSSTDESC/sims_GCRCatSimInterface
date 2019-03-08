import os

out_name = 'list_of_completed_obshistid.txt'

cat_dir = os.path.join('/global/projecta/projectdirs/lsst/production',
                       'DC2_ImSim/Run2.1i/instCat')

assert os.path.isdir(cat_dir)

contents = os.listdir(cat_dir)

with open(out_name,'w') as out_file:
    for sub_dir in contents:
        sorted_dir = os.path.join(cat_dir, sub_dir)
        if not os.path.isdir(sorted_dir):
            continue
        sub_contents = os.listdir(sorted_dir)
        for name in sub_contents:
            full_name = os.path.join(sorted_dir, name)
            if os.path.isfile(full_name):
                obsid = name.replace('.tar.gz','')
            else:
                obsid = name

            print(obsid)
            try:
                obsid_int = int(obsid)
                out_file.write('%d\n' % obsid_int)
            except ValueError:
                pass
