"""
This script will copy non-hdf5 files that the InstanceCatalog
needs to access into scratch space so that each InstanceCatalog
generation node can have access to its own supporting data (hopefully
limiting resource contention).

Call it as

python make_config.py $NEW_DIR_FOR_DATA name_of_config_file.json

The data will be copies into $NEW_DIR_FOR_DATA and a bespoke
config file will be written in name_of_config_file.json, referring
to the copied data
"""
import os
import shutil
import json
import time

import argparse


def copy_instcat_config(dest_dir, new_config_name):

    t_start = time.time()
    with open('config_file_2.1.static.json', 'rb') as in_file:
        new_config = json.load(in_file)

    with open('config_file_2.1.to_copy.json', 'rb') as in_file:
        to_copy_config = json.load(in_file)

    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    for name in to_copy_config:
        old_file = to_copy_config[name]
        if os.path.isfile(old_file):
            new_file = os.path.join(dest_dir, os.path.basename(old_file))
            shutil.copy(old_file, new_file)
        elif os.path.isdir(old_file):
            sub_dirs = old_file.split('/')
            sub_dirs.reverse()
            chosen_dir = None
            for sub in sub_dirs:
                if len(sub)>1:
                    chosen_dir = sub
                    break

            new_file = os.path.join(dest_dir, chosen_dir)
            shutil.copytree(old_file, new_file)

        new_config[name] = new_file

    with open(new_config_name, 'w') as out_file:
        json.dump(new_config, out_file)

    print('wrote %s after %e seconds' % (new_config_name, time.time()-t_start))

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', type=str, default=None)
    parser.add_argument('--file', type=str, default=None)
    args = parser.parse_args()

    assert args.dir is not None
    assert args.file is not None

    copy_instcat_config(args.dir, args.file)
