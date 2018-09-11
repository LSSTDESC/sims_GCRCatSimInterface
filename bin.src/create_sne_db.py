"""
This script will a directory containing csv files with the SNe parameters
for DC2 and load those parameters into a sqlite file that can be queried
by the InstanceCatalog generation code.
"""
import sqlite3
import os
from lsst.sims.utils import findHtmid
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--out_file', type=str, default=None,
                        help='name of the file to write')
    parser.add_argument('--in_dir', type=str, default=None,
                        help="name of directory containing the csv files "
                        "to be read into the database; all files ending in "
                        "'.csv' will be read.")

    args = parser.parse_args()
    if args.out_file is None:
        raise RuntimeError("Must specify out_file")
    if args.in_dir is None:
        raise RuntimeError("Must specify in_dir")

    data_dir_file_names = os.listdir(args.in_dir)

    if os.path.exists(args.out_file):
        raise RuntimeError("%s already exists" % args.out_file)

    htmid_level = 6

    with sqlite3.connect(args.out_file) as conn:
        cursor = conn.cursor()
        creation_cmd = """CREATE TABLE sne_params (
                          htmid_level_6 int,
                          galaxy_id int,
                          c_in real,
                          mB real,
                          t0_in real,
                          x0_in real,
                          x1_in real,
                          z_in real,
                          snid_in text,
                          snra_in real,
                          sndec_in real)"""

        cursor.execute(creation_cmd)
        conn.commit()
        for file_name in data_dir_file_names:
            if not file_name.endswith('csv'):
                continue
            full_name = os.path.join(args.in_dir, file_name)
            print('reading %s' % file_name)

            with open(full_name, 'r') as input_file:
                input_lines = input_file.readlines()[1:]

            params_list = []
            for line in input_lines:
                params = line.strip().split(',')
                assert len(params) == 23

                htmid = findHtmid(float(params[21]), float(params[22]),
                                  max_level=htmid_level)

                pv = [htmid, int(params[0])]
                pv += [float(pp) for pp in params[1:6]]
                pv += [float(params[11])]
                pv += [params[20]]
                pv += [float(params[21]), float(params[22])]
                assert len(pv) == 11
                params_list.append(tuple(pv))
            cursor.executemany('''INSERT INTO sne_params
                               VALUES(?,?,?,?,?,?,?,?,?,?,?)''', params_list)

            conn.commit()

        cursor.execute('''CREATE INDEX htmid_index ON sne_params (htmid_level_6)''')
        conn.commit()
