import sqlite3
import os
from lsst.sims.utils import findHtmid


if __name__ == "__main__":

    db_name = 'dc2_sne_params.db'

    data_dir = 'data'
    data_dir_file_names = os.listdir(data_dir)

    if os.path.exists(db_name):
        raise RuntimeError("%s already exists" % db_name)

    htmid_level = 6

    with sqlite3.connect(db_name) as conn:
        cursor = conn.cursor()
        creation_cmd = """CREATE TABLE sne_params (
                          htmid_level_6 int,
                          galaxy_id int,
                          c_in real,
                          mB real,
                          t0_in real,
                          x0_in real,
                          x1_in real,
                          rand_host real,
                          zbin int,
                          decJ2000_gal real,
                          'morphology/positionAngle' real,
                          raJ2000_gal real,
                          z_in real,
                          'morphology/spheroidHalfLightRadiusArcsec' real,
                          'morphology/diskHalfLightRadiusArcsec' real,
                          'morphology/spheroidMinorAxisArcsec' real,
                          'morphology/diskMinorAxisArcsec' real,
                          totalMassStellar real,
                          stellar_mass_bulge real,
                          diskMassStellar real,
                          zbin_gals int,
                          snid text,
                          snra_in real,
                          sndec_in real)"""

        cursor.execute(creation_cmd)
        conn.commit()
        for file_name in data_dir_file_names:
            if not file_name.endswith('csv'):
                continue
            full_name = os.path.join(data_dir, file_name)
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
                pv += [float(pp) for pp in params[1:7]]
                pv += [int(params[7])]
                pv += [float(pp) for pp in params[8:19]]
                pv += [int(params[19])]
                pv += [params[20]]
                pv += [float(params[21]), float(params[22])]
                params_list.append(tuple(pv))
            cursor.executemany('''INSERT INTO sne_params
                               VALUES(?,?,?,?,?,?,?,?,?,?,?,
                               ?,?,?,?,?,?,?,?,?,?,?,?,?)''', params_list)

            conn.commit()

        cursor.execute('''CREATE INDEX htmid_index ON sne_params (htmid_level_6)''')
        conn.commit()
