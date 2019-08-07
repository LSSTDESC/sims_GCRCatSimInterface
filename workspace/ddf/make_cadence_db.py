import numpy as np
import os
import sqlite3
import lsst.sims.utils as sims_utils

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--out_file', type=str, default=None)
    parser.add_argument('--in_file', type=str, default=None)
    args = parser.parse_args()

    if args.out_file is None:
        raise RuntimeError("Must specify out_file")

    if args.in_file is None:
        raise RuntimeError("Must specify in_file")

    cadence_file_name = args.out_file
    if os.path.exists(cadence_file_name):
        os.unlink(cadence_file_name)

    dtype = np.dtype([('id', int), ('ra', float), ('dec', float),
                      ('mjd', float), ('rawSeeing', float),
                      ('filter', str, 1), ('rotSkyPos', float)])

    data_file = args.in_file
    data = np.genfromtxt(data_file, dtype=dtype, delimiter=',',
                         skip_header=1)

    rotTelPos = []
    for ii in range(len(data)):
        obs = sims_utils.ObservationMetaData(pointingRA=np.degrees(data['ra'][ii]),
                                             pointingDec=np.degrees(data['dec'][ii]),
                                             mjd=data['mjd'][ii])
        rottel = sims_utils._getRotTelPos(data['ra'][ii], data['dec'][ii],
                                          obs, data['rotSkyPos'][ii])

        rotTelPos.append(rottel)

    with sqlite3.connect(cadence_file_name) as conn:
        cur = conn.cursor()
        cur.execute('''CREATE TABLE Summary(obsHistID int,
                       descDitheredRA real, descDitheredDec real,
                       expMJD real, rawSeeing real, filter text,
                       rotSkyPos real, descDitheredRotTelPos real,
                       fieldRA real, fieldDec real)''')
        conn.commit()
        values = ((int(data['id'][ii]),
                   data['ra'][ii], data['dec'][ii],
                   data['mjd'][ii], data['rawSeeing'][ii],
                   data['filter'][ii], data['rotSkyPos'][ii],
                   rotTelPos[ii],
                   data['ra'][ii], data['dec'][ii]) for ii in range(len(data)))

        cur.executemany('''INSERT INTO Summary VALUES(?,?,?,?,?,?,?,?,?,?)''',
                        values)

        conn.commit()
        cur.execute("CREATE INDEX obs_id ON Summary (obsHistID)")
        conn.commit()
