"""
This script makes sure that the metadata columns were correctly
transcribed to htmid_6_to_obsHistID_lookup.h5
"""
import numpy as np
import os
import sqlite3
import h5py

data_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth'
opsim_name = os.path.join(data_dir, 'minion_1016_desc_dithered_v4_sfd.db')
assert os.path.isfile(opsim_name)
h5_name = os.path.join(data_dir, 'htmid_6_to_obsHistID_lookup.h5')
assert os.path.isfile(h5_name)

with h5py.File(h5_name, 'r') as lookup_file:
    obsHistID = lookup_file['obsHistID'][()]
    ra = lookup_file['ra'][()]
    dec = lookup_file['dec'][()]
    rotTelPos = lookup_file['rotTelPos'][()]
    mjd = lookup_file['mjd'][()]
    filter_arr = lookup_file['filter'][()]

sorted_dex = np.argsort(obsHistID)
obsHistID = obsHistID[sorted_dex]
ra = ra[sorted_dex]
dec = dec[sorted_dex]
rotTelPos = rotTelPos[sorted_dex]
mjd = mjd[sorted_dex]

bp_to_int = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}

d_max = -1.0
with sqlite3.connect(opsim_name) as conn:
    cursor = conn.cursor()
    query = "SELECT obsHistID, descDitheredRA, descDitheredDec, "
    query += "descDitheredRotTelPos, expMJD, filter "
    query += "FROM Summary GROUP BY obsHistID"

    iterator = cursor.execute(query)
    chunk = iterator.fetchmany(100000)
    while len(chunk)>0:
        obs = np.array([int(c[0]) for c in chunk])
        r = np.array([float(c[1]) for c in chunk])
        d = np.array([float(c[2]) for c in chunk])
        t = np.array([float(c[3]) for c in chunk])
        m  = np.array([float(c[4]) for c in chunk])
        f = np.array([bp_to_int[c[5]] for c in chunk])
        chunk = iterator.fetchmany(100000)

        valid = np.isin(obs, obsHistID)
        if valid.sum()==0:
           continue
        obs = obs[valid]
        r = r[valid]
        d = d[valid]
        t = t[valid]
        m = m[valid]
        f = f[valid]

        valid = np.searchsorted(obsHistID, obs)
        np.testing.assert_array_equal(obsHistID[valid], obs)
        for test, control in zip([r,d,t,m, f],
                                 [ra,dec,rotTelPos,mjd,filter_arr]):
            delta = np.abs((test-control[valid])/(control[valid]+1.0e-6)).max()
            if delta>d_max:
                d_max = delta
                print('d_max %e' % (d_max))
