"""
This script will generate a lookup table that maps healpixel ID to obsHistID
in DC2
"""
import os
import numpy as np
import sqlite3
import h5py
import healpy

import time
import multiprocessing

from lsst.sims.utils import htmModule
from lsst.sims.utils import ra_dec_from_xyz
from lsst.sims.utils import xyz_from_ra_dec

from dc2_spatial_defintion import DC2_bounds

htmid_level = 6

def make_lookup(chunk, hs_list, my_lock, output_dict, mgr):
    """
    chunk will be the result of sqlite3.fetchmany() on the OpSim database;

    each row will be (obsHistID, descDitheredRA, descDitheredDec,
                      descDitheredRotTelPos, expMJD, filter)

    hs_list is a list of HalfSpaces defining DC2

    output_dict will need to be initialized with the valid values of htmid
    """
    print('looking up')
    bp_to_int = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}

    radius_deg = 2.1
    radius_rad = np.radians(radius_deg)
    local_lookup = {}
    local_dex = {}
    obsHistID_arr = []
    mjd_arr = []
    ra_arr = []
    dec_arr = []
    rotTelPos_arr = []
    filter_arr = []
    for pointing in chunk:
        obs_id = int(pointing[0])
        ra = float(pointing[1])
        dec = float(pointing[2])
        rotTelPos = float(pointing[3])
        mjd = float(pointing[4])
        bp = bp_to_int[pointing[5]]
        obs_xyz = xyz_from_ra_dec(np.degrees(ra), np.degrees(dec))

        obs_hs = htmModule.halfSpaceFromRaDec(np.degrees(ra),
                                              np.degrees(dec),
                                              radius_deg)

        is_valid_pointing = True
        for hs in hs_list:
            if not hs.intersects_circle(obs_hs.vector, obs_hs.phi):
                is_valid_pointing = False
                break

        if not is_valid_pointing:
            continue

        logged_metadata = False
        healpixel_list = healpy.query_disc(32, obs_xyz, radius_rad,
                                           nest=False, inclusive=True)
        for hpid in healpixel_list:
            if hpid not in local_lookup:
                local_lookup[hpid] = -1*np.ones(len(chunk), dtype=int)
                local_dex[hpid] = 0
            local_lookup[hpid][local_dex[hpid]] = obs_id
            local_dex[hpid] += 1

            if not logged_metadata:
                obsHistID_arr.append(obs_id)
                ra_arr.append(ra)
                dec_arr.append(dec)
                mjd_arr.append(mjd)
                rotTelPos_arr.append(rotTelPos)
                filter_arr.append(bp)
                logged_metadata = True

    obsHistID_arr = np.array(obsHistID_arr)
    ra_arr = np.array(ra_arr)
    dec_arr = np.array(dec_arr)
    mjd_arr = np.array(mjd_arr)
    rotTelPos_arr = np.array(rotTelPos_arr)
    filter_arr = np.array(filter_arr)

    with my_lock:
        for hpid in local_lookup:
            valid = np.where(local_lookup[hpid]>=0)
            if hpid not in output_dict:
                output_dict[hpid] = mgr.list()
            output_dict[hpid].append(local_lookup[hpid][valid])
            if len(output_dict['obsHistID'])>0:
                proto_obs = np.concatenate(output_dict['obsHistID'])
                valid = ~np.isin(obsHistID_arr, proto_obs)
            else:
                valid = np.array([True]*len(obsHistID_arr))
            output_dict['obsHistID'].append(obsHistID_arr[valid])
            output_dict['ra'].append(ra_arr[valid])
            output_dict['dec'].append(dec_arr[valid])
            output_dict['mjd'].append(mjd_arr[valid])
            output_dict['rotTelPos'].append(rotTelPos_arr[valid])
            output_dict['filter'].append(filter_arr[valid])
    print('leaving')

if __name__ == "__main__":

    data_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth'
    assert os.path.isdir(data_dir)
    opsim_name = os.path.join(data_dir, 'minion_1016_desc_dithered_v4_sfd.db')
    assert os.path.isfile(opsim_name)

    out_name = os.path.join(data_dir, 'hpid_to_obsHistID_lookup.h5')
    if os.path.isfile(out_name):
        raise RuntimeError("\n\n%s\nalready exists\n" % out_name)

    dc2b = DC2_bounds()
    hs_list = dc2b.hs_list

    mgr = multiprocessing.Manager()
    my_lock = mgr.Lock()
    out_dict = mgr.dict()
    out_dict['obsHistID'] = mgr.list()
    out_dict['mjd'] = mgr.list()
    out_dict['ra'] = mgr.list()
    out_dict['dec'] = mgr.list()
    out_dict['rotTelPos'] = mgr.list()
    out_dict['filter'] = mgr.list()

    n_threads = 40
    chunk_size = 100000
    query = "SELECT obsHistID, descDitheredRA, descDitheredDec, "
    query += "descDitheredRotTelPos, expMJD, filter "
    query += "FROM Summary GROUP BY obsHistID "

    t_start = time.time()
    p_list = []
    n_chunks = 0
    with sqlite3.connect(opsim_name) as connection:
        cursor = connection.cursor()
        iterator = cursor.execute(query)
        chunk = iterator.fetchmany(chunk_size)
        while len(chunk) > 0:
            p = multiprocessing.Process(
                      target=make_lookup,
                      args=(chunk, hs_list, my_lock, out_dict, mgr))
            p.start()
            p_list.append(p)
            chunk = iterator.fetchmany(chunk_size)
            while len(p_list)>=n_threads:
                exit_list = []
                for p in p_list:
                    exit_list.append(p.exitcode)
                was_popped = False
                for ii in range(len(p_list)-1,-1,-1):
                    if exit_list[ii] is not None:
                        p_list.pop(ii)
                        was_popped = True
                        n_chunks += 1
                if was_popped:
                    duration = (time.time()-t_start)/3600.0
                    per = duration/n_chunks
                    print(n_chunks,duration,per)

    for p in p_list:
        p.join()

    meta_key_list = ['obsHistID', 'ra', 'dec', 'mjd', 'rotTelPos', 'filter']

    with h5py.File(out_name, 'w') as out_file:
        for htmid in out_dict:
            if htmid in meta_key_list:
                continue
            data = np.sort(np.concatenate(out_dict[htmid])).astype(int)
            out_file.create_dataset('%d' % htmid, data=data)

        obs_id = np.concatenate(out_dict['obsHistID'])
        sorted_dex = np.argsort(obs_id)
        obs_id = obs_id[sorted_dex].astype(int)
        out_file.create_dataset('obsHistID', data=obs_id)
        out_file.create_dataset('filter',
               data=np.concatenate(out_dict['filter'])[sorted_dex].astype(int))

        for meta_key in ['mjd', 'ra', 'dec', 'rotTelPos']:
            out_file.create_dataset(meta_key,
                       data=np.concatenate(out_dict[meta_key])[sorted_dex])
