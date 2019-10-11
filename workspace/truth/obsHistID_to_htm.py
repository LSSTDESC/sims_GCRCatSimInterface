"""
This script will generate a lookup table that maps htmid_6 to obsHistID
in DC2
"""
import os
import numpy as np
import sqlite3
import h5py
from lsst.sims.utils import htmModule

import time
import multiprocessing

from lsst.sims.utils import ra_dec_from_xyz

htmid_level = 6

def make_lookup(chunk, hs_list, my_lock, output_dict):
    """
    chunk will be the result of sqlite3.fetchmany() on the OpSim database;

    each row will be (obsHistID, descDitheredRA, descDitheredDec,
                      descDitheredRotTelPos, expMJD)

    hs_list is a list of HalfSpaces defining DC2

    output_dict will need to be initialized with the valid values of htmid
    """
    radius = 2.1  # in degrees
    local_lookup = {}
    local_dex = {}
    obsHistID_arr = []
    mjd_arr = []
    ra_arr = []
    dec_arr = []
    rotTelPos_arr = []
    for pointing in chunk:
        obs_id = int(pointing[0])
        ra = float(pointing[1])
        dec = float(pointing[2])
        rotTelPos = float(pointing[3])
        mjd = float(pointing[4])
        obs_hs = htmModule.halfSpaceFromRaDec(np.degrees(ra),
                                              np.degrees(dec),
                                              radius)

        is_valid_pointing = True
        for hs in hs_list:
            if not hs.intersects_circle(obs_hs.vector, obs_hs.phi):
                is_valid_pointing = False
                break

        if not is_valid_pointing:
            continue

        logged_metadata = False
        obs_bounds = obs_hs.findAllTrixels(htmid_level)
        for bound in obs_bounds:
            for htmid in range(bound[0], bound[1]+1, 1):
                if htmid not in local_lookup:
                    local_lookup[htmid] = -1*np.ones(len(chunk), dtype=int)
                    local_dex[htmid] = 0
                local_lookup[htmid][local_dex[htmid]] = obs_id
                local_dex[htmid] += 1

                if not logged_metadata:
                    obsHistID_arr.append(obs_id)
                    ra_arr.append(ra)
                    dec_arr.append(dec)
                    mjd_arr.append(mjd)
                    rotTelPos_arr.append(rotTelPos)
                    logged_metadata = True

    obsHistID_arr = np.array(obsHistID_arr)
    ra_arr = np.array(ra_arr)
    dec_arr = np.array(dec_arr)
    mjd_arr = np.array(mjd_arr)
    rotTelPos_arr = np.array(rotTelPos_arr)

    with my_lock:
        for htmid in local_lookup:
            if htmid in output_dict:
                valid = np.where(local_lookup[htmid]>=0)
                output_dict[htmid].append(local_lookup[htmid][valid])
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

if __name__ == "__main__":

    data_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth'
    assert os.path.isdir(data_dir)
    opsim_name = os.path.join(data_dir, 'minion_1016_desc_dithered_v4_sfd.db')
    assert os.path.isfile(opsim_name)

    out_name = os.path.join(data_dir, 'htmid_6_to_obsHistID_lookup.h5')
    if os.path.isfile(out_name):
        raise RuntimeError("\n\n%s\nalready exists\n" % out_name)

    ne_corner = (71.46, -27.25)
    nw_corner = (52.25, -27.25)
    se_corner = (73.79, -44.33)
    sw_corner = (49.92, -44.33)
    pt_inside = (60.0, -30.0)

    hs_list = []
    hs_list.append(htmModule.halfSpaceFromPoints(ne_corner,
                                                 nw_corner,
                                                 pt_inside))

    hs_list.append(htmModule.halfSpaceFromPoints(ne_corner,
                                                 se_corner,
                                                 pt_inside))

    hs_list.append(htmModule.halfSpaceFromPoints(se_corner,
                                                 sw_corner,
                                                 pt_inside))

    hs_list.append(htmModule.halfSpaceFromPoints(sw_corner,
                                                 nw_corner,
                                                 pt_inside))

    dc2_bound = None
    for hs in hs_list:
        local_bound = hs.findAllTrixels(htmid_level)
        if dc2_bound is None:
            dc2_bound = local_bound
        else:
            dc2_bound = htmModule.HalfSpace.join_trixel_bound_sets(dc2_bound,
                                                                   local_bound)

    htmid_list = []
    for bound in dc2_bound:
        for ii in range(bound[0], bound[1]+1, 1):
            htmid_list.append(ii)

    print('N htmid: %d ' % len(htmid_list))
    mgr = multiprocessing.Manager()
    my_lock = mgr.Lock()
    out_dict = mgr.dict()
    for hh in htmid_list:
        assert hh not in out_dict
        out_dict[hh] = mgr.list()

    out_dict['obsHistID'] = mgr.list()
    out_dict['mjd'] = mgr.list()
    out_dict['ra'] = mgr.list()
    out_dict['dec'] = mgr.list()
    out_dict['rotTelPos'] = mgr.list()

    n_threads = 40
    chunk_size = 100000
    query = "SELECT obsHistID, descDitheredRA, descDitheredDec, "
    query += "descDitheredRotTelPos, expMJD "
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
                      args=(chunk, hs_list, my_lock, out_dict))
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

    meta_key_list = ['obsHistID', 'ra', 'dec', 'mjd', 'rotTelPos']

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

        for meta_key in ['mjd', 'ra', 'dec', 'rotTelPos']:
            out_file.create_dataset(meta_key,
                       data=np.concatenate(out_dict[meta_key])[sorted_dex])

    #with open('check_trixels.txt', 'w') as out_file:
    #    for hh in htmid_list:
    #        t = htmModule.trixelFromHtmid(hh)
    #        xyz = np.array(t.corners).transpose()
    #        ra, dec = ra_dec_from_xyz(xyz[0], xyz[1], xyz[2])
    #        for rr, dd in zip(ra, dec):
    #            out_file.write('%e %e\n' % (rr,dd))
