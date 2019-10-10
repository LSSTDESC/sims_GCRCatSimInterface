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

def make_lookup(chunk, my_lock, output_dict):
    """
    chunk will be the result of sqlite3.fetchmany() on the OpSim database;

    each row will be (obsHistID, descDitheredRA, descDitheredDec)

    output_dict will need to be initialized with the valid values of htmid
    """
    radius = 2.1  # in degrees
    ne_corner = (71.46, -27.25)
    nw_corner = (52.25, -27.25)
    se_corner = (73.79, -44.33)
    sw_corner = (49.92, -44.33)
    pt_inside = (60.0, -30.0)

    hs_list = []
    hs_list.append(htmModule.halfSpaceFromPts(ne_corner, nw_corner, pt_inside))
    hs_list.append(htmModule.halfSpaceFromPts(ne_corner, se_corner, pt_inside))
    hs_list.append(htmModule.halfSpaceFromPts(se_corner, sw_corner, pt_inside))
    hs_list.append(htmModule.halfSpaceFromPts(sw_corner, nw_corner, pt_inside))

    local_lookup = {}
    local_dex = {}
    for pointing in chunk:
        obs_id = int(pointing[0])
        obs_hs = htmModule.halfSpaceFromRaDec(np.degrees(float(pointing[1])),
                                              np.degrees(float(pointing[2])),
                                              radius)

        is_valid_pointing = True
        for hs in hs_list:
            if not hs.intersects_circle(obs_hs.vector, obs_hs.phi):
                is_valid_pointing = False
                break

        if not is_valid_pointing:
            continue

        obs_bounds = obs_hs.findAllTrixels(6)
        for bound in obs_bounds:
            for htmid in range(bound[0], bound[1]+1, 1):
                if htmid not in local_lookup:
                    local_lookup[htmid] = -1*np.ones(len(chunk), dtype=int)
                    local_dex[htmid] = 0
                local_lookup[htmid][local_dex[htmid]] = obs_id
                local_dex[htmid] += 1

    with my_lock:
        for htmid in local_lookup:
            if htmid in output_dict:
                valid = np.where(local_lookup[htmid]>=0)
                output_dict[htmid].append(local_lookup[htmid][valid])

if __name__ == "__main__":
    pass
