import os
import sqlite3
import numpy as np
from lsst.sims.utils import htmModule as htm

fname = 'minion_1016_desc_dithered_v4_sfd.db'
assert os.path.isfile(fname)

out_name = 'list_of_dc2_ddf_obsHistId.txt'

ddf_corners = {}
ddf_corners['ne'] = (53.764, -27.533)
ddf_corners['nw'] = (52.486, -27.533)
ddf_corners['se'] = (53.771, -28.667)
ddf_corners['sw'] = (52.479, -28.667)

ddf_center = (53.125, -28.100)
ra_rad = np.radians(ddf_center[0])
dec_rad = np.radians(ddf_center[1])
ddf_center_cart = np.array([np.cos(dec_rad)*np.cos(ra_rad),
                            np.cos(dec_rad)*np.sin(ra_rad),
                            np.sin(dec_rad)])

ddf_corners_cart = {}
for name in ddf_corners:
    ra_rad = np.radians(ddf_corners[name][0])
    dec_rad = np.radians(ddf_corners[name][1])
    vv = np.array([np.cos(dec_rad)*np.cos(ra_rad),
                   np.cos(dec_rad)*np.sin(ra_rad),
                   np.sin(dec_rad)])
    ddf_corners_cart[name] = vv

ddf_edges = []
ddf_edges.append((ddf_corners_cart['nw'], ddf_corners_cart['ne']))
ddf_edges.append((ddf_corners_cart['nw'], ddf_corners_cart['sw']))
ddf_edges.append((ddf_corners_cart['ne'], ddf_corners_cart['se']))
ddf_edges.append((ddf_corners_cart['se'], ddf_corners_cart['sw']))

ddf_halfspaces = []
ddf_halfspaces.append(htm.halfSpaceFromRaDec(0.0, -90.0, 62.467))
ddf_halfspaces.append(htm.halfSpaceFromRaDec(0.0, 90.0, 118.667))
ddf_halfspaces.append(htm.halfSpaceFromPoints(ddf_corners['ne'],
                                             ddf_corners['se'],
                                             ddf_center))
ddf_halfspaces.append(htm.halfSpaceFromPoints(ddf_corners['nw'],
                                              ddf_corners['sw'],
                                              ddf_center))

with open(out_name, 'w') as out_file:
    with sqlite3.connect(fname) as connection:
        c = connection.cursor()
        query = "SELECT obsHistID, descDitheredRA, descDitheredDec "
        query += "FROM Summary GROUP BY obsHistID ORDER BY obsHistID"
        results = c.execute(query).fetchall()

    ct_in = 0
    for i_ptng, ptng in enumerate(results):
        if i_ptng>0 and i_ptng%10000==0:
            print('done %d; inside %d' % (i_ptng, ct_in))
        obs_id = int(ptng[0])
        ra = float(ptng[1])
        dec = float(ptng[2])
        if np.abs(ra)>10.0 or np.abs(dec)>10.0:
            raise RuntimeError("I think these are in degrees")

        vv = np.array([np.cos(dec)*np.cos(ra),
                       np.cos(dec)*np.sin(ra),
                       np.sin(dec)])
        ra = np.degrees(ra)
        dec = np.degrees(dec)

        # check if center of pointing is in the DDF
        overlaps = True
        for hs in ddf_halfspaces:
            if not hs.contains_pt(vv):
                overlaps = False
                break

        if not overlaps:
            ptng_hs = htm.halfSpaceFromRaDec(ra, dec, 2.1)

            # see if center of ddf is in the pointing
            if ptng_hs.contains_pt(ddf_center_cart):
                overlaps = True

            # see if pointing HalfSpace contains one of the DDF corners
            if not overlaps:
                for name in ddf_corners_cart:
                    if ptng_hs.contains_pt(ddf_corners_cart[name]):
                        overlaps = True
                        break

            # see if pointing HalfSpace intersects an edge
            # of the ddf
            if not overlaps:
                for edge in ddf_edges:
                    if ptng_hs.intersects_edge(edge[0], edge[1]):
                        overlaps = True
                        break

        if overlaps:
            out_file.write('%d\n' % obs_id)
            ct_in += 1
