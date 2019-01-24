import os
import sqlite3
import gzip
import numpy as np

import lsst.sims.utils.htmModule as htm
from lsst.sims.utils import angularSeparation
from lsst.sims.catUtils.supernovae import SNObject

import argparse


def validate_sne(cat_dir, obsid, fov_deg=2.1):
    """
    Parameters
    ----------
    cat_dir is the parent dir of $obsid

    obsid is the obsHistID of the pointing (an int)

    fov_deg is the radius of the field of view in degrees
    """

    project_dir = os.path.join('/global/projecta/projectdirs',
                               'lsst/groups/SSim/DC2/cosmoDC2_v1.1.4')

    sne_db_name = os.path.join(project_dir, 'sne_cosmoDC2_v1.1.4_MS_DDF.db')
    if not os.path.isfile(sne_db_name):
        raise RuntimeError("\n%s\nis not a file" % sne_db_name)

    instcat_dir = os.path.join(cat_dir, '%.8d' % obsid)
    if not os.path.isdir(instcat_dir):
        raise RuntimeError("\n%s\nis not a dir" % instcat_dir)

    sne_name = os.path.join(instcat_dir, 'sne_cat_%d.txt.gz' % obsid)
    if not os.path.isfile(sne_name):
        raise RuntimeError("\n%s\nis not a file" % sne_name)

    phosim_name = os.path.join(instcat_dir, 'phosim_cat_%d.txt' % obsid)
    if not os.path.isfile(phosim_name):
        raise RuntimeError("\n%s\nis not a file" % phosim_name)

    sed_dir = os.path.join(instcat_dir, "Dynamic")
    if not os.path.isdir(sed_dir):
        raise RuntimeError("\n%s\nis not a dir" % sed_dir)

    opsim_db = os.path.join("/global/projecta/projectdirs",
                            "lsst/groups/SSim/DC2",
                            "minion_1016_desc_dithered_v4_sfd.db")
    if not os.path.isfile(opsim_db):
        raise RuntimeError("\n%s\nis not a file" % opsim_db)

    with sqlite3.connect(opsim_db) as conn:
        c = conn.cursor()
        r = c.execute("SELECT descDitheredRA, descDitheredDec FROM Summary "
                      "WHERE obsHistID==%d" % obsid).fetchall()
        pointing_ra = np.degrees(r[0][0])
        pointing_dec = np.degrees(r[0][1])

    with sqlite3.connect(sne_db_name) as conn:
        c = conn.cursor()
        query = "SELECT snid_in, snra_in, sndec_in, "
        query += "c_in, mB, t0_in, x0_in, x1_in, z_in "
        query += "FROM sne_params WHERE ("

        where_clause = None
        half_space = htm.halfSpaceFromRaDec(pointing_ra, pointing_dec, fov_deg)
        bounds = half_space.findAllTrixels(6)

        for bound in bounds:
            if where_clause is not None:
                where_clause += " OR ("
            else:
                where_clause = "("

            if bound[0] == bound[1]:
                where_clause += "htmid_level_6==%d" % bound[0]
            else:
                where_clause += "htmid_level_6>=%d AND htmid_level_6<=%d" % (bound[0],bound[1])

            where_clause += ")"

        query += where_clause+")"

        sn_params = c.execute(query).fetchall()
        print(len(sn_params))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--cat_dir', type=str, default=None,
                        help='parent dir of $obs')
    parser.add_argument('--obs', type=int, default=None,
                        help='obsHistID of the pointing')
    parser.add_argument('--fov_deg', type=float, default=2.1,
                        help='radius of field of view in degrees')

    args = parser.parse_args()

    validate_sne(args.cat_dir, args.obs, fov_deg=args.fov_deg)
