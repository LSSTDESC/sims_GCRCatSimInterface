import os
import sqlite3
import gzip
import numpy as np

import lsst.sims.utils.htmModule as htm
from lsst.sims.utils import angularSeparation
from lsst.sims.catUtils.supernovae import SNObject

import argparse


def validate_sne(cat_dir, obsid):
    """
    Parameters
    ----------
    cat_dir is the parent dir of $obsid

    obsid is the obsHistID of the pointing (an int)
    """

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

    print(pointing_ra, pointing_dec)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--cat_dir', type=str, default=None,
                        help='parent dir of $obs')
    parser.add_argument('--obs', type=int, default=None,
                        help='obsHistID of the pointing')

    args = parser.parse_args()

    validate_sne(args.cat_dir, args.obs)
