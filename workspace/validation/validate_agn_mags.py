import os
import sqlite3
import numpy as np
import GCRCatalogs
from GCR import GCRQuery

from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels

import argparse

if __name__ == "__main__":

    project_dir = os.path.join('/global/projecta/projectdirs',
                               'lsst/groups/SSim/DC2/cosmoDC2_v1.1.4')

    default_agn_db = os.path.join(project_dir,
                                  'agn_db_mbh7_mi30_sf4.db')

    parser = argparse.ArgumentParser()
    parser.add_argument('--cat_dir', type=str, default=None,
                        help='parent directory of $obsHistID/')
    parser.add_argument('--obs', type=int, default=None,
                        help='obsHistID of pointing')
    parser.add_argument('--agn_db', type=str,
                        default=default_agn_db,
                        help='Name of agn parameters db\n'
                        '(default %s)' % default_agn_db)

    args = parser.parse_args()

    if not os.path.isfile(args.agn_db):
        raise RuntimeError('\n%s\nis not a file\n' % args.agn_db)
