import os
import sqlite3
import numpy as np
import GCRCatalogs

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--cosmoDC2', type=str,
                        default='cosmoDC2_v1.1.4_image',
                        help='extragalactic catalog to load '
                        '(default: cosmoDC2_v1.1.4_image)')
    parser.add_argument('--sn_db', type=str,
                        default=None,
                        help='SNe parameter database to check '
                        'against')

    args = parser.parse_args()

    if not os.path.isfile(args.sn_db):
        raise RuntimeError('\n%s\nis not a file\n')
