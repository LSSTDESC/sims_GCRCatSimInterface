import os
import numpy as np
import sqlite3
import json
import copy

from lsst.sims.utils import findHtmid
from lsst.utils import getPackageDir
from lsst.sims.catUtils.supernovae import SNObject
from desc.sims.GCRCatSimInterface.TruthCatalogUtils import _truth_trixel_level


__all__ = ["add_unsprinkled_sne_params"]

_simsGCR_data_dir = os.path.join(getPackageDir('sims_GCRCatSimInterface'), 'data')

def add_unsprinkled_sne_params(sqlite_file,
                               list_of_sne_files=[os.path.join(_simsGCR_data_dir, 'uDDF_hostlessSN_trimmed.csv'),
                                                  os.path.join(_simsGCR_data_dir, 'MainSurvey_hostlessSN_trimmed.csv'),
                                                  os.path.join(_simsGCR_data_dir, 'MainSurvey_hostlessSN_highz_trimmed.csv'),
                                                  os.path.join(_simsGCR_data_dir, 'uDDFHostedSNPositions_trimmed.csv'),
                                                  os.path.join(_simsGCR_data_dir, 'MainSurveyHostedSNPositions_trimmed.csv')]):
    """
    Parameters
    ----------
    sqlite_file -- path to sqlite database that will contain the photometry
    parameters for the astrophysical objects being simulated (must already
    exist)

    list_of_sne_files -- a list of paths to csv files containing the
    supernova parameters.  This is the way supernova parameters were stored
    in protoDC2.  We will have to do something different for cosmoDC2.
    Note: these paths must be listed in the same order they are listed in
    InstanceCatalogWriter.py in order for the objects uniqueIDs to get
    set the same way they are set in the InstanceCatalogs (though this
    may not actually matter)

    Returns
    -------
    Nothing.  This method will just add rows to sqlite_file corresponding to
    the supernovae whose parameters are in list_of_sne_files
    """

    if not os.path.isfile(sqlite_file):
        raise RuntimeError("\n%s\nis not a file\n" % sqlite_file)

    sne_obj_id = 42  # offset for uniqueID

    dtype = np.dtype([('snid', int), ('x0', float), ('t0', float),
                      ('x1', float), ('c', float), ('z', float),
                      ('snra', float), ('sndec', float)])

    base_sn_obj = SNObject(0., 0.)

    with sqlite3.connect(sqlite_file) as conn:
        cursor = conn.cursor()
        for i_cat, csv_file_name in enumerate(list_of_sne_files):
            print('processing %s' % csv_file_name)
            sne_data = np.genfromtxt(csv_file_name, skip_header=1,
                                     delimiter=',', dtype=dtype)

            htmid_arr = findHtmid(sne_data['snra'], sne_data['sndec'],
                                  _truth_trixel_level)

            sn_param_arr = []
            for i_obj in range(len(sne_data)):
                sn_dict = copy.deepcopy(base_sn_obj.SNstate)
                sn_dict['_ra'] = np.radians(sne_data['snra'][i_obj])
                sn_dict['_dec'] = np.radians(sne_data['sndec'][i_obj])
                sn_dict['z'] = sne_data['z'][i_obj]
                sn_dict['c'] = sne_data['c'][i_obj]
                sn_dict['x0'] = sne_data['x0'][i_obj]
                sn_dict['x1'] = sne_data['x1'][i_obj]
                sn_dict['t0'] = sne_data['t0'][i_obj]
                sn_obj = SNObject.fromSNState(sn_dict)
                sn_obj.mwEBVfromMaps()
                sn_param_arr.append(json.dumps(sn_obj.SNstate))

            values = ((int(sne_data['snid'][i_obj]*1024+i_cat+sne_obj_id),
                       -1, int(htmid_arr[i_obj]),
                       sne_data['snra'][i_obj], sne_data['sndec'][i_obj],
                       sne_data['z'][i_obj], 'None', -1.0, 'None',
                       sn_param_arr[i_obj],
                       sne_data['t0'][i_obj],
                       0,1,0,0.0,0.0,0.0)
                      for i_obj in range(len(sne_data)))

            cursor.executemany('''INSERT INTO zpoint
                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',values)
            conn.commit()
