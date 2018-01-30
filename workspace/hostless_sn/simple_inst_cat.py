import numpy as np
import os

from lsst.sims.catalogs.db import fileDBObject
from lsst.sims.utils import ObservationMetaData
from lsst.utils import getPackageDir
from lsst.sims.catalogs.definitions import InstanceCatalog

class SimpleCatalogClass(InstanceCatalog):
    column_outputs = ['raJ2000', 'decJ2000', 'z']
    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees}


class SNFileDBObject(fileDBObject):
    columns = [('raJ2000', 'snra*PI()/180.'),
               ('decJ2000', 'sndec*PI()/180.'),
               ('Tt0', 't0'),
               ('Tx0', 'x0'),
               ('Tx1', 'x1'),
               ('Tc', 'c'),
               ('Tsnid', 'snid'),
               ('redshift', 'z')
              ]

    dbDefaultValues = {'varsimobjid':-1,
                       'runid':-1,
                       'ismultiple':-1,
                       'run':-1,
                       'runobjid':-1}

if __name__ == "__main__":

    database_file_name = os.path.join(getPackageDir('sims_GCRCatSimInterface'),
                                      'data', 'uDDF_hostlessSN.csv')

    assert os.path.exists(database_file_name)

    db = SNFileDBObject(database_file_name, idColKey='snid')
    db.raColName = 'snra'
    db.decColName = 'sndec'
    db.objectTypeId = 612

    obs = ObservationMetaData(pointingRA=55.0, pointingDec=-28.0,
                              boundType='circle', boundLength=0.5)

    cat = SimpleCatalogClass(db, obs_metadata=obs)
    cat.write_catalog('simple_catalog.txt', chunk_size=10000)
