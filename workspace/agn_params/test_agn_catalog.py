import os
import numpy as np

from lsst.sims.catalogs.definitions import InstanceCatalog
from desc.sims.GCRCatSimInterface import agnDESCQAObject_protoDC2
from lsst.sims.utils import ObservationMetaData

class VarParCat(InstanceCatalog):
    column_outputs = ['galaxy_id', 'raJ2000', 'decJ2000',
                      'magNorm', 'varParamStr']


    transformations = {'raJ2000': np.degrees,
                       'decJ2000': np.degrees}

if __name__ == "__main__":

    ra = 112.0
    dec = -75.0

    db = agnDESCQAObject_protoDC2(yaml_file_name='proto-dc2_v2.1.2')

    agn_db_name = os.path.join(os.environ['SCRATCH'], 'proto_dc2_agn',
                               'test_agn.db')

    assert os.path.exists(agn_db_name)

    db.agn_params_db = agn_db_name
    db.field_ra = ra
    db.field_dec = dec

    obs = ObservationMetaData(pointingRA=ra, pointingDec=dec,
                              boundType='circle',
                              boundLength=1.0)

    cat = VarParCat(db, obs_metadata=obs)
    cat.write_catalog('agn_junk.txt', chunk_size=1000000)
