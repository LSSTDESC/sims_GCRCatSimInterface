import os
import numpy as np

from desc.sims.GCRCatSimInterface import agnDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import PhoSimDESCQA_AGN
from lsst.sims.utils import ObservationMetaData

if __name__ == "__main__":

    ra = 112.0
    dec = -75.0

    db = agnDESCQAObject_protoDC2(yaml_file_name='protoDC2')

    agn_db_name = os.path.join(os.environ['SCRATCH'], 'proto_dc2_agn',
                               'test_agn.db')

    assert os.path.exists(agn_db_name)

    db.agn_params_db = agn_db_name
    db.field_ra = ra
    db.field_dec = dec

    obs1 = ObservationMetaData(pointingRA=ra, pointingDec=dec,
                               rotSkyPos=11.0,
                               mjd=59585.0,
                               bandpassName='g',
                               boundType='circle',
                               boundLength=0.1)

    obs2 =  ObservationMetaData(pointingRA=ra, pointingDec=dec,
                               rotSkyPos=11.0,
                               mjd=59785.0,
                               bandpassName='g',
                               boundType='circle',
                               boundLength=0.1)

    cat1 = PhoSimDESCQA_AGN(db, obs_metadata=obs1)
    cat1.phoSimHeaderMap = {}
    cat1.write_catalog('phosim_agn_cat_1.txt', chunk_size=1000000)

    cat2 = PhoSimDESCQA_AGN(db, obs_metadata=obs2)
    cat2.phoSimHeaderMap = {}
    cat2.write_catalog('phosim_agn_cat_2.txt', chunk_size=1000000)

