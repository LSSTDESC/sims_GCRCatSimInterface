import unittest
import os

from desc.sims.GCRCatSimInterface import diskDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import PhoSimDESCQA
from lsst.sims.utils import ObservationMetaData


class BulgePhoSimCatalogTestCase(unittest.TestCase):

    def setUp(self):
        self.cat_name = 'dummy_disk_catalog.txt'

    def tearDown(self):
        if os.path.exists(self.cat_name):
            os.unlink(self.cat_name)

    def test_dummy_disk_catalog(self):
        """
        Just try producing a PhoSim InstanceCatalog from a fake
        ObservationMetaData, using protoDC2 (to make sure we don't
        break the whole interface)
        """
        db = diskDESCQAObject_protoDC2(yaml_file_name='proto-dc2_v2.1.2')
        db.field_ra = 13.8
        db.field_dec = -45.6
        obs = ObservationMetaData(pointingRA=14.0, pointingDec=-45.0,
                                  mjd=59580.0, rotSkyPos=112.0,
                                  bandpassName='z',
                                  boundType='circle', boundLength=0.01)
        cat = PhoSimDESCQA(db, obs_metadata=obs, cannot_be_null=['hasDisk'])
        cat.phoSimHeaderMap = {}
        cat.write_catalog(self.cat_name)
        with open(self.cat_name, 'r') as in_file:
            cat_lines = in_file.readlines()

        self.assertGreater(len(cat_lines), 100)

if __name__ == "__main__":
    unittest.main()
