import unittest
import os
import tempfile
import shutil

from lsst.utils import getPackageDir
from desc.sims.GCRCatSimInterface import diskDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import PhoSimDESCQA
from lsst.sims.utils import ObservationMetaData

mag_grid = os.path.join(getPackageDir('sims_GCRCatSimInterface'), 'data',
                        'CatSimMagGrid.txt')


@unittest.skipIf(not os.path.exists(mag_grid),
                 'Have not created SED magnitude grid, yet')
class BulgePhoSimCatalogTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.out_dir = tempfile.mkdtemp(prefix='full_instanceCatalog')

    @classmethod
    def tearDown(cls):
        list_of_files = os.listdir(cls.out_dir)
        for file_name in list_of_files:
            os.unlink(os.path.join(cls.out_dir, file_name))
        shutil.rmtree(cls.out_dir)

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
        cat_name = os.path.join(self.out_dir, 'disk_phosim_cat.txt')
        cat.write_catalog(cat_name)
        with open(cat_name, 'r') as in_file:
            cat_lines = in_file.readlines()

        self.assertGreater(len(cat_lines), 100)

        if os.path.exists(cat_name):
            os.unlink(cat_name)


if __name__ == "__main__":
    unittest.main()
