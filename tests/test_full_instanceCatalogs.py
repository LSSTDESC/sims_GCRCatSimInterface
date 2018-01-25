import unittest
import os
import tempfile
import shutil

from lsst.utils import getPackageDir
from desc.sims.GCRCatSimInterface import diskDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import PhoSimDESCQA
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.definitions import InstanceCatalog

mag_grid = os.path.join(getPackageDir('sims_GCRCatSimInterface'), 'data',
                        'CatSimMagGrid.txt')


class BulgePhoSimCatalogTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.out_dir = tempfile.mkdtemp(prefix='full_instanceCatalog')
        cls.field_ra = 82.3
        cls.field_dec = -81.4

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.out_dir):
            list_of_files = os.listdir(cls.out_dir)
            for file_name in list_of_files:
                os.unlink(os.path.join(cls.out_dir, file_name))
            shutil.rmtree(cls.out_dir)

    @unittest.skipIf(not os.path.exists(mag_grid),
                     'Have not created SED magnitude grid, yet')
    def test_disk_phosim_catalog(self):
        """
        Just try producing a PhoSim InstanceCatalog from a fake
        ObservationMetaData, using protoDC2 (to make sure we don't
        break the whole interface)
        """
        db = diskDESCQAObject_protoDC2(yaml_file_name='proto-dc2_v2.1.2')
        db.field_ra = self.field_ra
        db.field_dec = self.field_dec
        obs = ObservationMetaData(pointingRA=self.field_ra+0.2,
                                  pointingDec=self.field_dec-0.2,
                                  mjd=59580.0, rotSkyPos=112.0,
                                  bandpassName='z',
                                  boundType='circle', boundLength=0.01)

        cat = PhoSimDESCQA(db, obs_metadata=obs, cannot_be_null=['hasDisk'])
        cat.phoSimHeaderMap = {}
        cat_name = os.path.join(self.out_dir, 'disk_phosim_cat.txt')
        self.assertTrue(os.path.exists(self.out_dir))
        cat.write_catalog(cat_name)
        with open(cat_name, 'r') as in_file:
            cat_lines = in_file.readlines()

        self.assertGreater(len(cat_lines), 50)

        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def test_default_varParamStr(self):
        """
        Test that DESCQAObjects now return varParamStr='None' by default
        """
        db = diskDESCQAObject_protoDC2(yaml_file_name='proto-dc2_v2.1.2')
        db.field_ra = self.field_ra
        db.field_dec = self.field_dec

        obs = ObservationMetaData(pointingRA=self.field_ra-0.7,
                                  pointingDec=self.field_dec+1.0)

        class VarParamStrTestClass(InstanceCatalog):
            column_outputs = ['raJ2000', 'decJ2000', 'varParamStr']

        cat = VarParamStrTestClass(db, obs_metadata=obs)
        cat_name = os.path.join(self.out_dir, 'varParamStr_cat.txt')
        cat.write_catalog(cat_name)
        line_ct = 0
        with open(cat_name, 'r') as in_file:
            for line in in_file:
                if line[0] == '#':
                    continue
                cols = line.strip().split()
                self.assertEqual(cols[2],'None')
                line_ct += 1
        self.assertGreater(line_ct, 100)

        if os.path.exists(cat_name):
            os.unlink(cat_name)

if __name__ == "__main__":
    unittest.main()
