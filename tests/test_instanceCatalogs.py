"""
Example unit tests for sims_GCRCatSimInterface package
"""
import os
import unittest
from lsst.sims.catUtils.baseCatalogModels import StarObj
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from desc.sims.GCRCatSimInterface import make_instcat_header

class InstanceCatalogTestCase(unittest.TestCase):
    def setUp(self):
        self.obsHistID = 1418971
        obs_gen = ObservationMetaDataGenerator(database=os.environ['OPSIMDB'],
                                               driver='sqlite')
        self.obs_md \
            = obs_gen.getObservationMetaData(obsHistID=self.obsHistID)[0]
        self.outfile = 'phosim_instcat_%i.txt' % self.obsHistID

    def tearDown(self):
        os.remove(self.outfile)

    def test_make_instcat_header(self):
        star_db = StarObj(database='LSSTCATSIM',
                          host='fatboy.phys.washington.edu',
                          port=1433, driver='mssql+pymssql')
        cat = make_instcat_header(star_db, self.obs_md, self.outfile)

if __name__ == '__main__':
    unittest.main()
