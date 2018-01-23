from desc.sims.GCRCatSimInterface import CompoundDESCQAInstanceCatalog
from desc.sims.GCRCatSimInterface import bulgeDESCQAObject
from desc.sims.GCRCatSimInterface import diskDESCQAObject
from desc.sims.GCRCatSimInterface import CompoundDESCQACatalogDBObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.utils import ObservationMetaData

class _testDESCQAObj(object):
    yaml_file_name = 'proto-dc2_v2.1.2'

class bulgeDESCQAObject_test(_testDESCQAObj, bulgeDESCQAObject):
    objid = 'bulge_descqa'

class diskDESCQAObject_test(_testDESCQAObj, diskDESCQAObject):
    objid = 'disk_descqa'

class CatForBulge(InstanceCatalog):
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'minorAxis', 'majorAxis']
    transformations = {'minorAxis': arcsecFromRadians,
                       'majorAxis': arcsecFromRadians}

class CatForDisk(InstanceCatalog):
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'minorAxis']
    transformations = {'minorAxis': arcsecFromRadians}

if __name__ == "__main__":

    obs = ObservationMetaData(pointingRA=0.0, pointingDec=0.0,
                              boundType='circle', boundLength=0.01)

    cat = CompoundDESCQAInstanceCatalog([CatForBulge, CatForDisk],
                                        [bulgeDESCQAObject_test,
                                         diskDESCQAObject_test],
                                        obs_metadata=obs,
                                        compoundDBclass=CompoundDESCQACatalogDBObject)

    cat.write_catalog('descqa_compound_cat.txt', chunk_size=10000)
