from lsst.sims.catalogs.definitions import CompoundInstanceCatalog
from desc.sims.GCRCatSimInterface import bulgeDESCQAObject
from desc.sims.GCRCatSimInterface import diskDESCQAObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.utils import ObservationMetaData

class CatForBulge(InstanceCatalog):
    cannot_be_null = ['sedFilename']
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'sedFilename', 'majorAxis']

class CatForDisk(InstanceCatalog):
    cannot_be_null = ['sedFilename']
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'sedFilename']

if __name__ == "__main__":

    obs = ObservationMetaData(pointingRA=0.0, pointingDec=0.0,
                              boundType='circle', boundLength=0.01)

    cat = CompoundInstanceCatalog([CatForBulge, CatForDisk],
                                  [bulgeDESCQAObject,
                                   diskDESCQAObject],
                                  obs_metadata=obs)

    cat.write_catalog('descqa_compound_cat.txt', chunk_size=10000)
