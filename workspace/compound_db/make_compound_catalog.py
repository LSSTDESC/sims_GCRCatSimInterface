from lsst.sims.catUtils.baseCatalogModels import GalaxyBulgeObj
from lsst.sims.catUtils.baseCatalogModels import GalaxyDiskObj
from lsst.sims.catUtils.baseCatalogModels import GalaxyTileCompoundObj
from lsst.sims.catalogs.definitions import CompoundInstanceCatalog
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.utils import ObservationMetaData

class _fatboy(object):
    database = 'LSSTCATSIM'
    host = 'fatboy.phys.washington.edu'
    port = 1433
    driver = 'mssql+pymssql'

class GalaxyBulgeObj_fatboy(_fatboy, GalaxyBulgeObj):
    pass

class GalaxyDiskObj_fatboy(_fatboy, GalaxyDiskObj):
    pass

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
                                  [GalaxyBulgeObj_fatboy,
                                   GalaxyDiskObj_fatboy],
                                  obs_metadata=obs,
                                  compoundDBclass=GalaxyTileCompoundObj)

    cat.write_catalog('compound_cat.txt', chunk_size=10000)
