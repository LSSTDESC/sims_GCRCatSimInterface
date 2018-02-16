import numpy as np
import os
from lsst.sims.photUtils import cache_LSST_seds
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import AstrometryStars, PhotometryStars
import copy

class Dc1RefCat(InstanceCatalog, AstrometryStars, PhotometryStars):
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z',
                      'lsst_y', 'isresolved', 'isvariable']
    default_columns = [('isresolved', 0, int), ('isvariable', 0, int)]
    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees}
    default_formats = {'S': '%s', 'f': '%.8f', 'i': '%i'}


output_dir = os.path.join("/astro", "users", "danielsf", "dc1_catalogs")

cache_LSST_seds()

from lsst.sims.utils import ObservationMetaData

obs = ObservationMetaData(pointingRA=np.degrees(1.641324),
                          pointingDec=np.degrees(-0.496321),
                          boundType='circle',
                          boundLength=8.0)

from lsst.sims.catUtils.baseCatalogModels import StarObj

star_db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

cat = Dc1RefCat(star_db, obs_metadata=obs)
file_name = os.path.join(output_dir, 'dc1_reference_catalog_8deg_radius.txt')
cat.write_catalog(file_name, chunk_size=10000)
