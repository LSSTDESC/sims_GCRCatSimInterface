import numpy as np
import os
from lsst.sims.photUtils import cache_LSST_seds
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import AstrometryStars, PhotometryStars
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import StarObj
from lsst.sims.utils import arcsecFromRadians

import copy
import argparse


class Dc2RefCatMixin(InstanceCatalog):
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'lsst_u', 'sigma_lsst_u',
                      'lsst_g', 'sigma_lsst_g',
                      'lsst_r', 'sigma_lsst_r',
                      'lsst_i', 'sigma_lsst_i',
                      'lsst_z', 'sigma_lsst_z',
                      'lsst_y', 'sigma_lsst_y',
                      'isresolved', 'isvariable',
                      'properMotionRa', 'properMotionDec', 'parallax',
                      'radialVelocity']

    default_columns = [('isresolved', 0, int), ('isvariable', 0, int)]

    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees,
                       'properMotionDec': arcsecFromRadians,
                       'properMotionRa': arcsecFromRadians,
                       'parallax': arcsecFromRadians}

    default_formats = {'S': '%s', 'f': '%.8f', 'i': '%i'}
    override_formats = {'properMotionRa': '%.8g',
                        'properMotionDec': '%.8g',
                        'parallax': '%.8g'}



class Dc2RefCatStars(AstrometryStars, PhotometryStars, Dc2RefCatMixin):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ra', type=float, default=55.064,
                        help="Center RA of area in degrees "
                        "(default = 55.064)")
    parser.add_argument('--dec', type=float, default=-29.783,
                        help="Center Dec of area in degrees "
                        "(default = -29.783)")

    parser.add_argument('--fov', type=float, default=2.5,
                        help="Field of view radius in degrees "
                        "(default = 2.5)")
    parser.add_argument('--out_dir', type=str, default='.',
                        help="Directory where file will be made "
                        "(default = '.')")

    args = parser.parse_args()

    cache_LSST_seds(wavelen_min=200.0,
                    wavelen_max=1600.0)


    obs = ObservationMetaData(pointingRA=args.ra,
                              pointingDec=args.dec,
                              boundType='circle',
                              boundLength=args.fov,
                              bandpassName = ['u','g','r','i','z','y'])

    obs.m5 = [23.68, 24.89, 24.43, 24.0, 24.45, 22.60]

    star_db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                      port=1433, driver='mssql+pymssql')

    cat = Dc2RefCatStars(star_db, obs_metadata=obs)
    file_name = os.path.join(args.out_dir, 'dc2_reference_catalog.txt')
    cat.write_catalog(file_name, chunk_size=10000)
