import numpy as np
import os
from lsst.sims.catalogs.decorators import compound
from lsst.sims.photUtils import cache_LSST_seds
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import AstrometryStars, PhotometryStars
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import StarObj
from lsst.sims.utils import arcsecFromRadians

import copy
import argparse


class Dc2RefCatMixin(object):
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'lsst_u', 'sigma_lsst_u',
                      'lsst_g', 'sigma_lsst_g',
                      'lsst_r', 'sigma_lsst_r',
                      'lsst_i', 'sigma_lsst_i',
                      'lsst_z', 'sigma_lsst_z',
                      'lsst_y', 'sigma_lsst_y',
                      'lsst_u_smeared',
                      'lsst_g_smeared',
                      'lsst_r_smeared',
                      'lsst_i_smeared',
                      'lsst_z_smeared',
                      'lsst_y_smeared',
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

    @compound('sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r',
              'sigma_lsst_i', 'sigma_lsst_z', 'sigma_lsst_y')
    def get_lsst_photometric_uncertainties(self):
        n_obj = len(self.column_by_name('uniqueId'))
        fiducial_val = 0.01
        arr = fiducial_val*np.ones(n_obj, dtype=float)
        return np.array([arr]*6)

    def _smear_photometry(self, mag_name, sigma_name):
        """
        Parameters
        ----------
        mag_name: str
            name of the column containing the true magnitude

        sigma_name: str
            name of the column containing the sigma in the magnitude

        Returns
        -------
        smeared_mag: nd.array
            array of magnitudes smeared by sigma
        """
        if not hasattr(self, '_phot_rng'):
            self._phot_rng = np.random.RandomState(self.db_obj.objectTypeId)
        mag = self.column_by_name(mag_name)
        sig = self.column_by_name(sigma_name)
        offset = self._phot_rng.normal(loc=0.0,
                                       scale=1.0,
                                       size=len(sig))

        return mag + offset*sig

    @compound('lsst_u_smeared', 'lsst_g_smeared', 'lsst_r_smeared',
              'lsst_i_smeared', 'lsst_z_smeared', 'lsst_y_smeared')
    def get_smeared_photometry(self):
        return np.array([self._smear_photometry('lsst_u', 'sigma_lsst_u'),
                         self._smear_photometry('lsst_g', 'sigma_lsst_g'),
                         self._smear_photometry('lsst_r', 'sigma_lsst_r'),
                         self._smear_photometry('lsst_i', 'sigma_lsst_i'),
                         self._smear_photometry('lsst_z', 'sigma_lsst_z'),
                         self._smear_photometry('lsst_y', 'sigma_lsst_y')])


class Dc2RefCatStars(Dc2RefCatMixin, AstrometryStars, PhotometryStars,
                     InstanceCatalog):
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

    cache_LSST_seds(wavelen_min=0.0,
                    wavelen_max=1600.0)


    obs = ObservationMetaData(pointingRA=args.ra,
                              pointingDec=args.dec,
                              boundType='circle',
                              boundLength=args.fov)

    star_db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                      port=1433, driver='mssql+pymssql')

    cat = Dc2RefCatStars(star_db, obs_metadata=obs)
    file_name = os.path.join(args.out_dir, 'dc2_reference_catalog.txt')
    cat.write_catalog(file_name, chunk_size=10000)
