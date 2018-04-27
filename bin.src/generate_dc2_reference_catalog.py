import sqlite3
import numpy as np
import os
from lsst.utils import getPackageDir
from lsst.sims.catalogs.decorators import compound, cached
from lsst.sims.photUtils import cache_LSST_seds
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import AstrometryStars, PhotometryStars
from lsst.sims.catUtils.mixins import AstrometryGalaxies, EBVmixin
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import StarObj
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.utils import cartesianFromSpherical
from lsst.sims.utils import sphericalFromCartesian

from lsst.sims.photUtils import Sed, BandpassDict
from lsst.sims.photUtils import getImsimFluxNorm
from lsst.sims.utils import defaultSpecMap
from lsst.sims.utils import radiansFromArcsec

from desc.sims.GCRCatSimInterface import DESCQACatalogMixin
from desc.sims.GCRCatSimInterface import DESCQAObject_protoDC2

import copy
import argparse


class DESCQAReferenceObject(DESCQAObject_protoDC2):
    objectTypeId = 231
    objid = 'ref_descqa'
    _columns_need_postfix = False


class Dc2RefCatMixin(object):
    column_outputs = ['uniqueId',
                      'raJ2000', 'decJ2000',
                      'sigma_raJ2000', 'sigma_decJ2000',
                      'raJ2000_smeared', 'decJ2000_smeared',
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
                      'isresolved', 'isagn',
                      'properMotionRa', 'properMotionDec', 'parallax',
                      'radialVelocity']


    transformations = {'raJ2000': np.degrees,
                       'decJ2000': np.degrees,
                       'sigma_raJ2000': np.degrees,
                       'sigma_decJ2000': np.degrees,
                       'raJ2000_smeared': np.degrees,
                       'decJ2000_smeared': np.degrees,
                       'properMotionDec': arcsecFromRadians,
                       'properMotionRa': arcsecFromRadians,
                       'parallax': arcsecFromRadians}

    default_formats = {'S': '%s', 'f': '%.8f', 'i': '%i'}

    override_formats = {'properMotionRa': '%.17f',
                        'properMotionDec': '%.17f',
                        'parallax': '%.17f',
                        'raJ2000': '%.17f',
                        'decJ2000': '%.17f',
                        'sigma_raJ2000': '%.17f',
                        'sigma_decJ2000': '%.17f',
                        'raJ2000_smeared': '%.17f',
                        'decJ2000_smeared': '%.17f'}

    @compound('sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r',
              'sigma_lsst_i', 'sigma_lsst_z', 'sigma_lsst_y')
    def get_lsst_photometric_uncertainties(self):
        n_obj = len(self.column_by_name('uniqueId'))
        fiducial_val = 0.001
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

    @compound('sigma_raJ2000', 'sigma_decJ2000')
    def get_astrometric_uncertainties(self):
        fiducial_val = radiansFromArcsec(0.0001)  # in radians
        n_obj = len(self.column_by_name('uniqueId'))
        return np.array([fiducial_val*np.ones(n_obj, dtype=float),
                         fiducial_val*np.ones(n_obj, dtype=float)])

    def _smear_astrometry(self, ra, dec, err):
        """
        Parameters
        ----------
        ra: ndarray
            in radians

        dec: ndarray
            in radians

        err: ndarray
            the angle uncertainty (in radians)

        Returns
        -------
        smeared ra/dec: ndarrays
            in radians
        """
        if not hasattr(self, '_astrometry_rng'):
            self._astrometry_rng = np.random.RandomState(self.db_obj.objectTypeId+288)

        xyz = cartesianFromSpherical(ra, dec)

        delta = self._astrometry_rng.normal(0.0, 1.0, len(ra))*err

        sin_d = np.sin(delta)
        cos_d = np.cos(delta)
        random_vec = self._astrometry_rng.normal(0.0,1.0,(len(xyz), 3))

        new_xyz = np.zeros(xyz.shape, dtype=float)
        for i_obj in range(len(xyz)):
            dot_product = np.dot(xyz[i_obj], random_vec[i_obj])
            random_vec[i_obj] -= dot_product*xyz[i_obj]
            norm = np.sqrt(np.dot(random_vec[i_obj], random_vec[i_obj]))
            random_vec[i_obj]/=norm

            new_xyz[i_obj] = cos_d[i_obj]*xyz[i_obj] + sin_d[i_obj]*random_vec[i_obj]

            norm = np.sqrt(np.dot(new_xyz[i_obj],new_xyz[i_obj]))

            new_xyz[i_obj]/=norm

        return sphericalFromCartesian(new_xyz)

    @compound('raJ2000_smeared', 'decJ2000_smeared')
    def get_smeared_astrometry(self):
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')
        err = self.column_by_name('sigma_raJ2000')
        return self._smear_astrometry(ra, dec, err)


class Dc2RefCatStars(Dc2RefCatMixin, AstrometryStars, PhotometryStars,
                     InstanceCatalog):

    default_columns = [('isresolved', 0, int),
                       ('isagn', 0, int)]

    @compound('properMotionRa', 'properMotionDec',
              'parallax', 'radialVelocity')
    def get_override_astrometry(self):
        n_obj = len(self.column_by_name('uniqueId'))
        return np.zeros((4,n_obj), dtype=float)


class Dc2RefCatGalaxies(Dc2RefCatMixin, DESCQACatalogMixin,
                        AstrometryGalaxies, EBVmixin,
                        InstanceCatalog):

    default_columns = [('isresolved', 1, int),
                       ('properMotionRa', 0.0, float),
                       ('properMotionDec', 0.0, float),
                       ('parallax', 0.0, float),
                       ('radialVelocity', 0.0, float),
                       ('galacticRv', 3.1, float)]

    _agn_param_db = None

    def _load_agn_id(self):
        if self._agn_param_db is None:
            self._agn_id_set = set()
            return

        with sqlite3.connect(self._agn_param_db) as connection:
            cc = connection.cursor()
            query = 'SELECT galaxy_id FROM agn_params'
            results = cc.execute(query).fetchall()
            self._agn_id_set = set([rr[0] for rr in results])

    @cached
    def get_isagn(self):
        galaxy_id = self.column_by_name('galaxy_id')
        output = np.zeros(len(self.column_by_name('uniqueId')), dtype=int)
        if len(output) == 0:
            return output

        if not hasattr(self, '_agn_id_set'):
            self._load_agn_id()

        for i_g, gg in enumerate(galaxy_id):
            if gg in self._agn_id_set:
                output[i_g] = 1
        return output

    def _calculate_fluxes(self, sedname, magnorm, redshift,
                          internal_av, internal_rv,
                          galactic_av, galactic_rv):
        """
        Parameters
        ----------
        sedname: ndarray of strings
            Name of SED files

        magnorm: ndarray of floats
            normalizing magnitudes

        redshift: ndarray of floats
            redshift of sources

        internal_av: ndarray of floats
            A_v due to internal dust

        internal_rv: ndarray of floats
            R_v due to internal dust

        galactic_av: ndarray of floats
            A_v due to Milky Way dust

        galactic_rv: ndarray of floats
            R_v due to Milky Way dust

        Returns
        -------
        flux_list: ndarray of floats
            flux_list[i][j] will be the flux of the jth object
            in the ith band of LSST.  Fluxes are in Janskys (see
            docstring of lsst.sims.photUtils.Sed.calcFlux)
        """

        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
            self.sed_dir = getPackageDir('sims_sed_library')

        flux_list = np.zeros((len(sedname), 6), dtype=float)
        for i_obj in range(len(sedname)):
            sed_obj = Sed()
            sed_obj.readSED_flambda(os.path.join(self.sed_dir,
                                                 defaultSpecMap[sedname[i_obj]]))

            fnorm = getImsimFluxNorm(sed_obj, magnorm[i_obj])
            sed_obj.multiplyFluxNorm(fnorm)
            a_int, b_int = sed_obj.setupCCMab()
            sed_obj.addCCMDust(a_int, b_int, A_v=internal_av[i_obj],
                               R_v=internal_rv[i_obj])

            sed_obj.redshiftSED(redshift[i_obj], dimming=True)

            a_gal, b_gal = sed_obj.setupCCMab()
            sed_obj.addCCMDust(a_gal, b_gal, A_v=galactic_av[i_obj],
                               R_v=galactic_rv[i_obj])

            obj_flux_list = self.lsstBandpassDict.fluxListForSed(sed_obj)
            flux_list[i_obj] = obj_flux_list

        return flux_list.transpose()

    @compound('lsst_u', 'lsst_g', 'lsst_r',
              'lsst_i', 'lsst_z', 'lsst_y')
    def get_reference_photometry(self):


        # NOTE: these are unlensed magnitudes.
        # protoDC2 does not include lensed magnitudes.
        # These columns should be replaced with mag_[ugrizY]_lsst
        # when we have cosmoDC2

        return np.array([self.column_by_name('mag_true_u_lsst'),
                         self.column_by_name('mag_true_g_lsst'),
                         self.column_by_name('mag_true_r_lsst'),
                         self.column_by_name('mag_true_i_lsst'),
                         self.column_by_name('mag_true_z_lsst'),
                         self.column_by_name('mag_true_Y_lsst')])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ra', type=float, default=55.064,
                        help="Center RA of area in degrees "
                        "(default = 55.064)")
    parser.add_argument('--dec', type=float, default=-29.783,
                        help="Center Dec of area in degrees "
                        "(default = -29.783)")

    parser.add_argument('--fov', type=float, default=3.6,
                        help="Field of view radius in degrees "
                        "(default = 3.6)")

    parser.add_argument('--out_dir', type=str, default='.',
                        help="Directory where file will be made "
                        "(default = '.')")

    parser.add_argument("--out_name", type=str,
                        default="dc2_reference_catalog.txt",
                        help="The name of the reference catalog file "
                        "to be written (default='dc2_reference_catalog.txt')")

    parser.add_argument("--cat", type=str, default="protoDC2",
                        help="The GCR handle for the extragalactic "
                        "catalog to be used (default='protoDC2')")

    parser.add_argument('--agn_db', type=str,
                        help='Path to AGN parameters database '
                             '(defaults to location at NERSC)',
                        default=os.path.join('/global', 'projecta',
                                             'projectdirs',
                                             'lsst', 'groups',
                                             'SSim', 'DC2',
                                             'agn_db_mbh_7.0_m_i_30.0.sqlite'))

    args = parser.parse_args()

    try:
        cache_LSST_seds(wavelen_min=0.0,
                        wavelen_max=1600.0)
    except PermissionError:
        pass

    obs = ObservationMetaData(pointingRA=args.ra,
                              pointingDec=args.dec,
                              boundType='circle',
                              boundLength=args.fov)


    file_name = os.path.join(args.out_dir, args.out_name)

    with open(file_name, 'w') as out_file:
        out_file.write('# reference catalog generated from: %s' % args.cat)

    star_db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                      port=1433, driver='mssql+pymssql')

    cat = Dc2RefCatStars(star_db, obs_metadata=obs)

    cat.write_catalog(file_name, chunk_size=10000, write_mode='a')

    gal_db = DESCQAReferenceObject(yaml_file_name=args.cat)
    gal_db.field_ra = obs.pointingRA
    gal_db.field_dec = obs.pointingDec

    cat = Dc2RefCatGalaxies(gal_db, obs_metadata=obs)
    cat._agn_param_db = args.agn_db

    cat.write_catalog(file_name, chunk_size=10000,
                      write_header=False,
                      write_mode='a')
