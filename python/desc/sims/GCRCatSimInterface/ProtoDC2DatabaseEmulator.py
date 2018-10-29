"""
This script will define classes that enable CatSim to interface with GCR
"""
import numpy as np
import os

from desc.sims.GCRCatSimInterface import DESCQAObject
from desc.sims.GCRCatSimInterface import deg2rad_double, arcsec2rad

from lsst.sims.utils import angularSeparation
from lsst.sims.catalogs.db import DBObject
from lsst.sims.utils import halfSpaceFromRaDec

__all__ = ["DESCQAObject_protoDC2",
           "bulgeDESCQAObject_protoDC2",
           "diskDESCQAObject_protoDC2",
           "knotsDESCQAObject_protoDC2",
           "agnDESCQAObject_protoDC2",
           "AGN_postprocessing_mixin",
           "FieldRotator"]


_ALPHA_Q_ADD_ON_IS_AVAILABLE = True
try:
    from GCRCatalogs.alphaq_addon import AlphaQAddonCatalog
except ImportError:
    _ALPHA_Q_ADD_ON_IS_AVAILABLE = False



_LSST_IS_AVAILABLE = True
try:
    from lsst.sims.utils import rotationMatrixFromVectors
    from lsst.sims.utils import cartesianFromSpherical, sphericalFromCartesian
except ImportError:
    _LSST_IS_AVAILABLE = False



class FieldRotator(object):

    def __init__(self, ra0, dec0, ra1, dec1):
        """
        Parameters
        ----------
        ra0, dec0 are the coordinates of the original field
        center in degrees

        ra1, dec1 are the coordinates of the new field center
        in degrees

        The transform() method of this class operates by first
        applying a rotation that carries the original field center
        into the new field center.  Points are then transformed into
        a basis in which the unit vector defining the new field center
        is the x-axis.  A rotation about the x-axis is applied so that
        a point that was due north of the original field center is still
        due north of the field center at the new location.  Finally,
        points are transformed back into the original x,y,z bases.
        """

        # do we actually need to do the rotation, or is the simulation
        # already in the right spot?
        self._needs_to_be_rotated = True
        rot_dist = angularSeparation(ra0, dec0, ra1, dec1)
        if rot_dist<1.0/3600.0:
            self._needs_to_be_rotated = False
            return

        # find the rotation that carries the original field center
        # to the new field center
        xyz = cartesianFromSpherical(np.radians(ra0), np.radians(dec0))
        xyz1 = cartesianFromSpherical(np.radians(ra1), np.radians(dec1))
        if np.abs(1.0-np.dot(xyz, xyz1))<1.0e-10:
            self._transformation = np.identity(3, dtype=float)
            return

        first_rotation = rotationMatrixFromVectors(xyz, xyz1)

        # create a basis set in which the unit vector
        # defining the new field center is the x axis
        xx = np.dot(first_rotation, xyz)
        rng = np.random.RandomState(99)
        mag = np.NaN
        while np.abs(mag)<1.0e-20 or np.isnan(mag):
            random_vec = rng.random_sample(3)
            comp = np.dot(random_vec, xx)
            yy = random_vec - comp*xx
            mag = np.sqrt((yy**2).sum())
            yy /= mag

        zz = np.cross(xx, yy)

        to_self_bases = np.array([xx,
                                  yy,
                                  zz])

        out_of_self_bases =to_self_bases.transpose()

        # Take a point due north of the original field
        # center.  Apply first_rotation to carry it to
        # the new field.  Transform it to the [xx, yy, zz]
        # bases and find the rotation about xx that will
        # make it due north of the new field center.
        # Finally, transform back to the original bases.
        d_dec = 0.1
        north = cartesianFromSpherical(np.radians(ra0),
                                       np.radians(dec0+d_dec))

        north = np.dot(first_rotation, north)

        #print(np.degrees(sphericalFromCartesian(north)))

        north_true = cartesianFromSpherical(np.radians(ra1),
                                            np.radians(dec1+d_dec))

        north = np.dot(to_self_bases, north)
        north_true = np.dot(to_self_bases, north_true)
        north = np.array([north[1], north[2]])
        north /= np.sqrt((north**2).sum())
        north_true = np.array([north_true[1], north_true[2]])
        north_true /= np.sqrt((north_true**2).sum())

        c = north_true[0]*north[0]+north_true[1]*north[1]
        s = north[0]*north_true[1]-north[1]*north_true[0]
        norm = np.sqrt(c*c+s*s)
        c = c/norm
        s = s/norm

        nprime = np.array([c*north[0]-s*north[1],
                           s*north[0]+c*north[1]])

        yz_rotation = np.array([[1.0, 0.0, 0.0],
                                [0.0, c, -s],
                                [0.0, s, c]])

        second_rotation = np.dot(out_of_self_bases,
                                 np.dot(yz_rotation,
                                        to_self_bases))

        self._transformation = np.dot(second_rotation,
                                      first_rotation)

    def transform(self, ra, dec):
        """
        ra, dec are in degrees; return the RA, Dec coordinates
        of the point about the new field center
        """
        xyz = cartesianFromSpherical(np.radians(ra), np.radians(dec)).transpose()
        xyz = np.dot(self._transformation, xyz).transpose()
        ra_out, dec_out = sphericalFromCartesian(xyz)
        return np.degrees(ra_out), np.degrees(dec_out)


class DESCQAObject_protoDC2(DESCQAObject):
    """
    This class is meant to mimic the CatalogDBObject usually used to
    connect CatSim to a database.
    """

    _cat_cache_suffix = '_rotated'
    tableid = 'galaxy'
    database = 'LSSTCATSIM'
    yaml_file_name = 'protoDC2'

    def _rotate_to_correct_field(self, ra_rad, dec_rad):
        """
        Takes arrays of RA and Dec (in radians) centered
        on RA=0, Dec=0 and rotates them to coordinates
        centered on self.field_ra, self.field_dec

        Returns the rotated RA, Dec in radians.
        """

        if not _LSST_IS_AVAILABLE:
            raise RuntimeError("\nCannot use DESCQAObject_protoDC2\n"
                               "The LSST simulations stack is not setup\n")

        if not hasattr(self, '_rotate_ra_in_cache'):
            self._rotate_ra_in_cache = None
            self._rotate_dec_in_cache = None

        if not hasattr(self, '_field_rotator'):
            self._field_rotator = FieldRotator(0.0, 0.0, self.field_ra, self.field_dec)

        if not self._field_rotator._needs_to_be_rotated:
            return ra_rad, dec_rad

        if self._rotate_ra_in_cache is None or \
           not np.array_equal(ra_rad, self._rotate_ra_in_cache) or \
           not np.array_equal(dec_rad, self._rotate_dec_in_cache):

            self._rotate_ra_in_cache = ra_rad
            self._rotate_dec_in_cache = dec_rad

            ra, dec = self._field_rotator.transform(np.degrees(ra_rad),
                                                    np.degrees(dec_rad))

            self._ra_rotated = np.radians(ra)
            self._dec_rotated = np.radians(dec)

        return self._ra_rotated, self._dec_rotated

    def _transform_ra(self, ra_deg, dec_deg):
        """
        Transform RA in degrees to RA_rot in radians
        where RA was from a catalog centered on
        RA=0, Dec=0 and RA_rot is from a catalog
        centered on RA=self.field_ra, Dec=self.field_dec
        """
        ra, dec = self._rotate_to_correct_field(deg2rad_double(ra_deg),
                                                deg2rad_double(dec_deg))
        return ra

    def _transform_dec(self, ra_deg, dec_deg):
        """
        Transform Dec in degrees to Dec_rot in radians
        where Dec was from a catalog centered on
        RA=0, Dec=0 and Dec_rot is from a catalog
        centered on RA=self.field_ra, Dec=self.field_dec
        """

        ra, dec = self._rotate_to_correct_field(deg2rad_double(ra_deg),
                                                deg2rad_double(dec_deg))
        return dec

    def _transform_object_coords(self, gc):
        """
        Apply transformations to the RA, Dec of astrophysical sources

        gc is a GCR catalog instance
        """
        gc.add_modifier_on_derived_quantities('raJ2000', self._transform_ra,
                                              'ra', 'dec')
        gc.add_modifier_on_derived_quantities('decJ2000', self._transform_dec,
                                              'ra', 'dec')

class bulgeDESCQAObject_protoDC2(DESCQAObject_protoDC2):
    # PhoSim uniqueIds are generated by taking
    # source catalog uniqueIds, multiplying by
    # 1024, and adding objectTypeId.  This
    # components of the same galaxy to have
    # different uniqueIds, even though they
    # share a uniqueId in the source catalog
    objectTypeId = 97
    objid = 'bulge_descqa'

    # some column names require an additional postfix
    _postfix = '::bulge'


class diskDESCQAObject_protoDC2(DESCQAObject_protoDC2):
    objectTypeId = 107
    objid = 'disk_descqa'
    _postfix = '::disk'


class knotsDESCQAObject_protoDC2(DESCQAObject_protoDC2):
    objectTypeId = 127
    objid = 'knots_descqa'
    _postfix = '::knots'


class AGN_postprocessing_mixin(object):

    def _do_agn_query(self, half_space):
        """
        Actually query the AGN parameter database for all AGN
        inside of the field of view specified by half_space
        (which is an lsst.sims.utils.HalfSpace)
        """

        if not os.path.exists(self.agn_params_db):
            raise RuntimeError('\n%s\n\ndoes not exist' % self.agn_params_db)

        if not hasattr(self, '_agn_dbo'):
            self._agn_dbo = DBObject(database=self.agn_params_db,
                                     driver='sqlite')

            self._agn_dtype = np.dtype([('galaxy_id', int),
                                        ('magNorm', float),
                                        ('varParamStr', str, 500)])



        self._cached_half_space = half_space
        trixel_bounds = half_space.findAllTrixels(8)

        query = 'SELECT galaxy_id, magNorm, varParamStr '
        query += 'FROM agn_params '
        query += 'WHERE '
        for i_bound, bound in enumerate(trixel_bounds):
            if i_bound>0:
                query += 'OR '
            if bound[0]==bound[1]:
                query += 'htmid_8 == %d ' % bound[0]
            else:
                query += '(htmid_8 >= %d AND htmid_8 <= %d) ' % (bound[0], bound[1])
        query += 'ORDER BY galaxy_id'

        self._agn_query_results = self._agn_dbo.execute_arbitrary(query,
                                                                  dtype=self._agn_dtype)


    def _prefilter_galaxy_id(self, obs_metadata):
        """
        Accept an ObservationMetaData characterizing
        the current pointing.

        Return a numpy array of galaxy_ids that are in the
        field of view and actually contain an AGN.
        """

        half_space = halfSpaceFromRaDec(obs_metadata.pointingRA,
                                        obs_metadata.pointingDec,
                                        obs_metadata.boundLength)

        if (not hasattr(self, "_agn_query_results") or
            half_space != self._cached_half_space):

            self._do_agn_query(half_space)

        return np.sort(self._agn_query_results['galaxy_id'])




    def _postprocess_results(self, master_chunk, obs_metadata):
        """
        query the database specified by agn_params_db to
        find the AGN varParamStr associated with each AGN
        """

        if self.agn_objid is None:
            gid_name = 'galaxy_id'
            varpar_name = 'varParamStr'
            magnorm_name = 'magNorm'
        else:
            gid_name = self.agn_objid + '_' + 'galaxy_id'
            varpar_name = self.agn_objid + '_' + 'varParamStr'
            magnorm_name = self.agn_objid + '_' + 'magNorm'

        if self.agn_params_db is None:
            return(master_chunk)

        half_space = halfSpaceFromRaDec(obs_metadata.pointingRA,
                                        obs_metadata.pointingDec,
                                        obs_metadata.boundLength)

        if (not hasattr(self, "_agn_query_results") or
            half_space != self._cached_half_space):

            self._do_agn_query(half_space)

        gid_arr = master_chunk[gid_name]
        m_sorted_dex = np.argsort(gid_arr)
        m_sorted_id = gid_arr[m_sorted_dex]
        valid_agn_dex = np.where(np.logical_and(self._agn_query_results['galaxy_id']>=gid_arr.min(),
                                                self._agn_query_results['galaxy_id']<=gid_arr.max()))

        valid_agn = self._agn_query_results[valid_agn_dex]

        # find the indices of the elements in master_chunk
        # that correspond to elements in agn_chunk
        m_elements = np.in1d(m_sorted_id, valid_agn['galaxy_id'], assume_unique=True)
        m_dex = m_sorted_dex[m_elements]

        # find the indices of the elements in agn_chunk
        # that correspond to elements in master_chunk
        a_dex = np.in1d(valid_agn['galaxy_id'], m_sorted_id, assume_unique=True)

        # make sure we have matched elements correctly
        np.testing.assert_array_equal(valid_agn['galaxy_id'][a_dex],
                                      master_chunk[gid_name][m_dex])

        if varpar_name in master_chunk.dtype.names:
            master_chunk[varpar_name][m_dex] = valid_agn['varParamStr'][a_dex]

        if magnorm_name in master_chunk.dtype.names:
            master_chunk[magnorm_name][m_dex] = valid_agn['magNorm'][a_dex]

        return self._final_pass(master_chunk)


class agnDESCQAObject_protoDC2(AGN_postprocessing_mixin, DESCQAObject_protoDC2):
    objectTypeId = 117
    objid = 'agn_descqa'
    _columns_need_postfix = False

    descqaDefaultValues = {'is_sprinkled': (0, int),
                           'varParamStr': (None, (str, 500)),
                           'magNorm': (np.NaN, np.float),
                           'sedFilename': ('agn.spec', (str, 500)),
                           'sn_truth_params': (None, (str, 500)),
                           'sn_t0': (0.0, np.float)}

    agn_params_db = None
    agn_objid = None
