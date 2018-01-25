"""
This script will define classes that enable CatSim to interface with GCR
"""
import numpy as np
import os

from desc.sims.GCRCatSimInterface import DESCQAObject
from desc.sims.GCRCatSimInterface import deg2rad_double, arcsec2rad

from lsst.sims.catalogs.db import DBObject

__all__ = ["DESCQAObject_protoDC2",
           "bulgeDESCQAObject_protoDC2",
           "diskDESCQAObject_protoDC2",
           "agnDESCQAObject_protoDC2",
           "FieldRotator"]


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

    def _transform_catalog(self, gc):
        """
        Accept a GCR catalog object and add transformations to the
        columns in order to get the quantities expected by the CatSim
        code

        Parameters
        ----------
        gc -- a GCRCatalog object;
              the result of calling GCRCatalogs.load_catalog()
        """

        gc.add_modifier_on_derived_quantities('raJ2000', self._transform_ra,
                                              'ra_true', 'dec_true')
        gc.add_modifier_on_derived_quantities('decJ2000', self._transform_dec,
                                              'ra_true', 'dec_true')

        gc.add_quantity_modifier('redshift', gc.get_quantity_modifier('redshift_true'), overwrite=True)
        gc.add_quantity_modifier('true_redshift', gc.get_quantity_modifier('redshift_true'))
        gc.add_quantity_modifier('gamma1', gc.get_quantity_modifier('shear_1'))
        gc.add_quantity_modifier('gamma2', gc.get_quantity_modifier('shear_2'))
        gc.add_quantity_modifier('kappa', gc.get_quantity_modifier('convergence'))

        gc.add_quantity_modifier('positionAngle', gc.get_quantity_modifier('position_angle_true'))

        gc.add_modifier_on_derived_quantities('majorAxis::disk', arcsec2rad, 'size_disk_true')
        gc.add_modifier_on_derived_quantities('minorAxis::disk', arcsec2rad, 'size_minor_disk_true')
        gc.add_modifier_on_derived_quantities('majorAxis::bulge', arcsec2rad, 'size_bulge_true')
        gc.add_modifier_on_derived_quantities('minorAxis::bulge', arcsec2rad, 'size_minor_bulge_true')

        gc.add_quantity_modifier('sindex::disk', gc.get_quantity_modifier('sersic_disk'))
        gc.add_quantity_modifier('sindex::bulge', gc.get_quantity_modifier('sersic_bulge'))

        return None


class bulgeDESCQAObject_protoDC2(DESCQAObject_protoDC2):
    # PhoSim uniqueIds are generated by taking
    # source catalog uniqueIds, multiplying by
    # 1024, and adding objectTypeId.  This
    # components of the same galaxy to have
    # different uniqueIds, even though they
    # share a uniqueId in the source catalog
    objectTypeId = 97

    # some column names require an additional postfix
    _postfix = '::bulge'


class diskDESCQAObject_protoDC2(DESCQAObject_protoDC2):
    objectTypeId = 107
    _postfix = '::disk'


class agnDESCQAObject_protoDC2(DESCQAObject_protoDC2):
    objectTypeId = 117
    _columns_need_postfix = False

    descqaDefaultValues = {'varParamStr': (None, (str, 500)),
                           'magNorm': (np.NaN, np.float)}

    agn_params_db = None

    def _postprocess_results(self, master_chunk):
        """
        query the database specified by agn_params_db to
        find the AGN varParamStr associated with each AGN
        """

        if self.agn_params_db is None:
            return(master_chunk)

        if not os.path.exists(self.agn_params_db):
            raise RuntimeError('\n%s\n\ndoes not exist' % self.agn_params_db)

        if not hasattr(self, '_agn_dbo'):
            self._agn_dbo = DBObject(database=self.agn_params_db,
                                     driver='sqlite')

            self._agn_dtype = np.dtype([('galaxy_id', int),
                                        ('magNorm', float),
                                        ('varParamStr', str, 500)])

        gid_min = master_chunk['galaxy_id'].min()
        gid_max = master_chunk['galaxy_id'].max()

        query = 'SELECT galaxy_id, magNorm, varParamStr '
        query += 'FROM agn_params '
        query += 'WHERE galaxy_id BETWEEN %d AND %d ' % (gid_min, gid_max)
        query += 'ORDER BY galaxy_id'

        agn_data_iter = self._agn_dbo.get_arbitrary_chunk_iterator(query,
                                                        dtype=self._agn_dtype,
                                                        chunk_size=1000000)


        m_sorted_dex = np.argsort(master_chunk['galaxy_id'])
        m_sorted_id = master_chunk['galaxy_id'][m_sorted_dex]
        for agn_chunk in agn_data_iter:

            # find the indices of the elements in master_chunk
            # that correspond to elements in agn_chunk
            m_elements = np.in1d(m_sorted_id, agn_chunk['galaxy_id'])
            m_dex = m_sorted_dex[m_elements]

            # find the indices of the elements in agn_chunk
            # that correspond to elements in master_chunk
            a_dex = np.in1d(agn_chunk['galaxy_id'], m_sorted_id)

            # make sure we have matched elements correctly
            np.testing.assert_array_equal(agn_chunk['galaxy_id'][a_dex],
                                          master_chunk['galaxy_id'][m_dex])

            master_chunk['varParamStr'][m_dex] = agn_chunk['varParamStr'][a_dex]
            master_chunk['magNorm'][m_dex] = agn_chunk['magNorm'][a_dex]

        return master_chunk
