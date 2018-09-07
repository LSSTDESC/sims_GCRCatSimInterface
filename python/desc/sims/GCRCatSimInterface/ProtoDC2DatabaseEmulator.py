"""
This script will define classes that enable CatSim to interface with GCR
"""
import numpy as np
import os

from desc.sims.GCRCatSimInterface import DESCQAObject
from desc.sims.GCRCatSimInterface import deg2rad_double, arcsec2rad
from desc.sims.GCRCatSimInterfae import agnDESCQAObject

from lsst.sims.catalogs.db import DBObject

__all__ = ["DESCQAObject_protoDC2",
           "bulgeDESCQAObject_protoDC2",
           "diskDESCQAObject_protoDC2",
           "knotsDESCQAObject_protoDC2",
           "agnDESCQAObject_protoDC2",
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
                                              'ra_true', 'dec_true')
        gc.add_modifier_on_derived_quantities('decJ2000', self._transform_dec,
                                              'ra_true', 'dec_true')

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


class agnDESCQAObject_protoDC2(agnDESCQAObject, DESCAObject_protoDC2):
    objectTypeId = 117
    objid = 'agn_descqa'
