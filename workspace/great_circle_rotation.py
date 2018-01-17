import numpy as np

from lsst.sims.utils import sphericalFromCartesian
from lsst.sims.utils import cartesianFromSpherical
from lsst.sims.utils import rotationMatrixFromVectors

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
        c = c/np.sqrt(c*c+s*s)
        s = s/np.sqrt(c*c+s*s)

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
        xyz = cartesianFromSpherical(np.radians(ra), np.radians(dec))
        xyz = np.dot(self._transformation, xyz)
        ra_out, dec_out = sphericalFromCartesian(xyz)
        return np.degrees(ra_out), np.degrees(dec_out)


ra_c_0 = 0.0
dec_c_0 =10.0

ra_c_1 = 14.36
dec_c_1= -37.75

rotator = FieldRotator(ra_c_0, dec_c_0,
                       ra_c_1, dec_c_1)

ra, dec = rotator.transform(ra_c_0,dec_c_0)
#print(ra,dec)
#exit()
try:
    assert np.abs(ra-ra_c_1)<1.0e-10
    assert np.abs(dec-dec_c_1)<1.0e-10
except:
    print('%.12e %.12e' % (ra,ra_c_1))
    print('%.12e %.12e' % (dec,dec_c_1))
    raise

ra_test = np.arange(ra_c_0-2.5,
                    ra_c_0+2.5, 0.1)



for ra_orig in ra_test:
    ra, dec = rotator.transform(ra_orig, 0.0)
    ra_north, dec_north = rotator.transform(ra_orig, 1.0e-1)
    ra_south, dec_south = rotator.transform(ra_orig, -1.0e-1)
    print(ra_orig-ra_c_0, ra_north-ra, dec_north-dec,
          ra_south-ra, dec_south-dec)


