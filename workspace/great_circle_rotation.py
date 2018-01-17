import numpy as np

from lsst.sims.utils import sphericalFromCartesian
from lsst.sims.utils import cartesianFromSpherical
from lsst.sims.utils import rotationMatrixFromVectors

class FieldRotator(object):

    def __init__(self, ra0, dec0, ra1, dec1):
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

        self._to_self_bases = np.array([xx,
                                        yy,
                                        zz])

        self._out_of_self_bases =self._to_self_bases.transpose()

        # find the transformation necessary to preserve
        # the due north relationship at field center
        d_dec = 0.1
        north = cartesianFromSpherical(np.radians(ra0),
                                       np.radians(dec0+d_dec))

        north = np.dot(first_rotation, north)

        #print(np.degrees(sphericalFromCartesian(north)))

        north_true = cartesianFromSpherical(np.radians(ra1),
                                            np.radians(dec1+d_dec))

        north = np.dot(self._to_self_bases, north)
        north_true = np.dot(self._to_self_bases, north_true)
        north = np.array([north[1], north[2]])
        north /= np.sqrt((north**2).sum())
        north_true = np.array([north_true[1], north_true[2]])
        north_true /= np.sqrt((north_true**2).sum())

        c = north_true[0]*north[0]+north_true[1]*north[1]
        s = north[0]*north_true[1]-north[1]*north_true[0]

        nprime = np.array([c*north[0]-s*north[1],
                           s*north[0]+c*north[1]])

        yz_rotation = np.array([[1.0, 0.0, 0.0],
                                [0.0, c, -s],
                                [0.0, s, c]])

        second_rotation = np.dot(self._out_of_self_bases,
                                 np.dot(yz_rotation,
                                        self._to_self_bases))

        self._transformation = np.dot(second_rotation,
                                      first_rotation)

    def transform(self, ra, dec):
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


