import numpy as np

from lsst.sims.utils import sphericalFromCartesian
from lsst.sims.utils import cartesianFromSpherical
from lsst.sims.utils import rotationMatrixFromVectors



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

ra, dec = rotator.transform(ra_test, np.zeros(len(ra_test),dtype=float))
ra_north, dec_north = rotator.transform(ra_test, 0.1*np.ones(len(ra_test)))
ra_south, dec_south = rotator.transform(ra_test, -0.1*np.ones(len(ra_test)))

assert len(ra) == len(ra_test)
assert len(dec) == len(ra)
assert len(ra_north) == len(ra_test)
assert len(dec_north) == len(ra_north)
assert len(ra_south) == len(ra_test)
assert len(dec_south) == len(ra_test)

for ii in range(len(ra_test)):
    print(ra_test[ii]-ra_c_0, ra_north[ii]-ra[ii],
          dec_north[ii]-dec[ii],
          ra_south[ii]-ra[ii], dec_south[ii]-dec[ii])

exit()

for ra_orig in ra_test:
    ra, dec = rotator.transform(ra_orig, 0.0)
    ra_north, dec_north = rotator.transform(ra_orig, 1.0e-1)
    ra_south, dec_south = rotator.transform(ra_orig, -1.0e-1)
    print(ra_orig-ra_c_0, ra_north-ra, dec_north-dec,
          ra_south-ra, dec_south-dec)
