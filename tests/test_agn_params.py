import unittest

import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, BandpassDict
from desc.sims.GCRCatSimInterface import k_correction
from desc.sims.GCRCatSimInterface import M_i_from_L_Mass

class M_i_test_case(unittest.TestCase):

    def test_m_i_abs(self):
        """
        Test the relationship between black hole mass,
        L/L_Eddington, and M_i by comparing to data points
        read off from MacLeod et al. 2010 Fig 15
        """
        l_edd = [-0.5, -0.5, -0.5,
                  0.15, 0.15, 0.1,
                 -0.9, -0.9, -0.9, -0.9,
                 -1.4, -1.4, -1.4,
                 -1.8, -1.8, -1.8]
        mbh = [7.7, 9.0, 9.6,
               8.4, 8.2, 7.9,
               9.6, 8.9, 8.4, 10.1,
               10.1, 9.6, 9.1,
               8.9, 9.2, 9.4]
        m_i = [-23.2, -26.5, -27.8,
               -26.8, -26.4, -26.0,
               -26.8, -25.2, -23.6, -28.5,
               -26.9, -25.2, -24.0,
               -23.2, -23.6, -24.0]

        l_edd = np.array(l_edd)
        mbh = np.array(mbh)
        m_i_test = M_i_from_L_Mass(l_edd, mbh)
        err = 0.0
        for ii in range(len(l_edd)):
            err += (m_i[ii]-m_i_test[ii])**2
        rms_err = np.sqrt(err/len(m_i_test))
        self.assertLess(rms_err, 0.5)


class K_correction_test_case(unittest.TestCase):

    def test_k_correction(self):
        """
        Test that the K correction correctly converts absolute magnitude
        to observed magnitude.
        """
        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        rng = np.random.RandomState(41321)
        sed_dir = os.path.join(getPackageDir('sims_sed_library'), 'galaxySED')
        list_of_sed_files = os.listdir(sed_dir)
        list_of_sed_files.sort()
        sed_to_check = rng.choice(list_of_sed_files, size=10)
        redshift_arr = rng.random_sample(len(sed_to_check))*2.0+0.1

        bp = bp_dict['g']
        for sed_name, zz in zip(sed_to_check, redshift_arr):
            full_name = os.path.join(sed_dir, sed_name)
            ss = Sed()
            ss.readSED_flambda(full_name)
            true_rest_mag = ss.calcMag(bp)
            ss.redshiftSED(zz, dimming=True)
            obs_mag = ss.calcMag(bp)
            k_corr = k_correction(ss, bp, zz)
            self.assertLess(np.abs(true_rest_mag-obs_mag+k_corr),
                            0.001)


if __name__ == "__main__":
    unittest.main()
