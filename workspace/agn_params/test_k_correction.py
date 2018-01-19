import unittest

import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, BandpassDict
from agn_param_module import k_correction

class K_correction_test_case(unittest.TestCase):

    def test_k_correction(self):

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
                            0.05)


if __name__ == "__main__":
    unittest.main()
