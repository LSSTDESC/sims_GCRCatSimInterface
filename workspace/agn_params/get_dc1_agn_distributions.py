"""
This script queries all of the AGN on fatboy in the protoDC2 field
of view and writes their variability parameters to a text file
expected by plot_params.py
"""

from lsst.sims.catalogs.db import DBObject
from lsst.sims.photUtils import BandpassDict, Sed, getImsimFluxNorm
import os
from lsst.utils import getPackageDir

from test_m_i import make_histogram
from test_m_i import plot_color_mesh
from agn_param_module import k_correction
from lsst.sims.photUtils import CosmologyObject
import numpy as np
import json

if __name__ == "__main__":


    sed_name = os.path.join(getPackageDir('sims_sed_library'), 'agnSED',
                            'agn.spec.gz')
    assert os.path.exists(sed_name)
    base_sed = Sed()
    base_sed.readSED_flambda(sed_name)

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    bp = bp_dict['i']

    z_grid = np.arange(0.0, 16.0, 0.01)
    k_grid = np.zeros(len(z_grid),dtype=float)

    for i_z, zz in enumerate(z_grid):
        ss = Sed(flambda=base_sed.flambda, wavelen=base_sed.wavelen)
        ss.redshiftSED(zz, dimming=True)
        k = k_correction(ss, bp, zz)
        k_grid[i_z] = k

    cosmo = CosmologyObject()


    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    query = 'SELECT magnorm_agn, redshift, varParamStr FROM '
    query += 'galaxy WHERE varParamStr IS NOT NULL '
    query += 'AND dec BETWEEN -2.5 AND 2.5 '
    query += 'AND (ra<2.5 OR ra>357.5)'

    dtype = np.dtype([('magnorm',float), ('redshift',float),
                       ('varParamStr', str, 400)])

    data_iter = db.get_arbitrary_chunk_iterator(query, dtype=dtype,
                                                chunk_size=10000)


    with open('data/dc1_agn_params.txt', 'w') as out_file:
        out_file.write('# z m_i M_i tau sfu sfg sfr sfi sfz sfy\n')
        for chunk in data_iter:
            DM = cosmo.distanceModulus(redshift=chunk['redshift'])
            k_corr = np.interp(chunk['redshift'], z_grid, k_grid)

            for i_row, agn in enumerate(chunk):
                ss = Sed(wavelen=base_sed.wavelen, flambda=base_sed.flambda)
                fnorm = getImsimFluxNorm(ss, agn['magnorm'])
                ss.multiplyFluxNorm(fnorm)
                ss.redshiftSED(agn['redshift'], dimming=True)
                mag = ss.calcMag(bp)
                abs_m_i = mag - DM[i_row] - k_corr[i_row]
                params = json.loads(agn['varParamStr'])['pars']
                out_file.write('%e %e %e %e %e %e %e %e %e %e\n' %
                (agn['redshift'], mag, abs_m_i,
                 params['agn_tau'],
                 params['agn_sfu'],
                 params['agn_sfg'],
                 params['agn_sfr'],
                 params['agn_sfi'],
                 params['agn_sfz'],
                 params['agn_sfy']))
