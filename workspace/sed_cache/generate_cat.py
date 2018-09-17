import numpy as np
import os
import healpy
from lsst.sims.catalogs.decorators import cached, compound
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from desc.sims.GCRCatSimInterface import diskDESCQAObject, bulgeDESCQAObject
from desc.sims.GCRCatSimInterface import PhoSimDESCQA

class DummyCat(PhoSimDESCQA):
    column_outputs = ['galaxy_id', 'healpixel',
                      'sedFilename_fitted', 'magNorm_fitted',
                      'A_v_comp', 'R_v_comp',
                      'redshift',
                      'mag_true_u_lsst', 'mag_true_g_lsst',
                      'mag_true_r_lsst', 'mag_true_i_lsst',
                      'mag_true_z_lsst', 'mag_true_y_lsst']

    @compound('A_v_comp', 'R_v_comp')
    def get_componentDust(self):
        if 'hasDisk' in self._cannot_be_null:
            av_name = 'A_v_disk'
            rv_name = 'R_v_disk'
        elif 'hasBulge' in self._cannot_be_null:
            av_name = 'A_v_bulge'
            rv_name = 'R_v_bulge'

        return np.array([self.column_by_name(av_name),
                         self.column_by_name(rv_name)])

    @cached
    def get_healpixel(self):
        rr = self.column_by_name('ra')
        dd = self.column_by_name('dec')
        return healpy.ang2pix(32, rr, dd, nest=False, lonlat=True)

project_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
assert os.path.isdir(project_dir)

minion = os.path.join(project_dir, 'minion_1016_desc_dithered_v4.db')
assert os.path.isfile(minion)

out_dir = os.path.join(os.environ['SCRATCH'], 'test_fitting')
assert os.path.isdir(out_dir)

obs_gen = ObservationMetaDataGenerator(minion)
obs = obs_gen.getObservationMetaData(obsHistID=2188, boundLength=2.1)[0]

vv = np.array([np.cos(obs._pointingDec)*np.cos(obs._pointingRA),
               np.cos(obs._pointingDec)*np.sin(obs._pointingRA),
               np.sin(obs._pointingDec)])

healpix_list = healpy.query_disc(32, vv, obs._boundLength,
                                 inclusive=True, nest=False)

print('\nhealpix list')
print(healpix_list)
print('\n')

obs.boundLength = 0.2

x_cat = 'cosmoDC2_v1.0_image_addon_knots'

diskDB = diskDESCQAObject(x_cat)
diskCat = DummyCat(diskDB, obs_metadata=obs, cannot_be_null=['hasDisk'])
diskCat.write_catalog(os.path.join(out_dir, 'disks.txt'),
                      write_mode='w', write_header=False,
                      chunk_size=100000)

del diskCat
del diskDB

bulgeDB = bulgeDESCQAObject(x_cat)
bulgeCat = DummyCat(bulgeDB, obs_metadata=obs, cannot_be_null=['hasBulge'])
bulgeCat.write_catalog(os.path.join(out_dir, 'bulges.txt'), write_mode='w',
                       write_header=False, chunk_size=100000)


