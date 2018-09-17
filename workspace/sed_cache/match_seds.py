import os
import h5py
import numpy as np
import healpy

import GCRCatalogs
from GCR import GCRQuery

from desc.sims.GCRCatSimInterface import sed_filter_names_from_catalog
from desc.sims.GCRCatSimInterface import sed_from_galacticus_mags

import argparse

_healpix_list = [10201, 10327, 10328, 10329, 10450,
                 10451, 10452, 10453, 10570, 10571,
                 10572, 10686, 10687]

_healpix_list = [10451]

def do_fitting(cat, component, healpix):

    filter_data = sed_filter_names_from_catalog(cat)
    filter_names = filter_data[component]['filter_name']
    wav_min = filter_data[component]['wav_min']
    wav_width = filter_data[component]['wav_width']

    H0 = cat.cosmology.H0.value
    Om0 = cat.cosmology.Om0

    healpix_query = GCRQuery('healpix_pixel==%d' % healpix)

    qties = cat.get_quantities(list(filter_names) + ['redshift', 'galaxy_id'],
                               native_filters=[healpix_query])

    mag_array = np.array([-2.5*np.log10(qties[ff]) for ff in filter_names])

    (sed_names,
     mag_norms) = sed_from_galacticus_mags(mag_array,
                                           qties['redshift'],
                                           H0, Om0,
                                           wav_min, wav_width)

    return qties['galaxy_id'], sed_names, mag_norms

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--healpix', type=int, default=None)
    args = parser.parse_args()
    assert args.healpix is not None

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image_addon_knots')

    gid, disk_sed, disk_mag = do_fitting(cat, 'disk', args.healpix)

    out_dir = os.path.join(os.environ['SCRATCH'], 'sed_cache')
    disk_file = h5py.File(os.path.join(out_dir, 'disk_%d.h5' % args.healpix), 'w')
    disk_file.create_dataset('galaxy_id', data=gid)
    disk_file.create_dataset('sed_name',
             data=np.array([dd.encode('utf-8') for dd in disk_sed]))
    disk_file.create_dataset('mag_norm', data=disk_mag)
    disk_file.close()

    del gid
    del disk_sed
    del disk_mag

    gid, bulge_sed, bulge_mag = do_fitting(cat, 'bulge', args.healpix)
    bulge_file = h5py.File(os.path.join(out_dir, 'bulge_%d.h5' % args.healpix), 'w')
    bulge_file.create_dataset('galaxy_id', data=gid)
    bulge_file.create_dataset('sed_name',
              data=np.array([bb.encode('utf-8') for bb in bulge_sed]))
    bulge_file.create_dataset('mag_norm', data=bulge_mag)
    bulge_file.close()
