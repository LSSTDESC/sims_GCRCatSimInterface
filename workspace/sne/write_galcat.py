import sys
import GCRCatalogs
import os

## check version
#print('GCRCatalogs =', GCRCatalogs.__version__, '|' ,'GCR =', GCRCatalogs.GCR.__version__)
#for key in GCRCatalogs.available_catalogs:
#    if 'cosmo' in key:
#        print(key)

import pandas as pd


def write_hdf(healpixel, out_dir):
    gc = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image',
                              {'healpix_pixels': [healpixel]})

    x = gc.get_quantities(['galaxy_id', 'stellar_mass', 'stellar_mass_disk', 'stellar_mass_bulge',
                       'redshift_true', 'size_disk_true', 'size_minor_disk_true', 'size_minor_bulge_true',
                       'size_bulge_true', 'ra', 'dec', 'position_angle_true'], return_iterator=True)
    df = pd.DataFrame(next(x))
    df.to_hdf(os.path.join(out_dir,'gals_%d_ra_dec.hdf' % healpixel),index=False, key='0')


if __name__ == "__main__":
    cat_config = GCRCatalogs.get_catalog_config('cosmoDC2_v1.0_image')
    healpix_list = cat_config['healpix_pixels']
    out_dir = os.path.join(os.environ['SCRATCH'], 'sne_production')
    assert os.path.isdir(out_dir)
    for hh in healpix_list:
        write_hdf(hh, out_dir)
