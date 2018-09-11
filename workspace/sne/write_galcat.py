import sys
# sys.path.insert(0, "/global/homes/k/kovacs/gcr-catalogs_cosmoDC2_v1.0/")
sys.path.insert(0, "/global/homes/r/rbiswas/src/DC2/gcr-catalogs")

import GCRCatalogs

## check version
print('GCRCatalogs =', GCRCatalogs.__version__, '|' ,'GCR =', GCRCatalogs.GCR.__version__)
for key in GCRCatalogs.available_catalogs:
    if 'cosmo' in key:
        print(key)

import pandas as pd


def write_hdf(healpixel):
    gc = GCRCatalogs.load_catalog('cosmoDC2_v0.4_test', 
                              {'healpix_pixels': [healpixel]})

    x = gc.get_quantities(['galaxy_id', 'stellar_mass', 'stellar_mass_disk', 'stellar_mass_bulge',
                       'redshift_true', 'size_disk_true', 'size_minor_disk_true', 'size_minor_bulge_true',
                       'size_bulge_true', 'ra', 'dec', 'position_angle_true'], return_iterator=True)
    df = pd.DataFrame(next(x))
    df.to_hdf('/global/cscratch1/sd/rbiswas/gals_{}_ra_dec.hdf'.format(healpixel),index=False, key='{}'.format(healpixel))

