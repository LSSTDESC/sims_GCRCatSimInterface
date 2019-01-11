import os
import numpy as np
import healpy as hp
import gzip
import argparse

import GCR
import GCRCatalogs

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--fov_deg', type=float, default=2.1)
    parser.add_argument('--obs', type=int, help='obsHistID')
    parser.add_argument('--cat_dir', type=str,
                        help='dir containing obsHistID/ dir')

    args = parser.parse_args()

    agn_cache_file = os.path.join(os.environ['TWINKLES_DIR'], 'data',
                                  'cosmoDC2_v1.1.4_agn_cache.csv')

    if not os.path.isfile(agn_cache_file):
        raise RuntimeError('\n%s\nis not a file\n' % agn_cache_file)

    sne_cache_file = os.path.join(os.environ['TWINKLES_DIR'], 'data',
                                  'cosmoDC2_v1.1.4_sne_cache.csv')

    if not os.path.isfile(sne_cache_file):
        raise RuntimeError('\n%s\nis not a file\n' % sne_cache_file)

    # get a list of all of the galaxies that are going to be replace
    # by the sprinkler
    sprinkled_gid = []
    for file_name in (agn_cache_file, sne_cache_file):
        with open(file_name, 'r') as in_file:
            for line in in_file:
                if line.startswith('galtileid'):
                    continue
                params = line.strip().split(',')
                gid = int(params[0])
                sprinkled_gid.append(gid)

    sprinkled_gid = np.array(sprinkled_gid)

    max_gcr_gid = int(1.5e10) # this should actually be 1.2e10;
                              # the sprinkler will not put anything
                              # that close to the cutoff, though

    inst_cat_dir = os.path.join(args.cat_dir, '%.8d' % args.obs)
    if not os.path.isdir(inst_cat_dir):
        raise RuntimeError('\n%s\nis not a dir\n' % inst_cat_dir)

    opsim_db = os.path.join('/global', 'projecta', 'projectdirs',
                            'lsst', 'groups', 'SSim', 'DC2',
                             'minion_1016_desc_dithered_v4_sfd.db')

    if not os.path.isfile(opsim_db):
        raise RuntimeError('\n%s\nis not a file\n' % opsim_db)

    obs_gen = ObservationMetaDataGenerator(database=opsim_db)
    obs = obs_gen.getObservationMetaData(obsHistID=args.obs,
                                         boundType='circle',
                                         boundLength=args.fov_deg)[0]

    ra = np.degrees(obs.OpsimMetaData['descDitheredRA'])
    dec = np.degrees(obs.OpsimMetaData['descDitheredDec'])

    obs.pointingRA = ra
    obs.pointingDec = dec

    boresite = np.array([np.cos(np.radians(dec))*np.cos(np.radians(ra)),
                         np.cos(np.radians(dec))*np.sin(np.radians(ra)),
                         np.sin(np.radians(dec))])
    healpix_list = hp.query_disc(32, boresite, np.radians(args.fov_deg),
                                 nest=False, inclusive=True)

    healpix_query = GCR.GCRQuery('healpix_pixel==%d' % healpix_list[0])
    for hh in healpix_list[1:]:
        healpix_query |= GCR.GCRQuery('healpix_pixel==%d' % hh)

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image_addon_knots')
    cat_q = cat.get_quantities(['galaxy_id', 'ra', 'dec',
                                'stellar_mass_bulge', 'stellar_mass_disk',
                                'sed_5467_339_disk', 'sed_5467_339_bulge'],
                               filters=[(lambda x: x<=29.0, 'mag_r_lsst')],
                               native_filters=[healpix_query])

    # only select those galaxies in our field of view
    cos_fov = np.cos(np.radians(args.fov_deg))
    cat_ra_rad = np.radians(cat_q['ra'])
    cat_dec_rad = np.radians(cat_q['dec'])
    xyz_cat = np.array([np.cos(cat_dec_rad)*np.cos(cat_ra_rad),
                        np.cos(cat_dec_rad)*np.sin(cat_ra_rad),
                        np.sin(cat_dec_rad)]).transpose()

    dot = np.dot(xyz_cat, boresite)
    cat_valid = np.where(dot>=cos_fov)

    for kk in cat_q.keys():
        cat_q[kk] = cat_q[kk][cat_valid]

    sorted_dex = np.argsort(cat_q['galaxy_id'])
    for kk in cat_q.keys():
        cat_q[kk] = cat_q[kk][sorted_dex]

    for component in ('bulge', 'disk'):
        file_name = '%s_gal_cat_%d.txt.gz' % (component, args.obs)
        full_name = os.path.join(inst_cat_dir, file_name)
        if not os.path.isfile(full_name):
            raise RuntimeError('\n%s\nis not a file\n' % full_name)

        instcat_gid = []
        instcat_ra = []
        instcat_dec = [] 
        print('opening %s' % file_name)
        with gzip.open(full_name, 'rb') as in_file:
            for line in in_file:
                params = line.strip().split(b' ')
                gid = int(params[1])//1024
                instcat_gid.append(gid)
                instcat_ra.append(float(params[2]))
                instcat_dec.append(float(params[3]))
        print('read it in')

        instcat_gid = np.array(instcat_gid)
        instcat_ra = np.array(instcat_ra)
        instcat_dec = np.array(instcat_dec)

        if instcat_gid.max() > max_gcr_gid:
            raise RuntimeError("\nmax gid in\n%s\n%e\n" % (full_name,
                                                           instcat_gid.max()))

        sorted_dex = np.argsort(instcat_gid)
        instcat_gid = instcat_gid[sorted_dex]
        instcat_ra = instcat_ra[sorted_dex]
        instcat_dec = instcat_dec[sorted_dex]

        component_name = 'stellar_mass_%s' % component
        component_sed_name = 'sed_5467_339_%s' % component
        has_component = np.where(np.logical_and(cat_q[component_name]>0.0,
                                                cat_q[component_sed_name]>0.0))
        gcr_gid = cat_q['galaxy_id'][has_component]
        gcr_ra = cat_q['ra'][has_component]
        gcr_dec = cat_q['dec'][has_component]
        gcr_comp = cat_q[component_name][has_component]

        if not np.array_equal(gcr_gid, instcat_gid):
            gcr_in_inst = np.in1d(gcr_gid, instcat_gid, assume_unique=True)
            inst_in_gcr = np.in1d(instcat_gid, gcr_gid, assume_unique=True)

            if not gcr_in_inst.all():
                violation = ~gcr_in_inst
                missing_gid = gcr_gid[violation]
                are_sprinkled = np.in1d(missing_gid, sprinkled_gid)
                if not are_sprinkled.all() or component != 'disk':
                    knots_name = os.path.join(inst_cat_dir,
                                              'knots_cat_%d.txt.gz' % args.obs)
                    assert os.path.isfile(knots_name)
                    knots_gid = []
                    with gzip.open(knots_name, 'rb') as in_file:
                        for line in in_file:
                            params = line.strip().split(b' ')
                            knots_gid.append(int(params[1])//1024)
                    knots_gid = np.array(knots_gid)
                    n_violation = len(np.where(violation)[0])
                    vio_in_knots = np.in1d(gcr_gid[violation], knots_gid)
                    print('WARNING: %d GCR galaxies were not in InstCat'
                          % n_violation)
                    print('d_len %d' % (len(instcat_gid)-len(gcr_gid)))
                    print('are sprinkled: %s (%d)' % (str(are_sprinkled.all()),
                                                      len(np.where(are_sprinkled)[0])))
                    print('in knots: ', vio_in_knots.all())
                    print('n: %d (%d)' % (len(np.where(vio_in_knots)[0]), len(knots_gid)))
                    print('\n')
            if not inst_in_gcr.all():
                violation = ~inst_in_gcr
                n_violation = len(np.where(violation)[0])
                print('WARNING: %d InstCat galaxies not in GCR'
                      % n_violation)
                print('\n')
