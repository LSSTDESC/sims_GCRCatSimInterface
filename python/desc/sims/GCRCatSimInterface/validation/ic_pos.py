import os
import numpy as np
import healpy as hp
import gzip
import argparse

import GCR
import GCRCatalogs

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

__all__ = ["validate_instance_catalog_positions"]

def validate_instance_catalog_positions(cat_dir, obsid, fov_deg):
    """
    Parameters
    ----------
    cat_dir is the parent dir of $obs

    obsid is the obsHistID of the pointing (an int)

    fov_deg is the radius of the field of view in degrees
    """
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

    inst_cat_dir = os.path.join(cat_dir, '%.8d' % obsid)
    if not os.path.isdir(inst_cat_dir):
        raise RuntimeError('\n%s\nis not a dir\n' % inst_cat_dir)

    # Make sure no sprinkled galaxies ended up in the knots catalog
    knots_gid = []
    knots_name = os.path.join(inst_cat_dir, 'knots_cat_%d.txt.gz' % obsid)
    if not os.path.isfile(knots_name):
        raise RuntimeError('\n%s\nis not a file\n' % knots_name)
    with gzip.open(knots_name, 'rb') as in_file:
        for line in in_file:
            params = line.strip().split(b' ')
            gid = int(params[1])//1024
            knots_gid.append(gid)
    knots_gid = np.array(knots_gid)
    if len(knots_gid) == 0:
        raise RuntimeError("There were no knots at all")

    sprinkled_knots = np.in1d(knots_gid, sprinkled_gid)
    if sprinkled_knots.any():
        raise RuntimeError("There were knots in sprinkled systems")

    opsim_db = os.path.join('/global', 'projecta', 'projectdirs',
                            'lsst', 'groups', 'SSim', 'DC2',
                             'minion_1016_desc_dithered_v4_sfd.db')

    if not os.path.isfile(opsim_db):
        raise RuntimeError('\n%s\nis not a file\n' % opsim_db)

    obs_gen = ObservationMetaDataGenerator(database=opsim_db)
    obs = obs_gen.getObservationMetaData(obsHistID=obsid,
                                         boundType='circle',
                                         boundLength=fov_deg)[0]

    ra = np.degrees(obs.OpsimMetaData['descDitheredRA'])
    dec = np.degrees(obs.OpsimMetaData['descDitheredDec'])

    obs.pointingRA = ra
    obs.pointingDec = dec

    boresite = np.array([np.cos(np.radians(dec))*np.cos(np.radians(ra)),
                         np.cos(np.radians(dec))*np.sin(np.radians(ra)),
                         np.sin(np.radians(dec))])
    healpix_list = hp.query_disc(32, boresite, np.radians(fov_deg),
                                 nest=False, inclusive=True)

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image_addon_knots')

    qty_names = ['galaxy_id', 'ra', 'dec',
                 'stellar_mass_bulge', 'stellar_mass_disk']

    final_q = {}
    for name in qty_names:
        final_q[name] = []

    final_q['disk_sed_filter'] = []
    final_q['bulge_sed_filter'] = []

    sed_qty_names = []
    for nn in cat.list_all_quantities():
        if nn.startswith('sed') and (nn.endswith('bulge') or nn.endswith('disk')):
            sed_qty_names.append(nn)

    for hh in healpix_list:
        healpix_query = GCR.GCRQuery('healpix_pixel==%d' % hh)
        cat_q = cat.get_quantities(qty_names,
                                   filters=[(lambda x: x<=29.0, 'mag_r_lsst')],
                                   native_filters=[healpix_query])

        # only select those galaxies in our field of view
        cos_fov = np.cos(np.radians(fov_deg))
        cat_ra_rad = np.radians(cat_q['ra'])
        cat_dec_rad = np.radians(cat_q['dec'])
        xyz_cat = np.array([np.cos(cat_dec_rad)*np.cos(cat_ra_rad),
                            np.cos(cat_dec_rad)*np.sin(cat_ra_rad),
                            np.sin(cat_dec_rad)]).transpose()

        dot = np.dot(xyz_cat, boresite)
        cat_valid = np.where(dot>=cos_fov)

        sed_filter = {}
        sed_filter['disk'] = None
        sed_filter['bulge'] = None
        for component in ('disk', 'bulge'):
            for nn in sed_qty_names:
                if component not in nn:
                    continue
                sed_data = cat.get_quantities([nn],
                                   filters=[(lambda x: x<=29.0, 'mag_r_lsst')],
                                   native_filters=[healpix_query])

                if sed_filter[component] is None:
                    sed_filter[component] = (sed_data[nn]>0.0)
                else:
                    sed_filter[component] &= (sed_data[nn]>0.0)

        for kk in cat_q.keys():
            if kk in final_q:
                final_q[kk].append(cat_q[kk][cat_valid])

        final_q['disk_sed_filter'].append(sed_filter['disk'][cat_valid])
        final_q['bulge_sed_filter'].append(sed_filter['bulge'][cat_valid])

    cat_q = {}
    list_of_keys = list(final_q.keys())
    for kk in list_of_keys:
        data = final_q.pop(kk)
        cat_q[kk] = np.concatenate(data)

    sorted_dex = np.argsort(cat_q['galaxy_id'])
    for kk in cat_q.keys():
        cat_q[kk] = cat_q[kk][sorted_dex]

    for component in ('bulge', 'disk'):
        file_name = '%s_gal_cat_%d.txt.gz' % (component, obsid)
        full_name = os.path.join(inst_cat_dir, file_name)
        if not os.path.isfile(full_name):
            raise RuntimeError('\n%s\nis not a file\n' % full_name)

        instcat_gid = []
        instcat_ra = []
        instcat_dec = []
        instcat_magnorm = []

        with gzip.open(full_name, 'rb') as in_file:
            for line in in_file:
                params = line.strip().split(b' ')
                gid = int(params[1])//1024
                instcat_gid.append(gid)
                instcat_ra.append(float(params[2]))
                instcat_dec.append(float(params[3]))
                instcat_magnorm.append(float(params[4]))

        instcat_gid = np.array(instcat_gid)
        instcat_ra = np.array(instcat_ra)
        instcat_dec = np.array(instcat_dec)
        instcat_magnorm = np.array(instcat_magnorm)

        if instcat_gid.max() > max_gcr_gid:
            raise RuntimeError("\nmax gid in\n%s\n%e\n" % (full_name,
                                                           instcat_gid.max()))

        sorted_dex = np.argsort(instcat_gid)
        instcat_gid = instcat_gid[sorted_dex]
        instcat_ra = instcat_ra[sorted_dex]
        instcat_dec = instcat_dec[sorted_dex]
        instcat_magnorm = instcat_magnorm[sorted_dex]

        has_component = cat_q['stellar_mass_%s' % component]>0.0
        has_component &= cat_q['%s_sed_filter' % component]

        # the sprinkler will delete disks from the InstanceCatalog
        if component == 'disk':
            has_component &= ~np.in1d(cat_q['galaxy_id'], sprinkled_gid)

        has_component = np.where(has_component)

        gcr_gid = cat_q['galaxy_id'][has_component]
        gcr_ra = cat_q['ra'][has_component]
        gcr_dec = cat_q['dec'][has_component]
        gcr_comp = cat_q['stellar_mass_%s' % component][has_component]

        if not np.array_equal(gcr_gid, instcat_gid):
            gcr_in_inst = np.in1d(gcr_gid, instcat_gid)
            inst_in_gcr = np.in1d(instcat_gid, gcr_gid)

            if not gcr_in_inst.all():
                violation = ~gcr_in_inst
                n_violation = len(np.where(violation)[0])
                msg = '\nNot all GCR galaxies in InstanceCatalog\n'
                msg += '%d GCR galaxies were not in InstCat\n' % n_violation
                msg += 'd_len %d\n' % (len(instcat_gid)-len(gcr_gid))

                raise RuntimeError(msg)

            if not inst_in_gcr.all():
                violation = ~inst_in_gcr
                if instcat_magnorm[violation].min()<50.0:
                    n_violation = len(np.where(violation)[0])

                    msg = '\nNot all InstanceCatalog galaxies in GCR\n'
                    msg += 'WARNING: %d InstCat galaxies not in GCR' % n_violation
                    msg += '%s\n' % str(instcat_gid[violation])

                    raise RuntimeError(msg)

                # if this effects a change, it is because there were
                # some galaxies with magNorm>50 in the InstanceCatalog;
                # this is going to trim them out
                instcat_ra = instcat_ra[inst_in_gcr]
                instcat_dec = instcat_dec[inst_in_gcr]

        # now that we have verified all of the galaxies that should be in
        # the catalog are in the catalog, we will verify their positions by
        # converting all of the galaxy positions into Cartesian vectors,
        # choosing three points to define the orientation, and computing
        # the dot product of the whole set of Cartesian vectors with
        # those three.  These dot products should not differ between
        # the GCR objects and the InstanceCatalog objects

        gcr_ra_rad = np.radians(gcr_ra)
        gcr_dec_rad = np.radians(gcr_dec)
        gcr_cos_dec = np.cos(gcr_dec_rad)
        gcr_xyz = np.array([gcr_cos_dec*np.cos(gcr_ra_rad),
                            gcr_cos_dec*np.sin(gcr_ra_rad),
                            np.sin(gcr_dec_rad)]).transpose()

        instcat_ra_rad = np.radians(instcat_ra)
        instcat_dec_rad = np.radians(instcat_dec)
        instcat_cos_dec = np.cos(instcat_dec_rad)
        instcat_xyz = np.array([instcat_cos_dec*np.cos(instcat_ra_rad),
                                instcat_cos_dec*np.sin(instcat_ra_rad),
                                np.sin(instcat_dec_rad)]).transpose()

        rng = np.random.RandomState(77123)
        dx = rng.choice(np.arange(len(instcat_xyz), dtype=int), replace=False, size=3)

        for idx in dx:
            gcr_dot = np.dot(gcr_xyz, gcr_xyz[idx])
            instcat_dot = np.dot(instcat_xyz, instcat_xyz[idx])
            delta = np.abs(gcr_dot-instcat_dot)
            if delta.max() > 1.0e-6:
                raise RuntimeError("dot products were off %e %e %e" %
                                   (delta.min(), np.median(delta), delta.max()))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--fov_deg', type=float, default=2.1,
                        help='radius of field of view in degrees '
                        '(default: 2.1)')
    parser.add_argument('--obs', type=int, help='obsHistID')
    parser.add_argument('--cat_dir', type=str,
                        help='parent directory of $obsHistID/')

    args = parser.parse_args()

    validate_instance_catalog_positions(args.cat_dir, args.obs, args.fov_deg)
