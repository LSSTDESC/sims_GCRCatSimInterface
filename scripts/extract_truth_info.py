import os
import glob
import sys
from collections import defaultdict
import multiprocessing
import numpy as np
import pandas as pd
from lsst.sims.photUtils import BandpassDict
import desc.imsim
from desc.imsim.imSim import find_file_path
from get_config import get_config

bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

try:
    config_file = sys.argv[1]
except IndexError:
    config_file = 'instcat_validation_config.ini'
config = get_config(config_file)

def get_image_dirs():
    """Dummy function needed for FitsImage entries"""
    pass

def extract_sed(object_line, sed_dirs, gamma2_sign=-1,
                disable_MW=False):
    params = object_line.strip().split()
    unique_id = params[1]
    ra_phosim = float(params[2])
    dec_phosim = float(params[3])
    mag_norm = float(params[4])
    sed_name = params[5]
    redshift = float(params[6])
    gamma1 = float(params[7])
    gamma2 = gamma2_sign*float(params[8])
    kappa = float(params[9])
    internal_av = 0
    internal_rv = 0
    galactic_av = 0
    galactic_rv = 0
    fits_file = None
    pixel_scale = 0
    rotation_angle = 0
    npoints = 0
    semi_major_arcsec = 0
    semi_minor_arcsec = 0
    position_angle_degrees = 0
    sersic_index = 0
    if params[12].lower() == 'point':
        gs_type = 'pointSource'
        i_gal_dust_model = 14
        if params[13].lower() != 'none':
            i_gal_dust_model = 16
            internal_av = float(params[14])
            internal_rv =float(params[15])
        if params[i_gal_dust_model].lower() != 'none':
            galactic_av = float(params[i_gal_dust_model+1])
            galactic_rv = float(params[i_gal_dust_model+2])
    elif params[12].lower() == 'sersic2d':
        gs_type = 'sersic'
        semi_major_arcsec = float(params[13])
        semi_minor_arcsec = float(params[14])
        position_angle_degrees = float(params[15])
        sersic_index = float(params[16])
        i_gal_dust_model = 18
        if params[17].lower() != 'none':
            i_gal_dust_model = 20
            internal_av = float(params[18])
            internal_rv = float(params[19])
        if params[i_gal_dust_model].lower() != 'none':
            galactic_av = float(params[i_gal_dust_model+1])
            galactic_rv = float(params[i_gal_dust_model+2])
    elif params[12].lower() == 'knots':
        gs_type = 'RandomWalk'
        semi_major_arcsec = float(params[13])
        semi_minor_arcsec = float(params[14])
        position_angle_degrees = float(params[15])
        npoints = int(params[16])
        i_gal_dust_model = 18
        if params[17].lower() != 'none':
            i_gal_dust_model = 20
            internal_av = float(params[18])
            internal_rv = float(params[19])
        if params[i_gal_dust_model].lower() != 'none':
            galactic_av = float(params[i_gal_dust_model+1])
            galactic_rv = float(params[i_gal_dust_model+2])
    elif (params[12].endswith('.fits') or params[12].endswith('.fits.gz')):
        gs_type = 'FitsImage'
        fits_file = find_file_path(params[12], get_image_dirs())
        pixel_scale = float(params[13])
        rotation_angle = float(params[14])
        i_gal_dust_model = 16
        if params[15].lower() != 'none':
            i_gal_dust_model = 18
            internal_av = float(params[16])
            internal_rv = float(params[17])
        if params[i_gal_dust_model].lower() != 'none':
            galactic_av = float(params[i_gal_dust_model+1])
            galactic_rv = float(params[i_gal_dust_model+2])

    mu = 1./((1. - kappa)**2 - (gamma1**2 + gamma2**2))
    sed_file = find_file_path(sed_name, sed_dirs)
    if disable_MW:
        galactic_av = 0
    return unique_id, ra_phosim, dec_phosim, mu, gs_type, \
        desc.imsim.SedWrapper(sed_file, mag_norm, redshift,
                              internal_av, internal_rv,
                              galactic_av, galactic_rv, bp_dict)


def extract_truth_info(instcat, sed_dirs):
    global config
    obj_ids = []
    gal_ids = []
    ras, decs = [], []
    gs_types = []
    fluxes = defaultdict(list)
    with open(instcat) as fd:
        for line in fd:
            if not line.startswith('object'):
                continue
            obj_id, ra, dec, mu, gs_type, sed \
                = extract_sed(line, sed_dirs, disable_MW=config['disable_MW'])
            try:
                gal_id = int(obj_id) >> 10
            except ValueError:
                gal_id = obj_id
            obj_ids.append(obj_id)
            gal_ids.append(gal_id)
            ras.append(ra)
            decs.append(dec)
            gs_types.append(gs_type)
            for band in 'ugrizy':
                fluxes[band].append(mu*sed.sed_obj.calcFlux(bp_dict[band]))

    data = {'obj_id': obj_ids, 'gal_id': gal_ids, 'ra': ras, 'dec': decs,
            'gs_type': gs_types}
    for band in 'ugrizy':
        data[f'flux_{band}'] = fluxes[band]
    return pd.DataFrame(data=data)


def make_gal_truthcat(df_arg, exclude_gals_w_knots=False):
    if exclude_gals_w_knots:
        df = df_arg.query('gs_type != "RandomWalk"')
    else:
        df = df_arg
    gal_ids = set(df['gal_id'])
    fluxes = defaultdict(lambda: np.zeros(len('ugrizy')))
    coords = defaultdict(lambda: np.zeros(2))
    for iloc in range(len(df)):
        row = df.iloc[iloc]
        gal_id = row['gal_id']
        for i, band in enumerate('ugrizy'):
            fluxes[gal_id][i] += row[f'flux_{band}']
        coords[gal_id][0] = row['ra']
        coords[gal_id][1] = row['dec']
    gal_fluxes = np.array(list(fluxes.values())).transpose()
    data = {'gal_id': list(fluxes.keys())}
    data['ra'] = [_[0] for _ in coords.values()]
    data['dec'] = [_[1] for _ in coords.values()]
    for band, flux in zip('ugrizy', gal_fluxes):
        data[f'flux_{band}'] = flux
    return pd.DataFrame(data=data)


def make_truth_cats(det_name, sed_dirs, visit_dir):
    global config
    instcat_dir = f'instcats_{visit_dir}'
    os.makedirs(instcat_dir, exist_ok=True)
    my_truthcat_dir = config['ic_truthcat_prefix'] + f'_{visit_dir}'
    os.makedirs(my_truthcat_dir, exist_ok=True)

    instcat = f'{instcat_dir}/{det_name}_instcat.txt'
    truthcat = f'{my_truthcat_dir}/{det_name}_truthcat.pickle'
    gal_truthcat = f'{my_truthcat_dir}/{det_name}_gal_truthcat.pickle'

    if not os.path.isfile(truthcat):
        print(f'generating {det_name} truthcat')
        df = extract_truth_info(instcat, sed_dirs)
        df.to_pickle(truthcat)
    else:
        df = pd.read_pickle(truthcat)

    if not os.path.isfile(gal_truthcat):
        print(f'generating {det_name} galaxy truthcat')
        gal_df = make_gal_truthcat(df, exclude_gals_w_knots=False)
        gal_df.to_pickle(gal_truthcat)
    else:
        gal_df = pd.read_pickle(gal_truthcat)

    return df, gal_df

if __name__ == '__main__':
    sed_dirs = [os.path.dirname(_) for _ in glob.glob('00*/phosim_cat*.txt')]
    sed_dirs.append(os.environ['SIMS_SED_LIBRARY_DIR'])

#    visits = dict(u=2338, g=159521, r=40325, i=479028, z=8005, y=5883)
    visits = dict(i=6824)

    processes = config['processes']
    bands = config['bands']

    for band in bands:
        visit = visits[band]
        visit_dir = f'{visit:08d}'
        det_names = []
        with open(f'sensor_list_{visit_dir}.txt', 'r') as fd:
            for line in fd:
                det_name = 'R{}{}_S{}{}'.format(*[_ for _ in line if _.isdigit()])
                if os.path.isfile(f'instcats_{visit_dir}/{det_name}_instcat.txt'):
                    det_names.append(det_name)
        print(det_names)
        with multiprocessing.Pool(processes=processes) as pool:
            workers = []
            for det_name in det_names:
                workers.append(pool.apply_async(make_truth_cats,
                                                (det_name, sed_dirs, visit_dir)))
            pool.close()
            pool.join()
            [_.get() for _ in workers]
