import os
import glob
from collections import defaultdict
import multiprocessing
import numpy as np
import pandas as pd
from lsst.sims.photUtils import BandpassDict
import desc.imsim
from desc.imsim.imSim import find_file_path

bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

def extract_sed(object_line, sed_dirs, gamma2_sign=-1):
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

    sed_file = find_file_path(sed_name, sed_dirs)
    return unique_id, desc.imsim.SedWrapper(sed_file, mag_norm, redshift,
                                            internal_av, internal_rv,
                                            galactic_av, galactic_rv, bp_dict)


def photParams(instcat):
    return desc.imsim.photometricParameters(
        desc.imsim.metadata_from_file(instcat))


def extract_truth_info(instcat, sed_dirs):
    obj_ids = []
    gal_ids = []
    fluxes = defaultdict(list)
    phot_params = photParams(instcat)
    with open(instcat) as fd:
        for line in fd:
            if not line.startswith('object'):
                continue
            obj_id, sed = extract_sed(line, sed_dirs)
            try:
                gal_id = int(obj_id) >> 10
            except ValueError:
                gal_id = obj_id
            obj_ids.append(obj_id)
            gal_ids.append(gal_id)
            for band in 'ugrizy':
                fluxes[band].append(sed.sed_obj.calcFlux(bp_dict[band]))

    data = {'obj_id': obj_ids, 'gal_id': gal_ids}
    for band in 'ugrizy':
        data[f'flux_{band}'] = fluxes[band]
    return pd.DataFrame(data=data)


def make_gal_truthcat(df):
    gal_ids = set(df['gal_id'])
    fluxes = defaultdict(lambda: np.zeros(len('ugrizy')))
    for iloc in range(len(df)):
        row = df.iloc[iloc]
        for i, band in enumerate('ugrizy'):
            fluxes[row['gal_id']][i] += row[f'flux_{band}']
    gal_fluxes = np.array(list(fluxes.values())).transpose()
    data = {'gal_id': list(fluxes.keys())}
    for band, flux in zip('ugrizy', gal_fluxes):
        data[f'flux_{band}'] = flux
    return pd.DataFrame(data=data)


def make_truth_cats(det_name):
    sed_dirs = (os.environ['SIMS_SED_LIBRARY_DIR'], '00479028')

    instcat = f'instcats/{det_name}_instcat.txt'
    truthcat = f'{det_name}_truthcat.pickle'
    gal_truthcat = f'{det_name}_gal_truthcat.pickle'

    if not os.path.isfile(truthcat):
        print(f'generating {det_name} truthcat')
        df = extract_truth_info(instcat, sed_dirs)
        df.to_pickle(truthcat)
    else:
        df = pd.read_pickle(truthcat)

    if not os.path.isfile(gal_truthcat):
        print(f'generating {det_name} galaxy truthcat')
        gal_df = make_gal_truthcat(df)
        gal_df.to_pickle(gal_truthcat)
    else:
        gal_df = pd.read_pickle(gal_truthcat)

    return df, gal_df

if __name__ == '__main__':
    det_names = []
    with open('sensor_list.txt') as fd:
        for line in fd:
            det_names.append('R{}{}_S{}{}'.format(*[_ for _ in line
                                                    if _.isdigit()]))

    with multiprocessing.Pool(processes=2) as pool:
        workers = []
        for det_name in det_names[:2]:
            workers.append(pool.apply_async(make_truth_cats, (det_name,)))
        pool.close()
        pool.join()
        [_.get() for _ in workers]
