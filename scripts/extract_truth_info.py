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


    mu = 1./((1. - kappa)**2 - (gamma1**2 + gamma2**2))
    sed_file = find_file_path(sed_name, sed_dirs)
    return unique_id, ra_phosim, dec_phosim, mu, \
        desc.imsim.SedWrapper(sed_file, mag_norm, redshift,
                              internal_av, internal_rv,
                              galactic_av, galactic_rv, bp_dict)


def photParams(instcat):
    return desc.imsim.photometricParameters(
        desc.imsim.metadata_from_file(instcat))


def extract_truth_info(instcat, sed_dirs):
    obj_ids = []
    gal_ids = []
    ras, decs = [], []
    fluxes = defaultdict(list)
    phot_params = photParams(instcat)
    with open(instcat) as fd:
        for line in fd:
            if not line.startswith('object'):
                continue
            obj_id, ra, dec, mu, sed = extract_sed(line, sed_dirs)
            try:
                gal_id = int(obj_id) >> 10
            except ValueError:
                gal_id = obj_id
            obj_ids.append(obj_id)
            gal_ids.append(gal_id)
            ras.append(ra)
            decs.append(dec)
            for band in 'ugrizy':
                fluxes[band].append(mu*sed.sed_obj.calcFlux(bp_dict[band]))

    data = {'obj_id': obj_ids, 'gal_id': gal_ids, 'ra': ras, 'dec': decs}
    for band in 'ugrizy':
        data[f'flux_{band}'] = fluxes[band]
    return pd.DataFrame(data=data)


def make_gal_truthcat(df):
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


def make_truth_cats(det_name, sed_dirs):
    instcat = f'instcats/{det_name}_instcat.txt'
    truthcat = f'my_truth_cats/{det_name}_truthcat.pickle'
    gal_truthcat = f'my_truth_cats/{det_name}_gal_truthcat.pickle'

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
    sed_dirs = (os.environ['SIMS_SED_LIBRARY_DIR'],
                '/home/instcats/Run2.2i/00479028')

    det_names = []
    with open('det_names.txt') as fd:
        for line in fd:
            det_name = 'R{}{}_S{}{}'.format(*[_ for _ in line if _.isdigit()])
            if os.path.isfile(f'instcats/{det_name}_instcat.txt'):
                det_names.append(det_name)

    processes = 7
    with multiprocessing.Pool(processes=processes) as pool:
        workers = []
        for det_name in det_names[2:]:
            workers.append(pool.apply_async(make_truth_cats,
                                            (det_name, sed_dirs)))
        pool.close()
        pool.join()
        [_.get() for _ in workers]
