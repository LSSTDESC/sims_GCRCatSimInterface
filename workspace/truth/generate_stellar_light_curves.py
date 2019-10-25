import os
import sqlite3
import h5py
import numpy as np

from lsst.sims.utils import angularSeparation
from lsst.sims.utils import getRotSkyPos
from lsst.sims.utils import defaultSpecMap
import lsst.sims.photUtils as sims_photUtils
from lsst.sims.utils import ObservationMetaData
import lsst.sims.catUtils.mixins.VariabilityMixin as variability
from lsst.sims.coordUtils import chipNameFromRaDecLSST

from dc2_spatial_definition import DC2_bounds
from lsst.sims.utils import xyz_from_ra_dec

import multiprocessing
import time


def create_out_name(out_dir, hpid):
    out_file_name = os.path.join(out_dir, 'stellar_lc_%d.h5' % hpid)
    return out_file_name

class VariabilityGenerator(variability.StellarVariabilityModels,
                           variability.MLTflaringMixin,
                           variability.ParametrizedLightCurveMixin):
    """
    This is a class that mimics the API of an InstanceCatalog. This makes it
    easier for us to generate light curves from sources with differing
    variability models.
    """

    def __init__(self, chunk):
        """
        chunk = (simobjid, hpid, sedFilename, magNorm, ebv, varParamStr,
                 parallax, ra, dec)
        """

        # initialize some member variables expected by the
        # stellar variability models
        self.photParams = sims_photUtils.PhotometricParameters(nexp=1,
                                                              exptime=30)
        self.lsstBandpassDict = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()
        self._actually_calculated_columns = []
        for bp in 'ugrizy':
            self._actually_calculated_columns.append('lsst_%s' % bp)

        sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']
        self.quiescent_mags = {}
        for bp in 'ugrizy':
            self.quiescent_mags[bp] = np.zeros(len(chunk), dtype=float)

        self.parallax = np.zeros(len(chunk), dtype=float)
        self.ebv = np.zeros(len(chunk), dtype=float)
        self.varParamStr = np.empty(len(chunk), dtype=(str,300))
        self.simobjid = np.empty(len(chunk), dtype=int)
        ccm_wav = None

        # find the quiescent magnitudes of the stars
        for i_star, star in enumerate(chunk):
            self.simobjid[i_star] = int(star[0])
            full_file_name = os.path.join(sed_dir,
                                          defaultSpecMap[star[2]])
            spec = sims_photUtils.Sed()
            spec.readSED_flambda(full_file_name)
            fnorm = sims_photUtils.getImsimFluxNorm(spec, float(star[3]))
            spec.multiplyFluxNorm(fnorm)
            if ccm_wav is None or not np.array_equal(spec.wavelen, ccm_wav):
                ccm_wav = np.copy(spec.wavelen)
                a_x, b_x = spec.setupCCM_ab()
            ebv_val = float(star[4])
            spec.addDust(a_x, b_x, ebv=ebv_val, R_v=3.1)
            mag_list = self.lsstBandpassDict.magListForSed(spec)
            for i_bp, bp in enumerate('ugrizy'):
                self.quiescent_mags[bp][i_star] = mag_list[i_bp]
            self.ebv[i_star] = ebv_val
            self.varParamStr[i_star] = star[5]
            self.parallax[i_star] = float(star[6])

        self.parallax *= (np.pi/648000000.0) # convert from mas to radians

    def column_by_name(self, name):
        if name.startswith('quiescent'):
            bp = name[-1]
            return np.copy(self.quiescent_mags[bp])
        elif name == 'ebv':
            return np.copy(self.ebv)
        elif name=='parallax':
            return np.copy(self.parallax)
        elif name=='simobjid':
            return np.copy(self.simobjid)
        raise RuntimeError("\n\nCannot get column %s\n\n" % name)


def do_photometry(chunk,
                  obs_lock, obs_metadata_dict,
                  star_lock, star_data_dict,
                  job_lock, job_dict,
                  out_dir):
    """
    -Take a chunk from the stellar database
    -Calculate the light curves for all of the stars in the chunk
    -Write the light curves to disk
    -Save the metadata and static data to dicts to be output later

    Parameters
    ----------
    chunk -- a set of stars returned by sqlite3.fetchmany()
             the columns in chunk are:
             simobjid, hpid, sedFilename, magNorm, ebv, varParamStr,
             parallax, ra, decl
             Note: chunk must be from one hpid (healpix pixel ID)

    obs_lock -- the lock corresponding to the observation metadata
    for this healpixel

    obs_metadata_dict -- the dict where observation metadata is stored

    star_lock -- the lock corresponding to the stellar metadata for this
    healpixel

    star_data_dict -- the dict where stellar metadata (i.e. static data) is
    stored

    job_lock -- the lock used to prevent too many processes from entering
    light curve generation at once

    job_dict -- the dict keeping track of how many jobs are doing light curve
    generation

    out_dir -- directory where light curve files are to be stored

    Returns
    -------
    None
    """
    dummy_sed = sims_photUtils.Sed()

    # find the lookup file associating healpixel with obsHistID;
    # this should probably actually be passed into the method
    data_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth'
    assert os.path.isdir(data_dir)
    hpid_lookup_name = os.path.join(data_dir, 'hpid_to_obsHistID_lookup.h5')
    assert os.path.isfile(hpid_lookup_name)

    metadata_dict = {}
    metadata_keys = ['obsHistID', 'ra', 'dec', 'rotTelPos', 'mjd', 'filter']
    hpid = int(chunk[0][1])
    with h5py.File(hpid_lookup_name, 'r') as in_file:
        valid_obsid = in_file['%d' % hpid][()]
        for k in metadata_keys:
            metadata_dict[k] = in_file[k][()]

    valid_obsid = np.sort(valid_obsid)
    valid_dex = np.searchsorted(metadata_dict['obsHistID'], valid_obsid)
    np.testing.assert_array_equal(metadata_dict['obsHistID'][valid_dex],
                                  valid_obsid)

    for k in metadata_dict.keys():
        metadata_dict[k] = metadata_dict[k][valid_dex]

    # generate delta_magnitude light curves
    has_dmag = False
    while not has_dmag:

        # make sure no more than 5 processes are doing this at once
        with job_lock:
            if job_dict['running_dmag']>=5:
                continue
            else:
                job_dict['running_dmag'] += 1
                #print('running dmag %d' % job_dict['running_dmag'])

        t_start = time.time()
        var_gen = VariabilityGenerator(chunk)
        dmag_raw = var_gen.applyVariability(var_gen.varParamStr,
                                            expmjd=metadata_dict['mjd'])

        dmag_raw = dmag_raw.transpose([1,2,0])
        # transpose so the columns are (star, mjd, filter)
        assert dmag_raw.shape == (len(chunk), len(metadata_dict['mjd']), 6)

        dmag = np.zeros((len(chunk), len(metadata_dict['mjd'])),
                        dtype=float)

        for i_mjd in range(len(metadata_dict['mjd'])):
            dmag[:,i_mjd] = dmag_raw[:,i_mjd,metadata_dict['filter'][i_mjd]]

        del dmag_raw
        has_dmag = True
        with job_lock:
            job_dict['running_dmag'] -= 1
            #print('running dmag %d' % job_dict['running_dmag'])

    quiescent_fluxes = np.zeros((len(chunk),6), dtype=float)
    for i_bp, bp in enumerate('ugrizy'):
        quiescent_fluxes[:, i_bp] = dummy_sed.fluxFromMag(var_gen.column_by_name('quiescent_%s' % bp))

    t_dmag = time.time()-t_start

    star_ra = np.array([c[7] for c in chunk])
    star_dec = np.array([c[8] for c in chunk])

    # Now figure out which obsHistID are observing which stars
    obs_mask = np.zeros((len(chunk), len(metadata_dict['mjd'])),
                        dtype=bool)

    fov_radius = 2.1 # in degrees
    t_start = time.time()
    for i_obs in range(len(metadata_dict['mjd'])):
        ra = np.degrees(metadata_dict['ra'][i_obs])
        dec = np.degrees(metadata_dict['dec'][i_obs])

        ## back-of-the-envelope implementation
        ## (does not use focal plane model with chip gaps, etc.)
        #ang_dist = angularSeparation(ra, dec, star_ra, star_dec)
        #valid = ang_dist<fov_radius

        rotTel = np.degrees(metadata_dict['rotTelPos'][i_obs])
        obs_md = ObservationMetaData(pointingRA=ra, pointingDec=dec,
                                     mjd=metadata_dict['mjd'][i_obs])
        rotSky = getRotSkyPos(ra, dec, obs_md, rotTel)
        obs_md.rotSkyPos = rotSky
        chip_name = chipNameFromRaDecLSST(star_ra, star_dec,
                                          obs_metadata=obs_md).astype(str)
        valid = np.char.find(chip_name, 'None')<0

        obs_mask[:,i_obs] = valid

    # verify that any stars with empty light curves are, in fact
    # outside of DC2
    region_of_interest = DC2_bounds()
    xyz_offenders = []
    for i_star in range(len(chunk)):
        if obs_mask[i_star].sum()==0:
            xyz_offenders.append(xyz_from_ra_dec(star_ra[i_star],
                                                 star_dec[i_star]))

    if len(xyz_offenders)>0:
        xyz_offenders = np.array(xyz_offenders)
        v0 = region_of_interest.hs_list[0].contains_many_pts(xyz_offenders)
        v1 = region_of_interest.hs_list[1].contains_many_pts(xyz_offenders)
        v2 = region_of_interest.hs_list[2].contains_many_pts(xyz_offenders)
        v3 = region_of_interest.hs_list[3].contains_many_pts(xyz_offenders)
        if (v0&v1&v2&v3).sum()>0:
            msg = "\n\nsome stars in healpixel %d unobserved\n\n" % hpid
            raise RuntimeError(msg)

    t_flux = time.time()
    dflux_arr = []
    for i_star in range(len(chunk)):
        star_filter = metadata_dict['filter'][obs_mask[i_star]]
        q_flux = quiescent_fluxes[i_star][star_filter]
        assert len(q_flux) == obs_mask[i_star].sum()
        q_mag = dummy_sed.magFromFlux(q_flux)
        tot_mag = q_mag + dmag[i_star,obs_mask[i_star]]
        tot_flux = dummy_sed.fluxFromMag(tot_mag)
        dflux = tot_flux-q_flux
        dflux_arr.append(dflux)
        assert len(dflux) == obs_mask[i_star].sum()

    # store metadata in the output dicts
    with obs_lock:
        obs_metadata_dict['mjd'].append(metadata_dict['mjd'])
        obs_metadata_dict['obsHistID'].append(metadata_dict['obsHistID'])
        obs_metadata_dict['filter'].append(metadata_dict['filter'])

    with star_lock:
        local_simobjid = var_gen.column_by_name('simobjid')
        out_file_name = create_out_name(out_dir, hpid)

        # write the light curves to disk now; otherwise, the memory
        # footprint of the job becomes too large
        with h5py.File(out_file_name, 'a') as out_file:
            for i_star in range(len(local_simobjid)):
                if len(dflux_arr[i_star])==0:
                    continue
                simobjid = local_simobjid[i_star]
                out_file.create_dataset('%d_flux' % simobjid,
                                        data=dflux_arr[i_star])
                out_file.create_dataset('%d_obsHistID' % simobjid,
                              data=metadata_dict['obsHistID'][obs_mask[i_star]])

        star_data_dict['simobjid'].append(local_simobjid)
        star_data_dict['ra'].append(star_ra)
        star_data_dict['dec'].append(star_dec)
        for i_bp, bp in enumerate('ugrizy'):
            star_data_dict['quiescent_%s' % bp].append(quiescent_fluxes[:,i_bp])


def write_light_curve_metadata(out_dir, hpid, obs_dict, star_dict):
    """
    Write the metadata for a healpixel to disk

    Parameters
    ----------
    out_dir -- the directory where the output file will be written
    hpid -- the healpix pixel ID of
    obs_dict -- the dict of observation metadata corresponding to this hpid
    star_dict -- the dict of stellar metadata corresponding to this hpid

    Returns
    -------
    Nothing. Writes the metadata to out_dir/stellar_lc_{hpid}.h5

    Note 'stellar metadata' means static information about the stars
    in the healpixel
    """
    out_file_name = create_out_name(out_dir, hpid)
    print('writing: %s' % out_file_name)
    with h5py.File(out_file_name, 'a') as out_file:
        obsid = np.concatenate(obs_dict['obsHistID'])
        sorted_dex = np.argsort(obsid)
        for k in obs_dict:
            out_file.create_dataset(k,
                               data=np.concatenate(obs_dict[k][sorted_dex]))
        for k in star_dict:
            out_file.create_dataset(k, data=np.concatenate(star_dict[k]))


def find_lc_to_write(completed_hpid_list,
                     p_hpid_list,
                     obs_lock_hpid, obs_dict_hpid,
                     star_lock_hpid, star_dict_hpid):
    """
    Look for healpixels whose metadata can be written to disk because they are
    done being processed

    Parameters
    ----------
    completed_hpid_list -- a list of the healpixels whose entire field of view have
    been sent to do_photometry

    p_hpid_list -- the list of healpixels corresponding to the currently running
    do_photometry processes

    obs_lock_hpid -- the dict of locks for observation metadata (keyed on hpid)

    obs_dict_hpid -- the dict of output dicts for observation metadata (keyed on hpid)

    star_lock_hpid -- the dict of locks for stellar metadata (keyed on hpid)

    star_dict_hpid -- the dict of outputs dicts for stellar metdata (keyed on hpid)

    Returns
    --------
    a list of hpid that can be written to disk
    a list of corresponding observation metadata output dicts
    a list of corresponding stellar metadata output dicts
    """
    hpid_out = []
    star_dict_out = []
    obs_dict_out = []

    still_viable = list(star_lock_hpid.keys())

    # scan through the completed hpid;
    # if there are any who are not being run by a process
    # listed in p_hpid_list, add that healpix and its
    # output data to the list of outputs
    for hpid in completed_hpid_list:
        if not hpid in still_viable:
            continue

        is_done = True
        for hh in p_hpid_list:
            if hh == hpid:
                is_done = False
                break
        if not is_done:
            continue

        star_lock_hpid.pop(hpid)
        obs_lock_hpid.pop(hpid)
        star_dict = {}
        star_root = star_dict_hpid.pop(hpid)
        for k in list(star_root.keys()):
            star_dict[k] = list(star_root.pop(k))
        obs_dict = {}
        obs_root = obs_dict_hpid.pop(hpid)
        for k in list(obs_root.keys()):
            obs_dict[k] = list(obs_root.pop(k))


        hpid_out.append(hpid)
        star_dict_out.append(star_dict)
        obs_dict_out.append(obs_dict)

    return hpid_out, obs_dict_out, star_dict_out


if __name__ == "__main__":

    cache_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth'
    assert os.path.isdir(cache_dir)
    out_dir = os.path.join(cache_dir, 'dc2_stellar_lc')  # output directory

    # The file generated by obsHistID_to_healpix.py
    lookup_name = os.path.join(cache_dir, 'hpid_to_obsHistID_lookup.h5')
    assert os.path.isfile(lookup_name)

    # The database of stellar parameters
    stellar_db_name = os.path.join(cache_dir, 'dc2_stellar_healpixel.db')
    assert os.path.isfile(stellar_db_name)

    list_out_dir = os.listdir(out_dir)
    if len(list_out_dir)>0:
        msg = '\n\n%s\nis not empty\n\n'
        raise RuntimeError(msg)

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)


    # cache_LSST_seds will load (or create, if it doesn't exist)
    # a pickled file containing the entire LSST sims SED library.
    # It is important to call it with wavelen_min and wavelen_max set
    # because some of the SEDs extend far into the infrared, and
    # loading and storing them can be very resource intensive.
    # The pickle file will be saved in cache_dir for later use
    # (it takes about 15 minutes to create the file the first time)
    sims_photUtils.cache_LSST_seds(wavelen_min=0.0,
                                   wavelen_max=1500.0,
                                   cache_dir=cache_dir)

    # make sure that the Kepler light curves are loaded into memory
    kplr_dummy = variability.ParametrizedLightCurveMixin()
    kplr_dummy.load_parametrized_light_curves()

    # load the mapping between healpix pixel and obsHistID
    hpid_to_ct = {}
    with h5py.File(lookup_name, 'r') as in_file:
        for k in in_file.keys():
            try:
                hpid = int(k)
            except ValueError:
                continue

            n_obs = len(in_file['%d' % hpid][()])
            hpid_to_ct[hpid] = n_obs

    mgr = multiprocessing.Manager()

    # job_dict and job_lock exist to make sure that no more than
    # 5 processes are running the "generate delta magnitue" part
    # of the pipeline at once. This is the most memory intensive
    # part of the pipeline (requiring about 7 GB per process); it
    # is not actually the most time-consuming part, so it is alright
    # to limit the number of processes that can be running that step
    # at once.
    job_dict = mgr.dict()
    job_lock = mgr.Lock()

    print('created manager')
    time.sleep(30)
    print('moving on')

    # Each healpixel will have its own dict storing observation metadata
    # and stellar metadata; the locks prevent multiple processed from writing
    # to the same dict at once.
    obs_dict_hpid = {}
    obs_lock_hpid = {}
    star_dict_hpid = {}
    star_lock_hpid = {}

    job_dict['running_dmag'] = 0

    n_threads = 20

    p_list = []  # list of do_photometry processes
    p_hpid_list = []
    w_list = []  # list of processes writing metadata to disk
    completed_hpid_list = []
    with sqlite3.connect(stellar_db_name) as conn:
        cursor = conn.cursor()

        ct =0
        t_start = time.time()

        # sort healpix pixels so that you start with the pixels
        # that are visited most often
        hpid_to_ct_items = hpid_to_ct.items()
        hpid_key_list = np.array([itm[0] for itm in hpid_to_ct_items])
        hpid_val_list = np.array([itm[1] for itm in hpid_to_ct_items])
        sorted_dex = np.argsort(-1*hpid_val_list)
        hpid_key_list = hpid_key_list[sorted_dex]

        for hpid in hpid_key_list:

            # create the output dicts and locks for this healpixel
            obs_lock = mgr.Lock()
            star_lock = mgr.Lock()
            obs_metadata_dict = mgr.dict()
            star_data_dict = mgr.dict()

            obs_metadata_dict['mjd'] = mgr.list()
            obs_metadata_dict['obsHistID'] = mgr.list()
            obs_metadata_dict['filter'] = mgr.list()

            star_data_dict['simobjid'] = mgr.list()
            star_data_dict['ra'] = mgr.list()
            star_data_dict['dec'] = mgr.list()
            for bp in 'ugrizy':
                star_data_dict['quiescent_%s' % bp] = mgr.list()

            obs_lock_hpid[hpid] = obs_lock
            obs_dict_hpid[hpid] = obs_metadata_dict
            star_lock_hpid[hpid] = star_lock
            star_dict_hpid[hpid] = star_data_dict

            # limit chunksize so there aren't more than
            # 50 million (star,mjd) pairs being processed at once
            chunk_size = 50000000//hpid_to_ct[hpid]
            print('chunk_size %d' % chunk_size)

            query = "SELECT simobjid, hpid, sedFilename, magNorm, ebv, "
            query += "varParamStr, parallax, ra, decl FROM stars "
            query += "WHERE hpid=%d" % hpid

            data_iterator = cursor.execute(query)
            chunk = data_iterator.fetchmany(chunk_size)
            while len(chunk)>0:
                p = multiprocessing.Process(target=do_photometry,
                                      args=(chunk,
                                            obs_lock_hpid[hpid],
                                            obs_dict_hpid[hpid],
                                            star_lock_hpid[hpid],
                                            star_dict_hpid[hpid],
                                            job_lock, job_dict,
                                            out_dir))
                p.start()
                p_list.append(p)
                p_hpid_list.append(hpid)

                # if all of the desired threads are taken, wait for a job
                # to exit
                while len(p_list) + len(w_list)>=n_threads:
                    e_list = list([p.exitcode for p in p_list])
                    was_popped = False
                    for ii in range(len(e_list)-1,-1,-1):
                        if e_list[ii] is not None:
                            assert p_list[ii].exitcode is not None
                            p_list.pop(ii)
                            p_hpid_list.pop(ii)
                            was_popped = True
                            ct += chunk_size

                    e_list = list([w.exitcode for w in w_list])
                    for ii in range(len(e_list)-1,-1,-1):
                        if e_list[ii] is not None:
                            assert w_list[ii].exitcode is not None
                            w_list.pop(ii)

                    if was_popped:
                        duration = (time.time()-t_start)/3600.0
                        per = duration/ct
                        print('ran %d in %e hrs (%e)' % (ct, duration, per))

                chunk = data_iterator.fetchmany(chunk_size)

            completed_hpid_list.append(hpid)

            # look to see if any healpixels are ready to be written
            # to disk
            (hpid_to_write,
             obs_to_write,
             star_to_write) = find_lc_to_write(completed_hpid_list,
                                               p_hpid_list,
                                               obs_lock_hpid,
                                               obs_dict_hpid,
                                               star_lock_hpid,
                                               star_dict_hpid)

            for hh, oo, ss in zip(hpid_to_write,
                                  obs_to_write,
                                  star_to_write):

                w = multiprocessing.Process(target=write_light_curve_metadata,
                                            args=(out_dir, hh, oo, ss))
                w.start()
                w_list.append(w)


    for p in p_list:
        p.join()
    p_hpid_list = []

    (hpid_to_write,
     obs_to_write,
     star_to_write) = find_lc_to_write(completed_hpid_list,
                                       p_hpid_list,
                                       obs_lock_hpid,
                                       obs_dict_hpid,
                                       star_lock_hpid,
                                       star_dict_hpid)

    for hh, oo, ss in zip(hpid_to_write,
                          obs_to_write,
                          star_to_write):

        w = multiprocessing.Process(target=write_light_curve_metadata,
                                    args=(out_dir, hh, oo, ss))
        w.start()
        w_list.append(w)

    for w in w_list:
        w.join()

    print('that took %e hrs ' % ((time.time()-t_start)/3600.0))
