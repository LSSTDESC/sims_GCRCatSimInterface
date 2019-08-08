import os
import numpy as np
import h5py
import time

import multiprocessing

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.utils import _getRotSkyPos
from lsst.sims.utils import angularSeparation

from lsst.sims.coordUtils import chipNameFromRaDecLSST

def find_lightcurve_times(gid, ra, dec, obs_list, out_dir, my_lock_dict,
                          d_dex, log_file_name):

    #print('starting %d %d %d -- %d' %
    #(os.getpid(),gid.min(),gid.max(),len(gid)))

    bp_to_int = {}
    for i_bp, bp in enumerate('ugrizy'):
        bp_to_int[bp] = i_bp
    mjd_arr = np.array([obs.mjd.TAI for obs in obs_list])
    bp_arr = np.array([bp_to_int[obs.bandpass] for obs in obs_list])

    local_dict = {}
    for ii in gid:
        local_dict[ii] = -1*np.ones(len(obs_list))

    t_start = time.time()
    ct = 0
    for i_obs, obs in enumerate(obs_list):
        dist = angularSeparation(ra, dec,
                                 obs.pointingRA, obs.pointingDec)

        remotely_possible = np.where(dist<2.2)
        remotely_possible_gid = gid[remotely_possible]

        name_arr = chipNameFromRaDecLSST(ra[remotely_possible],
                                         dec[remotely_possible],
                                         obs_metadata=obs,
                                         band=obs.bandpass)

        valid = np.where(np.char.find(name_arr.astype(str), 'None')<0)
        for ii in gid[valid]:
            local_dict[ii][i_obs] = 1

    min_dex = gid.min()//d_dex
    max_dex = gid.max()//d_dex

    for name_dex in range(min_dex, max_dex+1):
        min_gid = name_dex*d_dex
        max_gid = (name_dex+1)*d_dex
        small_name = 'agn_lc_%d_%d.h5' % (min_gid,max_gid-1)

        out_name = os.path.join(out_dir, small_name)

        if not os.path.isfile(out_name):
            open_mode = 'w'
        else:
            open_mode = 'a'

        valid_gid = np.where(np.logical_and(gid>=min_gid, gid<max_gid))

        if len(valid_gid[0])==0:
            continue

        with my_lock_dict[small_name] as locked_context:
            with h5py.File(out_name, open_mode) as out_file:

                #print('writing %d -- %s' % (os.getpid(), small_name))
                if open_mode == 'w':
                    out_file.create_dataset('mjd', data=mjd_arr)
                    out_file.create_dataset('bp', data=bp_arr)

                for ii in gid[valid_gid]:
                    valid = np.where(local_dict[ii]>0)[0]
                    if len(valid)>0:
                        try:
                            out_file.create_dataset('%d' % ii, data=valid)
                        except RuntimeError:
                            msg = 'out_dir: %s\n' % out_dir
                            msg += 'small_name: %s\n' % small_name
                            msg += 'already has: %d\n\n' % ii
                            raise RuntimeError(msg)

    with my_lock_dict[log_file_name] as locked_context:
        with open(log_file_name, 'a') as out_file:
            for ii in gid:
                out_file.write('%d\n' % ii)


if __name__ == "__main__":
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'

    out_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth/agn_190802_test'
    assert os.path.isdir(out_dir)

    log_file_name = os.path.join(out_dir, 'agn_gid_lc_log.txt')

    d_dex = 500000000

    n_proc = 20

    data_dir = '/astro/store/pogo3/danielsf/desc_dc2_truth'
    assert os.path.isdir(data_dir)

    agn_data_dir = os.path.join(data_dir, 'agn_190730')
    assert os.path.isdir(agn_data_dir)

    wfd_file_name = '/local/lsst/danielsf/sims_GCRCatSimInterface'
    wfd_file_name = os.path.join(wfd_file_name, 'workspace/run2.1/data')
    wfd_file_name = os.path.join(wfd_file_name, 'wfd_obshistid_list.txt')
    assert os.path.isfile(wfd_file_name)

    minion_name = os.path.join(data_dir, 'minion_1016_desc_dithered_v4_sfd.db')
    assert os.path.isfile(minion_name)

    agn_support_data_name = os.path.join(agn_data_dir, 'agn_gid_z_ra_dec.h5')
    assert os.path.isfile(agn_support_data_name)

    gid_per_proc = 10000

    with h5py.File(agn_support_data_name, 'r') as agn_data:
        agn_ra = agn_data['ra'][()]
        agn_dec = agn_data['dec'][()]
        agn_gid = agn_data['galaxy_id'][()]

    sorted_dex = np.argsort(agn_gid)
    agn_gid = agn_gid[sorted_dex]
    agn_ra = agn_ra[sorted_dex]
    agn_dec = agn_dec[sorted_dex]

    shuffle_dexes = []
    chunk_dexes = []
    ct = 0
    for i_start in range(0,agn_gid.max(),d_dex):
        valid = np.where(np.logical_and(agn_gid>=i_start,
                                        agn_gid<i_start+d_dex))
        chunk_dexes.append(valid[0])
        ct += len(chunk_dexes[-1])

    final_n_visits = ct

    shuffle_ct = 0
    done_shuffling = False
    while not done_shuffling:
        done_shuffling = True
        for i_list, chunk_list in enumerate(chunk_dexes):
            if len(chunk_list) == 0:
                continue

            shuffle_ct += len(chunk_list[:gid_per_proc])
            shuffle_dexes.append(np.copy(chunk_list[:gid_per_proc]))
            chunk_list = chunk_list[gid_per_proc:]
            chunk_dexes[i_list] = chunk_list
            if len(chunk_list)>0:
                done_shuffling = False
            if(shuffle_ct>final_n_visits+1000000):
                raise RuntimeError("too many shuffle dexes")

    shuffle_dexes = np.concatenate(shuffle_dexes)
    valid = np.where(shuffle_dexes<len(agn_gid))
    shuffle_dexes = shuffle_dexes[valid]
    np.testing.assert_array_equal(np.sort(shuffle_dexes),
                                  np.arange(len(agn_gid)))
    agn_gid = agn_gid[shuffle_dexes]
    agn_ra = agn_ra[shuffle_dexes]
    agn_dec = agn_dec[shuffle_dexes]

    assert len(np.unique(agn_gid)) == len(agn_gid)

    obshistid_list = []
    with open(wfd_file_name, 'r') as in_file:
        for line in in_file:
            obshistid_list.append(int(line))

    n_obs = len(obshistid_list)
    ct_obs = 0

    obs_gen = ObservationMetaDataGenerator(minion_name)
    t_start = time.time()
    obs_list = []
    print('total %e' % len(obshistid_list))

    for ii in obshistid_list:
        ct_obs += 1
        obs = obs_gen.getObservationMetaData(obsHistID=ii)[0]

        obs.pointingRA = np.degrees(obs.OpsimMetaData['descDitheredRA'])
        obs.pointingDec = np.degrees(obs.OpsimMetaData['descDitheredDec'])
        obs.OpsimMetaData['rotTelPos'] = obs.OpsimMetaData['descDitheredRotTelPos']
        rotsky = np.degrees(_getRotSkyPos(obs._pointingRA, obs._pointingDec,
                                          obs,
                                          obs.OpsimMetaData['rotTelPos']))
        obs.rotSkyPos = rotsky
        obs_list.append(obs)

        if len(obs_list)%1000 == 0:
            print(len(obs_list)) 


    p_list = []
    mgr = multiprocessing.Manager()
    my_lock_dict = mgr.dict()

    for i_start in range(0, agn_gid.max(), d_dex):
        out_name = 'agn_lc_%d_%d.h5' % (i_start,i_start+d_dex-1)
        my_lock_dict[out_name] = mgr.Lock()

    my_lock_dict[log_file_name] = mgr.Lock()

    # for testing

    #test_dexes = np.concatenate([np.arange(0,10000),
    #                             np.arange(110000,120000),
    #                             np.arange(210000,220000),
    #                             np.arange(310000,310000),
    #                             np.arange(420000,430000),
    #                             np.arange(550000,560000),
    #                             np.arange(660000,670000),
    #                             np.arange(680000,4000000)])
    #agn_gid = agn_gid[test_dexes]
    #agn_ra = agn_ra[test_dexes]
    #agn_dec = agn_dec[test_dexes]

    #q = slice(0,-1,31)
    #agn_gid = agn_gid[q]
    #agn_ra = agn_ra[q]
    #agn_dec = agn_dec[q]
    #rng = np.random.RandomState(81623)
    #chosen_dex = rng.choice(np.arange(0,len(agn_gid),dtype=int),
    #                        size=1000000, replace=False)
    #agn_gid = agn_gid[chosen_dex]
    #agn_ra = agn_ra[chosen_dex]
    #agn_dec = agn_dec[chosen_dex]
    #sorted_dex = np.argsort(agn_gid)
    #agn_gid = agn_gid[sorted_dex]
    #agn_ra = agn_ra[sorted_dex]
    #agn_dec = agn_dec[sorted_dex]

    for i_start in range(0, len(agn_gid), gid_per_proc):
        local_gid = agn_gid[i_start:i_start+gid_per_proc]
        local_ra = agn_ra[i_start:i_start+gid_per_proc]
        local_dec = agn_dec[i_start:i_start+gid_per_proc]

        p = multiprocessing.Process(target=find_lightcurve_times,
                                    args=[local_gid, local_ra, local_dec,
                                          obs_list, out_dir, my_lock_dict,
                                          d_dex, log_file_name])
        p.start()
        p_list.append(p)

        while len(p_list)>=n_proc:
            exit_code_list = []
            for p in p_list:
                exit_code_list.append(p.exitcode)
            for i_p in range(len(exit_code_list)-1, -1, -1):
                if exit_code_list[i_p] is not None:
                    p_list.pop(i_p)

    for p in p_list:
        p.join()
    print('all done')
    time.sleep(30)
