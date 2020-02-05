from builtins import range
from builtins import object
import numpy as np
import linecache
import math
import os
import gzip
import numbers
import multiprocessing
import json as json
import numpy as np
from lsst.sims.catalogs.decorators import register_method, compound
from lsst.sims.catUtils.mixins import Variability

__all__ = ["ExtraGalacticVariabilityModels", "VariabilityAGN"]

class ExtraGalacticVariabilityModels(Variability):
    """
    A mixin providing the model for AGN variability.
    """

    _agn_walk_start_date = 58580.0
    _agn_threads = 1

    @register_method('applyAgn')
    def applyAgn(self, valid_dexes, params, expmjd,
                 variability_cache=None, redshift=None):

        if redshift is None:
            redshift_arr = self.column_by_name('redshift')
        else:
            redshift_arr = redshift

        if len(params) == 0:
            return np.array([[],[],[],[],[],[]])

        if isinstance(expmjd, numbers.Number):
            dMags = np.zeros((6, self.num_variable_obj(params)))
            max_mjd = expmjd
            min_mjd = expmjd
            mjd_is_number = True
        else:
            dMags = np.zeros((6, self.num_variable_obj(params), len(expmjd)))
            max_mjd = max(expmjd)
            min_mjd = min(expmjd)
            mjd_is_number = False

        seed_arr = params['seed']
        # tau_u_arr = params['agn_tau_u'].astype(float)
        # tau_g_arr = params['agn_tau_g'].astype(float)
        # tau_r_arr = params['agn_tau_r'].astype(float)
        # tau_i_arr = params['agn_tau_i'].astype(float)
        # tau_z_arr = params['agn_tau_z'].astype(float)
        # tau_y_arr = params['agn_tau_y'].astype(float)
        # sfu_arr = params['agn_sfu'].astype(float)
        # sfg_arr = params['agn_sfg'].astype(float)
        # sfr_arr = params['agn_sfr'].astype(float)
        # sfi_arr = params['agn_sfi'].astype(float)
        # sfz_arr = params['agn_sfz'].astype(float)
        # sfy_arr = params['agn_sfy'].astype(float)

        duration_observer_frame = max_mjd - self._agn_walk_start_date

        if duration_observer_frame < 0 or min_mjd < self._agn_walk_start_date:
            raise RuntimeError("WARNING: Time offset greater than minimum epoch.  " +
                               "Not applying variability. "+
                               "expmjd: %e should be > start_date: %e  " % (min_mjd, self._agn_walk_start_date) +
                               "in applyAgn variability method")

        if self._agn_threads == 1 or len(valid_dexes[0])==1:
            for filt_num, filt_name in list(enumerate(['u', 'g', 'r', 'i', 'z', 'y'])):
                tau_arr = params['agn_tau_%s' % filt_name].astype(float)
                sf_arr = params['agn_sf_%s' % filt_name].astype(float)
                for i_obj in valid_dexes[0]:
                    seed = seed_arr[i_obj]
                    tau_filt = tau_arr[i_obj]
                    time_dilation = 1.0+redshift_arr[i_obj]
                    sf_filt = sf_arr[i_obj]
                    dMags[filt_num][i_obj] = self._simulate_agn(expmjd, tau_filt, time_dilation, sf_filt, seed)
        else:
            p_list = []

            mgr = multiprocessing.Manager()
            if mjd_is_number:
                out_struct = mgr.Array('d', [0]*len(valid_dexes[0]))
            else:
                out_struct = mgr.dict()

            #################
            # Try to subdivide the AGN into batches such that the number
            # of time steps simulated by each thread is close to equal
            tot_steps = 0
            n_steps = []
            for tt, zz in zip(tau_arr[valid_dexes], redshift_arr[valid_dexes]):
                dilation = 1.0+zz
                dt = tt/100.0
                dur = (duration_observer_frame/dilation)
                nt = dur/dt
                tot_steps += nt
                n_steps.append(nt)

            batch_target = tot_steps/self._agn_threads
            i_start_arr = [0]
            i_end_arr = []
            current_batch = n_steps[0]
            for ii in range(1,len(n_steps),1):
                current_batch += n_steps[ii]
                if ii == len(n_steps)-1:
                    i_end_arr.append(len(n_steps))
                elif len(i_start_arr)<self._agn_threads:
                    if current_batch>=batch_target:
                        i_end_arr.append(ii)
                        i_start_arr.append(ii)
                        current_batch = n_steps[ii]

            if len(i_start_arr) != len(i_end_arr):
                raise RuntimeError('len i_start %d len i_end %d; dexes %d' %
                                   (len(i_start_arr),
                                    len(i_end_arr),
                                    len(valid_dexes[0])))
            assert len(i_start_arr) <= self._agn_threads
            ############

            # Actually simulate the AGN on the the number of threads allotted
            for i_start, i_end in zip(i_start_arr, i_end_arr):
                dexes = valid_dexes[0][i_start:i_end]
                if mjd_is_number:
                    out_dexes = range(i_start,i_end,1)
                else:
                    out_dexes = dexes
                p = multiprocessing.Process(target=self._threaded_simulate_agn,
                                            args=(expmjd, tau_arr[dexes],
                                                  1.0+redshift_arr[dexes],
                                                  sfu_arr[dexes],
                                                  seed_arr[dexes],
                                                  out_dexes,
                                                  out_struct))
                p.start()
                p_list.append(p)
            for p in p_list:
                p.join()

            if mjd_is_number:
                dMags[0][valid_dexes] = out_struct[:]
            else:
                for i_obj in out_struct.keys():
                    dMags[0][i_obj] = out_struct[i_obj]

        # for i_filter, filter_name in enumerate(('g', 'r', 'i', 'z', 'y')):
        #     for i_obj in valid_dexes[0]:
        #         dMags[i_filter+1][i_obj] = dMags[0][i_obj]*params['agn_sf%s' % filter_name][i_obj]/params['agn_sfu'][i_obj]

        return dMags

    def _threaded_simulate_agn(self, expmjd, tau_arr,
                               time_dilation_arr, sf_u_arr,
                               seed_arr, dex_arr, out_struct):

        if isinstance(expmjd, numbers.Number):
            mjd_is_number = True
        else:
            mjd_is_number = False

        for tau, time_dilation, sf_u, seed, dex in \
        zip(tau_arr, time_dilation_arr, sf_u_arr, seed_arr, dex_arr):
            out_struct[dex] = self._simulate_agn(expmjd, tau, time_dilation,
                                                 sf_u, seed)

    def _simulate_agn(self, expmjd, tau, time_dilation, sf_u, seed):
            """
            Simulate the u-band light curve for a single AGN

            Parameters
            ----------
            expmjd -- a number or numpy array of dates for the light curver

            tau -- the characteristic timescale of the AGN in days

            time_dilation -- (1+z) for the AGN

            sf_u -- the u-band structure function of the AGN

            seed -- the seed for the random number generator

            Returns
            -------
            a numpy array (or number) of delta_magnitude in the u-band at expmjd
            """

            if not isinstance(expmjd, numbers.Number):
                d_m_out = np.zeros(len(expmjd))
                duration_observer_frame = max(expmjd) - self._agn_walk_start_date
            else:
                duration_observer_frame = expmjd - self._agn_walk_start_date


            rng = np.random.RandomState(seed)
            dt = tau/100.
            duration_rest_frame = duration_observer_frame/time_dilation
            nbins = int(math.ceil(duration_rest_frame/dt))+1

            time_dexes = np.round((expmjd-self._agn_walk_start_date)/(time_dilation*dt)).astype(int)
            time_dex_map = {}
            ct_dex = 0
            if not isinstance(time_dexes, numbers.Number):
                for i_t_dex, t_dex in enumerate(time_dexes):
                    if t_dex in time_dex_map:
                        time_dex_map[t_dex].append(i_t_dex)
                    else:
                        time_dex_map[t_dex] = [i_t_dex]
                time_dexes = set(time_dexes)
            else:
                time_dex_map[time_dexes] = [0]
                time_dexes = set([time_dexes])

            dx2 = 0.0
            x1 = 0.0
            x2 = 0.0

            dt_over_tau = dt/tau
            es = rng.normal(0., 1., nbins)*math.sqrt(dt_over_tau)
            for i_time in range(nbins):
                #The second term differs from Zeljko's equation by sqrt(2.)
                #because he assumes stdev = sf_u/sqrt(2)
                dx1 = dx2
                dx2 = -dx1*dt_over_tau + sf_u*es[i_time] + dx1
                x1 = x2
                x2 += dt

                if i_time in time_dexes:
                    if isinstance(expmjd, numbers.Number):
                        dm_val = ((expmjd-self._agn_walk_start_date)*(dx1-dx2)/time_dilation+dx2*x1-dx1*x2)/(x1-x2)
                        d_m_out = dm_val
                    else:
                        for i_time_out in time_dex_map[i_time]:
                            local_end = (expmjd[i_time_out]-self._agn_walk_start_date)/time_dilation
                            dm_val = (local_end*(dx1-dx2)+dx2*x1-dx1*x2)/(x1-x2)
                            d_m_out[i_time_out] = dm_val

            return d_m_out


class VariabilityAGN(ExtraGalacticVariabilityModels):
    """
    This is a mixin which wraps the methods from the class
    ExtraGalacticVariabilityModels into getters for InstanceCatalogs
    of AGN.  Getters in this method should define columns named like

    delta_columnName

    where columnName is the name of the baseline (non-varying) magnitude
    column to which delta_columnName will be added.  The getters in the
    photometry mixins will know to find these columns and add them to
    columnName, provided that the columns here follow this naming convention.

    Thus: merely including VariabilityStars in the inheritance tree of
    an InstanceCatalog daughter class will activate variability for any column
    for which delta_columnName is defined.
    """

    @compound('delta_lsst_u', 'delta_lsst_g', 'delta_lsst_r',
             'delta_lsst_i', 'delta_lsst_z', 'delta_lsst_y')
    def get_stellar_variability(self):
        """
        Getter for the change in magnitudes due to stellar
        variability.  The PhotometryStars mixin is clever enough
        to automatically add this to the baseline magnitude.
        """

        varParams = self.column_by_name('varParamStr')
        dmag = self.applyVariability(varParams)
        if dmag.shape != (6, len(varParams)):
            raise RuntimeError("applyVariability is returning "
                               "an array of shape %s\n" % dmag.shape
                               + "should be (6, %d)" % len(varParams))
        return dmag
