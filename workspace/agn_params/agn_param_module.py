import numpy as np
import os
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, BandpassDict

def log_Eddington_ratio(bhmass, accretion_rate):
    """
    Parameters
    ----------
    bhmass is in solar masses

    accretion_rate is in solar masses per Gyr

    Returns
    -------
    log10(L/L_Eddington)
    """
    # reference for expressions defining Eddington luminosity
    # http://www-astro.physics.ox.ac.uk/~garret/teaching/lecture7-2012.pdf

    log_m_sun = np.log10(1.98855) + 33.0  # in grams
    log_G = np.log10(6.674) - 8.0  # in cgs units
    log_m_proton = log_m_sun - np.log10(1.189) - 57.0  # see A.2 of Kolb and Turner
    log_sigma_T = np.log10(6.6524) - 25.0  # in cm^2 -- A.1.2 of Kolb and Turner
    log_c = np.log10(2.9979) + 10.0  # in cm/sec
    log_sec_per_yr = np.log10(3600.0*24.0*365.25)

    log_Const = np.log10(4.0) + np.log10(np.pi)
    log_Const += log_G + log_c + log_m_proton - log_sigma_T
    log_L_Eddington = np.log10(bhmass) + log_m_sun + log_Const

    log_epsilon = -1.0
    log_L = np.log10(accretion_rate)
    log_L += log_m_sun - 9.0 - log_sec_per_yr
    log_L += log_epsilon
    log_L += 2.0*log_c

    output = log_L-log_L_Eddington
    return output

def M_i_from_L_Mass(Ledd_ratio, bhmass):
    """
    Parameters
    ----------
    Ledd_ratio is the log10(L/L_Eddington) ratio

    bhmass is the log10(mass of the blackhole in solar masses)

    Returns
    -------
    Absolute i-band magnitude.  This will be read off from
    the apparent relationships in Figure 15 of MacLeod et al 2010
    (ApJ, 721, 1014)
    """

    if not hasattr(M_i_from_L_Mass, '_initialized'):
        M_i_from_L_Mass._initialized = True

        # example points taken from Figure 15 of MacLeod et al (2010)
        l_edd = [-0.5, -0.5,
                 -0.1, -0.1,
                 -1.1, -1.1,
                 -1.5, -1.5]
        mbh = [9.8, 7.8,
               9.0, 7.7,
               10.1, 8.3,
               10.1, 8.85]
        m_i = [-28.3, -23.2,
               -27.6, -24.4,
               -27.7, -23.2,
               -26.3, -23.1]

        l_edd = np.array(l_edd)
        mbh = np.array(mbh)
        m_i = np.array(m_i)

        theta_best = None
        l_edd_0_best = None
        mbh_0_best = None
        err_best = None
        mm = np.zeros((2,2), dtype=float)
        bb = np.zeros(2, dtype=float)
        nn = len(m_i)
        for theta in np.arange(0.0, 2.0*np.pi, 0.01):
            mm[0][0] = -1.0*nn*np.cos(theta)**2
            mm[0][1] = nn*np.sin(theta)*np.cos(theta)
            mm[1][0] = -1.0*nn*np.sin(theta)*np.cos(theta)
            mm[1][1] = nn*np.sin(theta)**2

            bb_sum = (m_i-l_edd*np.cos(theta)+mbh*np.sin(theta)).sum()
            bb[0] = np.cos(theta)*bb_sum
            bb[1] = np.sin(theta)*bb_sum

            try:
                vv = np.linalg.solve(mm, bb)
            except np.linalg.LinAlgError:
                continue
            err = ((m_i-np.cos(theta)*(l_edd-vv[0])+np.sin(theta)*(mbh-vv[1]))**2).sum()
            #print('err %e' % (np.sqrt(err/nn)))
            if err_best is None or err<err_best:
                err_best = err
                theta_best=theta
                l_edd_0_best = vv[0]
                mbh_0_best = vv[1]

        M_i_from_L_Mass._c_theta = np.cos(theta_best)
        M_i_from_L_Mass._s_theta = np.sin(theta_best)
        M_i_from_L_Mass._ledd_0 = l_edd_0_best
        M_i_from_L_Mass._mbh_0 = mbh_0_best
        print('err_best %e' % (np.sqrt(err_best/len(l_edd))))
        print('theta best %e slope %e\n' % (theta_best, np.tan(theta_best)))

    return (M_i_from_L_Mass._c_theta*(Ledd_ratio-M_i_from_L_Mass._ledd_0) -
            M_i_from_L_Mass._s_theta*(bhmass-M_i_from_L_Mass._mbh_0))


def k_correction(sed_obj, bp, redshift):
    """
    Parameters
    ----------
    sed_obj is an instantiation of Sed representing the observed
    spectral energy density of the source

    bp is an instantiation of Bandpass representing the bandpass
    in which we are calculating the magnitudes

    redshift is a float representing the redshift of the source

    Returns
    -------
    K correction in magnitudes according to equation (12) of
    Hogg et al. 2002 (arXiv:astro-ph/0210394)
    """
    if bp.phi is None:
        bp.sbToPhi()
    if sed_obj.fnu is None:
        sed_obj.flambdaTofnu()

    dilation = 1.0 + redshift

    restframe_wavelen_grid = bp.wavelen*dilation

    valid_bp_dex = np.where(np.abs(bp.sb)>0.0)
    valid_restframe_wavelen = restframe_wavelen_grid[valid_bp_dex]
    restframe_min_wavelen = valid_restframe_wavelen.min()
    restframe_max_wavelen = valid_restframe_wavelen.max()
    try:
        assert restframe_min_wavelen > sed_obj.wavelen.min()
        assert restframe_max_wavelen < sed_obj.wavelen.max()
    except:
        msg = '\nBP/(1+z) range '
        msg += '%.6e < lambda < %.6e\n' % (restframe_min_wavelen,
                                           restframe_max_wavelen)
        msg += 'SED range '
        mst += '%.6e < lambda < %.6e\n' % (sed_obj.wavelen.min(),
                                           sed_obj.wavelen.max())

        raise RuntimeError(msg)


    restframe_fnu = np.interp(restframe_wavelen_grid,
                              sed_obj.wavelen,
                              sed_obj.fnu,
                              left=0.0,
                              right=0.0)

    observed_fnu = np.interp(bp.wavelen,
                             sed_obj.wavelen,
                             sed_obj.fnu,
                             left=0.0,
                             right=0.0)

    restframe_integral = (0.5*(bp.sb[1:]*restframe_fnu[1:]/bp.wavelen[1:] +
                               bp.sb[:-1]*restframe_fnu[:-1]/bp.wavelen[:-1]) *
                              (bp.wavelen[1:]-bp.wavelen[:-1])).sum()

    observer_integral = (0.5*(bp.sb[1:]*observed_fnu[1:]/bp.wavelen[1:] +
                              bp.sb[:-1]*observed_fnu[:-1]/bp.wavelen[:-1]) *
                             (bp.wavelen[1:]-bp.wavelen[:-1])).sum()

    return -2.5*np.log10((1.0+redshift)*observer_integral/restframe_integral)


if __name__ == "__main__":

    # below is the code I used to test the log_Eddington_ratio method
    #################################################################
    """
    dtype = np.dtype([('bhmass', float), ('accretion_rate', float)])
    data = np.genfromtxt('data/proto_dc2_bh_params.txt', dtype=dtype)
    valid = np.where(np.logical_and(data['bhmass']!=0.0,
                                    data['accretion_rate']!=0.0))

    data = data[valid]
    log_rat = log_Eddington_ratio(data['bhmass'], data['accretion_rate'])

    #valid = np.where(np.logical_not(np.isinf(log_rat)))
    #log_rat = log_rat[valid]
    print(len(log_rat),
          len(np.where(np.isnan(log_rat))[0]),
          len(np.where(np.isinf(log_rat))[0]))

    sorted_dex = np.argsort(log_rat)
    log_rat = log_rat[sorted_dex]
    data = data[sorted_dex]
    n_rat = len(log_rat)
    print(log_rat[0])
    print(log_rat[n_rat//4])
    print(log_rat[n_rat//2])
    print(log_rat[(3*n_rat)//4])
    print(log_rat[-1])
    print(n_rat, len(np.where(log_rat>-2.0)[0]))
    """

    # below is the code used to test the K correction
    #################################################
    """
    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    rng = np.random.RandomState(41321)
    sed_dir = os.path.join(getPackageDir('sims_sed_library'), 'galaxySED')
    list_of_sed_files = os.listdir(sed_dir)
    list_of_sed_files.sort()
    sed_to_check = rng.choice(list_of_sed_files, size=10)
    redshift_arr = rng.random_sample(len(sed_to_check))*2.0+0.1

    bp = bp_dict['g']
    for sed_name, zz in zip(sed_to_check, redshift_arr):
        full_name = os.path.join(sed_dir, sed_name)
        ss = Sed()
        ss.readSED_flambda(full_name)
        true_rest_mag = ss.calcMag(bp)
        ss.redshiftSED(zz, dimming=True)
        obs_mag = ss.calcMag(bp)
        k_corr = k_correction(ss, bp, zz)
        print(true_rest_mag, obs_mag, k_corr, obs_mag-k_corr, zz)
    """

    print(M_i_from_L_Mass(-0.5, 8.7))
