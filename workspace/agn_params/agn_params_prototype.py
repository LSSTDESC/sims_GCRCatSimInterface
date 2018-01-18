import numpy as np
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

    restrame_wavelen_grid = bp.wavelen/dilation

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
