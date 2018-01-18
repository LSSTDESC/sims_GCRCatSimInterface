import numpy as np

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


if __name__ == "__main__":

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
