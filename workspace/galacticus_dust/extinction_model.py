""" Module implementing the extinction_curve function
approximating the treatment of dust in Galacticus.
"""
import numpy as np
from scipy.stats import norm

global_logx0 = 3.75
global_Acut = 0.98


def extinction_curve(wavelength, A, logx0=global_logx0,
            C=0.9, alpha=0.0002, Acut=global_Acut, bump_height=None):
    """ Fraction of the restframe flux absorbed by dust as a function of wavelength.

    extincted_flux = extinction_curve * non_extincted_flux

    Parameters
    ----------
    wavelength : ndarray of shape (nbins, )

    A : float
        Parameter encoding strength of extinction.
        Allowed values are in the interval (0, 1], with larger values
        corresponding to less extinction.

    Returns
    -------
    extinction : ndarray of shape (nbins, )
        All entries will be in the range [0, 1].
        Smaller values correspond to stronger extinction.
    """
    x = np.atleast_1d(wavelength)

    if A >= Acut:
        return np.ones_like(x)
    result = np.zeros_like(x)

    beta = calculate_beta_from_A(A)
    xp = calculate_powerlaw_intercept(A, beta, C, logx0)
    asymptotic_mask = x >= xp

    D = 1. - C
    if alpha is None:
        alpha = powerlaw_derivative(xp, beta, A)/D

    result[asymptotic_mask] = asymptotic_behavior(x[asymptotic_mask], xp, alpha, C)

    result[~asymptotic_mask] = powerlaw(x[~asymptotic_mask], beta, A)

    if bump_height is None:
        bump_height = calculate_bump_height(A)
    bump = dust_grain_bump(x, bump_height)

    return result*bump


def powerlaw(x, beta, A, logx0=global_logx0, Acut=global_Acut):
    logx = np.log10(x)
    log_result = np.log(A) - beta*(logx-logx0)
    if A >= Acut:
        return np.ones_like(x)
    else:
        return np.exp(log_result)


def calculate_powerlaw_params(w, data, logx0=global_logx0):
    A = calculate_A(w, data, logx0)
    beta = calculate_beta_from_A(A)
    return A, beta


def calculate_A(w, data, logx0=global_logx0):
    logw = np.log10(w)
    mask = (logw > logx0-0.2) & (logw < logx0+0.2)
    c1, A = np.polyfit(logw[mask]-logx0, data[mask], deg=1)
    return A


def calculate_beta_from_A(A, c0=-2.333, c1=2.155):
    return c0 + c1*A


def powerlaw_derivative(x, beta, A, logx0=global_logx0):
    x0 = 10**logx0
    return (-beta*A/x0)*((x/x0)**(-beta-1.))


def calculate_powerlaw_intercept(A, beta, C, logx0):
    return 10**((np.log(A/C)/beta) + logx0)


def asymptotic_behavior(x, xp, alpha, C):
    D = 1. - C
    return 1. - D*(np.exp(-alpha*(x-xp)))


def calculate_bump_height(A):
    bump_height = 100.
    return bump_height


def dust_grain_bump(w, bump_height=100., scale=0.02, loc=3.125):
    logw = np.log10(w)
    s = np.linspace(logw.min(), logw.max(), len(w))
    y = norm.pdf(s, loc=loc, scale=scale)
    z = 1.-y/float(bump_height)
    return z
