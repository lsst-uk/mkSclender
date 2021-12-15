import numpy as np
import math

def ranrm(theta):
    '''
    Normalise an angle expressed in radians into the range 0 to 2pi.
    '''
    w = math.fmod(theta, 2.0 * math.pi)
    if w < 0.0: w = w + 2.0 * math.pi
    return w

def tangent_plane_to_spherical(xi, eta, raz, dcz):
    '''
    Takes the local plane coordinates (xi, eta) and de-projects back onto the celestial sphere using the given
    tangent point origin (raz, dcz). Algorithm cribbed from Pat Wallaces' dtp2s.c in SLALIB.

    Returns spherical coordinates RA, Dec.

    All arguments and results in RADIANS.
    '''
    sdecz = math.sin (dcz)
    cdecz = math.cos (dcz)
    denom = cdecz - eta * sdecz
    ra = ranrm ( math.atan2 ( xi, denom ) + raz )
    dc = math.atan2 ( sdecz + (eta * cdecz), math.sqrt( (xi * xi) + (denom * denom) ))
    return ra, dc

deg2rad = math.pi / 180.0
mas2rad = deg2rad / (3600.0 * 1000.0)

def correct(ra, dc, pmra, pmdc, tdelta):
    '''
    Proper motion corrects the coordinates [degrees] by the given proper motion components [mas/yr] and
    time [years].

    Returns proper motion corrected coordinates [degrees].
    '''

    # Offset from a (0,0) origin for local plane coordinates:
    xi = tdelta * pmra * mas2rad
    eta = tdelta * pmdc * mas2rad

    # project back onto celestial sphere given the tangent point (answer in rad)
    rar, dcr =  tangent_plane_to_spherical(xi, eta, ra * deg2rad, dc * deg2rad)

    # return coords in original units (degrees)
    return rar / deg2rad, dcr / deg2rad




#def correct_gband (bp_rp , astrometric_params_solved , phot_g_mean_mag , phot_g_mean_flux ):
def correct_gband (bp_rp , astrometric_params_solvedOrig , phot_g_mean_magOrig):
    """
    Correct the G-band fluxes and magnitudes for the input list of Gaia EDR3 data.
    Parameters
    ----------
    bp_rp: float , array_like
    The (BP-RP) colour listed in the Gaia EDR3 archive.
    astrometric_params_solved: int, array_like
    The astrometric solution type listed in the Gaia EDR3 archive.
    phot_g_mean_mag: float , array_like
    The G-band magnitude as listed in the Gaia EDR3 archive.
    phot_g_mean_flux: float , array_like
    The G-band flux as listed in the Gaia EDR3 archive.
    Returns
    -------
    The corrected G-band magnitudes and fluxes. The corrections are only applied to
    sources with a 6-parameter astrometric solution fainter than G=13, for which a
    (BP-RP) colour is available.
    sources with a 6-parameter astrometric solution fainter than G=13, for which a
    (BP-RP) colour is available.
    Example
    gmag_corr , gflux_corr = correct_gband(bp_rp , astrometric_params_solved ,
    phot_g_mean_mag , phot_g_mean_flux)
    """

    
    if np.isscalar(bp_rp) or np.isscalar(astrometric_params_solvedOrig) or \
                    np.isscalar(phot_g_mean_magOrig): # or np.isscalar(phot_g_mean_flux):
        bp_rp = np.float64(bp_rp)
        astrometric_params_solvedOrig = np.int64(astrometric_params_solvedOrig)
        phot_g_mean_magOrig = np.float64(phot_g_mean_magOrig)
#        phot_g_mean_flux = np.float64(phot_g_mean_flux)

    if not (bp_rp.shape == astrometric_params_solvedOrig.shape \
                        == phot_g_mean_magOrig.shape):# == phot_g_mean_flux.shape):
        raise ValueError('Function parameters must be of the same shape')

    phot_g_mean_mag = np.copy(phot_g_mean_magOrig)
    phot_g_mean_mag[np.isnan(phot_g_mean_mag)]=-999.99999
    astrometric_params_solved=np.copy(np.copy(astrometric_params_solvedOrig))
    astrometric_params_solved[np.isnan(astrometric_params_solved)]=-999
    # added  np.isnan(phot_g_mean_mag)
    print phot_g_mean_mag,bp_rp,astrometric_params_solved
    do_not_correct = np.isnan(bp_rp) | (phot_g_mean_mag<=13) | \
                                       (astrometric_params_solved != 95)
    print do_not_correct,'DNC'
    bright_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>13) & \
                                   (phot_g_mean_mag<=16)
    faint_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>16)
    bp_rp_c = np.clip(bp_rp, 0.25, 3.0)

    correction_factor = np.ones_like(phot_g_mean_mag)
    correction_factor[faint_correct] = 1.00525 - 0.02323*bp_rp_c[faint_correct] + \
        0.01740*np.power(bp_rp_c[faint_correct],2) - \
        0.00253*np.power(bp_rp_c[faint_correct],3)
    correction_factor[bright_correct] = 1.00876 - 0.02540*bp_rp_c[bright_correct] + \
        0.01747*np.power(bp_rp_c[bright_correct],2) - \
        0.00277*np.power(bp_rp_c[bright_correct],3)

    gmag_corrected = phot_g_mean_magOrig - 2.5*np.log10(correction_factor)
#    gflux_corrected = phot_g_mean_flux * correction_factor

    return gmag_corrected #, gflux_corrected


def correct_flux_excess_factor (bp_rpOrig , phot_bp_rp_excess_factor ):
    """
    Calculate the corrected flux excess factor for the input Gaia EDR3 data.
    Parameters
    ----------
    bp_rp: float , array_like
    The (BP-RP) colour listed in the Gaia EDR3 archive.
    phot_bp_rp_flux_excess_factor: float , array_like
    The flux excess factor listed in the Gaia EDR3 archive.
    Returns
    -------
    The corrected value for the flux excess factor , which is zero for "normal" stars.
    Example
    -------
    phot_bp_rp_excess_factor_corr = correct_flux_excess_factor(bp_rp,
    phot_bp_rp_flux_excess_factor)
    """
    bp_rp = np.copy(bp_rpOrig)
    if np.isscalar(bp_rp) or np.isscalar(phot_bp_rp_excess_factor):
        bp_rp = np.float64(bp_rp)
        phot_bp_rp_excess_factor = np.float64(phot_bp_rp_excess_factor)
        
    if bp_rp.shape != phot_bp_rp_excess_factor.shape:
        raise ValueError('Function parameters must be of the same shape!')
    
#    print bp_rp,"SUB"

    do_not_correct = np.isnan(bp_rp)
    bp_rp[np.isnan(bp_rp)]=-999.99999
#    print bp_rp,"SUB"
    bluerange = np.logical_not(do_not_correct) & (bp_rp < 0.5)
    greenrange = np.logical_not(do_not_correct) & (bp_rp >= 0.5) & (bp_rp < 4.0)
    redrange = np.logical_not(do_not_correct) & (bp_rp > 4.0)

    correction = np.zeros_like(bp_rp)
    correction[bluerange] = 1.154360 + 0.033772*bp_rp[bluerange] + \
                                       0.032277*np.power(bp_rp[bluerange],2)
    correction[greenrange] = 1.162004 + 0.011464*bp_rp[greenrange] + \
                                        0.049255*np.power(bp_rp[greenrange],2) \
                                      - 0.005879*np.power(bp_rp[greenrange],3)
    correction[redrange] = 1.057572 + 0.140537*bp_rp[redrange]
    
    return phot_bp_rp_excess_factor - correction


def testF(myvar):
    myvar=myvar+10
    print myvar


def testNP(myArr):
    myArr +=10
    print myArr
