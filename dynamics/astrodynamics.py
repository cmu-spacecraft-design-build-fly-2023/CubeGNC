import numpy as np
import brahe




def keplerian_central_body(params, x, u):
    r = x[0:3] 
    r_ddot = (1 / params['m']) * u - (params['mu'] / np.linalg.norm(r) ** 3) * r
    return np.hstack((x[3:6], r_ddot))




def semi_major_axis(rp, ra):
    """
    rp: periapsis -> radius from the central body to the nearest point on the orbital path
    ra: apoapsis -> radius from the central to the farthest point on the orbital path
    semi major axis: longest semi-diameter of an ellipse. 
    """
    return (rp + ra) / 2

def orbital_period(sma, mu):
    """
    sma: semi major axis -> longest semi-diameter of an ellipse.
    mu: standard gravitational parameter of the central body
    """
    return 2.0 * np.pi * np.sqrt((sma ** 3) / mu)

def eccentricity(rp, ra):
    """
    rp: periapsis -> radius from the central body to the nearest point on the orbital path
    ra: apoapsis -> radius from the central to the farthest point on the orbital path
    """
    semi_major_axis_val = semi_major_axis(rp, ra)
    ecc = 1 - rp / semi_major_axis_val
    return ecc




def get_CART_from_OSC(x_oe, degree=False):
    """
    Return the cartesian state (position and velocity, ECI) given the corresponding set of osculating orbital elements.
    The vector must contain in order (SatelliteDynamics.jl):
        1. a, Semi-major axis [m]
        2. e, Eccentricity [dimensionless]
        3. i, Inclination [rad]
        4. Ω, Right Ascension of the Ascending Node (RAAN) [rad]
        5. ω, Argument of Perigee [rad]
        6. M, Mean anomaly [rad]
    """
    return brahe.sOSCtoCART(x_oe, degrees=degree)




def accel_gravity():
    pass

def accel_J2():
    pass

