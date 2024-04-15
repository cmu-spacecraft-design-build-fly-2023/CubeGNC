import numpy as np
from numpy.random import uniform, random, rand
from scipy.spatial.transform import Rotation
import CubeGNC.dynamics.astrodynamics as astro
import brahe
def mc_initial_orbital_elements(h_range, e_range, i_range, Ω_range, ω_range, M_range):
    """
    Generates a set of random orbital elements within specified ranges.

    Args:
        h_range: Tuple of (minimum altitude, maximum altitude) in meters.
        e_range: Tuple of (minimum eccentricity, maximum eccentricity).
        i_range: Tuple of (minimum inclination, maximum inclination) in degrees.
        Ω_range: Tuple of (minimum longitude of ascending node, maximum longitude of ascending node) in degrees.
        ω_range: Tuple of (minimum argument of perigee, maximum argument of perigee) in degrees.
        M_range: Tuple of (minimum mean anomaly, maximum mean anomaly) in degrees.

    Returns:
        List of orbital elements: [altitude, eccentricity, inclination, longitude of ascending node, argument of perigee, mean anomaly]
    """

    earth_radius = 6371000  # Replace with appropriate value if needed

    return [
        uniform(h_range[0], h_range[1]) + earth_radius,
        uniform(e_range[0], e_range[1]),
        np.deg2rad(uniform(i_range[0], i_range[1])),  # Convert inclination to radians
        np.deg2rad(uniform(Ω_range[0], Ω_range[1])),  # Convert longitude of ascending node to radians
        np.deg2rad(uniform(ω_range[0], ω_range[1])),  # Convert argument of perigee to radians
        np.deg2rad(uniform(M_range[0], M_range[1]))   # Convert mean anomaly to radians
    ]

def mc_initial_attitude():
    """
    Generates a random initial attitude quaternion.

    Returns:
        NumPy array representing a quaternion.
    """

    ϕ = random(3)
    ϕ /= np.linalg.norm(ϕ)
    θ = rand() * 2 * np.pi
    rotation_vector = θ * ϕ

    return Rotation.from_rotvec(rotation_vector).as_quat()

def mc_initial_angular_velocity(ω_magnitude_range):
    """
    Generates a random initial angular velocity vector.

    Args:
        ω_magnitude_range: Tuple of (minimum angular velocity magnitude, maximum angular velocity magnitude) in radians/second.

    Returns:
        NumPy array representing the angular velocity vector.
    """

    ω_magnitude = uniform(ω_magnitude_range[0], ω_magnitude_range[1])
    ω_direction = random(3)
    ω_direction /= np.linalg.norm(ω_direction)

    return ω_magnitude * ω_direction

def state_from_osc(xosc, q, ω):
    """
    Computes a state vector from orbital elements, attitude, and angular velocity.

    Args:
        xosc: Vector of orbital elements (altitude, eccentricity, inclination, longitude of ascending node, argument of perigee, mean anomaly) in radians.
        q: Vector representing the attitude quaternion.
        ω: Vector representing the angular velocity.

    Returns:
        NumPy array of the state vector.
    """
    rv = astro.get_CART_from_OSC(xosc, degrees=False)
    print('rv', rv)
    x = np.concatenate([rv[:3], rv[3:], q, ω])
    return x

def mc_setup_get_initial_state(h_range, e_range, i_range, Ω_range, ω_range, M_range, angular_rate_magnitude_range):
    """
    Sets up a function that generates random initial states for a satellite.

    Args:
        h_range: Tuple of (minimum altitude, maximum altitude) in meters.
        e_range: Tuple of (minimum eccentricity, maximum eccentricity).
        i_range: Tuple of (minimum inclination, maximum inclination) in degrees.
        Ω_range: Tuple of (minimum longitude of ascending node, maximum longitude of ascending node) in degrees.
        ω_range: Tuple of (minimum argument of perigee, maximum argument of perigee) in degrees.
        M_range: Tuple of (minimum mean anomaly, maximum mean anomaly) in degrees.
        angular_rate_magnitude_range: Tuple of (minimum angular velocity magnitude, maximum angular velocity magnitude) in radians/second.

    Returns:
        A function that, when called, returns a random initial state.
    """

    x0_osc = mc_initial_orbital_elements(h_range, e_range, i_range, Ω_range, ω_range, M_range)
    print('x0_osc', x0_osc)
    x0_osc[3] = brahe.sun_sync_inclination(x0_osc[0], x0_osc[1], use_degrees=False)
    q0 = mc_initial_attitude()
    ω0 = mc_initial_angular_velocity(angular_rate_magnitude_range)
    return state_from_osc(x0_osc, q0, ω0)

from numpy import deg2rad


if __name__ == "__main__":
    get_initial_state = mc_setup_get_initial_state(
        (450e3, 700e3),  # h_range
        (0.0, 0.0),  # e_range
        (deg2rad(20), deg2rad(160)),  # i_range
        (0, 2 * np.pi),  # Ω_range
        (0, 2 * np.pi),  # ω_range
        (0, 2 * np.pi),  # M_range
        (deg2rad(30), deg2rad(30)),  # angular_rate_magnitude_range
    )

    x0 = get_initial_state()
    print(x0)