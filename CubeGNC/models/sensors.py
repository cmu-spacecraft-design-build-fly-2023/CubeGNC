import numpy as np
import brahe
from math import radians

from CubeGNC.utils.transformations import *


def apply_SO3_noise(vec, std):
    """
    Adds S03-style noise to a vector by rotating it about a random axis.

    Args:
        vec: The vector to which noise will be added.
        std: The standard deviation of the noise.

    Returns:
        The noisy vector.
    """
    noise = std * np.random.randn(3)
    return dcm_from_phi(noise) @ vec


class SunSensor():

    def __init__(self, std_deg, offset=None):
        self.std = np.deg2rad(std_deg)
        self.offset = offset

    def measure(self, spacecraft):
        r_sun = brahe.sun_position(spacecraft.epoch)
        vec_sun_body = self.get_sun_vec_body(spacecraft.get_eci_state(), r_sun, spacecraft.get_attitude_state())
        return apply_SO3_noise(vec_sun_body, self.std)


    def get_sun_vec_body(self, x_eci, sun_eci, attitude_q):
        # Vector from the spacecraft to the Sun in ECI frame
        vec_sun_eci = sun_eci - x_eci
        # Vector from the spacecraft to the Sun in body frame
        vec_sun_body = dcm_from_q(attitude_q) @ (vec_sun_eci / np.linalg.norm(vec_sun_eci)) # normalized
        return vec_sun_body
    

    ## TODO Flux values 





class GPS():

    def __init__(self, std):
        self.std = std

    def measure(self, spacecraft):
        pass




class Gyroscope():

    def __init__(self, rotate_std_deg, noise_std_degps, initial_bias_deg):
        # TODO bias dynamics
        self.offset = dcm_from_phi(radians(rotate_std_deg) * np.random.randn(3))
        self.bias = radians(initial_bias_deg) * np.random.randn(3) / np.linalg.norm(np.random.randn(3))
        self.noise = radians(noise_std_degps) * np.random.randn(3)

    def measure(self, spacecraft):
        measured_value = self.offset * spacecraft.get_state()[10:13] + self.bias + self.noise
        return measured_value
    

    
class Magnetometer():
    """
    def __init__(self, rotate_std_deg, noise_std_deg):
        self.offset = dcm_from_phi(radians(rotate_std_deg) * np.random.randn(3))
        self.noise_std = radians(noise_std_deg)

    def measure(self, spacecraft):
        true_B_body = ...
        measured_value = apply_SO3_noise(self.offset * true_B_body, self.noise_std)
        return measured_value
    """
    

class Vision():

    def __init__(self):
        pass

    def measure(self, spacecraft):
        pass
    
