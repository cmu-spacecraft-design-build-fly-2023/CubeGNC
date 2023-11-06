import numpy as np
import brahe


from CubeGNC.utils.transformations import *


def apply_SO3_noise(vec, std):
    # Axis-angle vector noise 
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

    def __init__(self, std, bias):
        self.std = std
        self.bias = bias
        # TODO bias dynamics

    def measure(self, spacecraft):
        pass
    

    
class Magnetometer():

    def __init__(self):
        pass

    def measure(self, spacecraft):
        pass
    

class Vision():

    def __init__(self):
        pass

    def measure(self, spacecraft):
        pass
    
