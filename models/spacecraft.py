import numpy as np 


class Spacecraft:

    def __init__(self, configuration):
        # TODO handle configuration
        pass


    def advance(self, u, dt):
        # Need three versions? (orbit, attitude, both)
        pass



"""

config =    {
    "mode" = "orbit-only", "attitude-only", "both",
    "mass" = 0.5,
    "inertia" = [],
    "fleible" = ..,
    "dt" = ....,
    "gyroscope" = True,
    }

std_noise = ... 
bias = 0.1 
Gyroscope(std_noise, bias)
spacecraft = Cubesat(config)

params["inertia"]

"""