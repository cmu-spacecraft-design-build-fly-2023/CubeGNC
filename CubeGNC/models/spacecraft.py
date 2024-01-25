import numpy as np 
from scipy.linalg import expm

import brahe
from brahe import frames
from brahe.epoch import Epoch
from brahe.orbit_dynamics.gravity import accel_gravity

from CubeGNC.dynamics.drag import *
from CubeGNC.utils.transformations import *
import CubeGNC.dynamics.astrodynamics as astro
from datetime import datetime
from spaceweather import sw_daily

class Spacecraft:

    REQUIRED_KEYS = ["mass", "dt"]
    H = np.vstack([np.zeros(3), np.eye(3)])

    DEFAULTS = {
        "flexible": False,
        "drag": False,
        "third-body": False,
        "solar-radiation": False,
        "gravity-gradient": False,
        "inertia": np.array([[0.0033, 0., 0.],[0., 0.0033, 0.],[0., 0., 0.0033]]),
        "Cd": 2,
        "crossA": 0.01 # m^2
    }

    µ = 3.986004418e14


    def __init__(self, configuration):
        
        # Check if all required keys are provided in the configuration
        missing_keys = [key for key in self.REQUIRED_KEYS if key not in configuration]
        if missing_keys:
            raise ValueError(f"Missing required keys: {', '.join(missing_keys)}")
        

        ## Initial state
        if "initial_attitude" not in configuration:
            raise ValueError(f"Missing required 'initial_attitude'.")
        
        if "initial_orbit_eci" not in configuration and "initial_orbit_oe" not in configuration:
            raise ValueError(f"Missing required initial orbit: 'initial_orbit_eci' or 'initial_orbit_oe'.")
        
        # TODO Docs
        orbit = None
        if "initial_orbit_oe" in configuration:
            orbit = astro.get_CART_from_OSC(configuration["initial_orbit_oe"])
        else:
            orbit = np.array(configuration["initial_orbit_eci"])

        # TODO check sizes 
        self._state = np.concatenate((orbit, np.array(configuration["initial_attitude"])))
        


        ##  Mass
        if configuration["mass"] >= 0.0:
            self._mass = configuration["mass"]
        else:
            raise ValueError("Negative mass value")
        

        ## Timestep
        if configuration["dt"] >= 0.0:
            self._dt = configuration["dt"]
        else:
            raise ValueError("Negative dt value")
        

        ## Inertia 
        if "inertia" in configuration:

            if len(configuration["inertia"]) == 6:

                Iv = configuration["inertia"]
                inertia = np.array([[Iv[0], Iv[3], Iv[4],],[Iv[3], Iv[1], Iv[5]],[Iv[4], Iv[5], Iv[2]]])

                # Check positive-definiteness
                if np.all(np.linalg.eigvals(inertia) > 0):
                    self.J = np.array(inertia)
                    self.invJ = np.linalg.inv(self.J)
                else:
                    raise ValueError("Inertia is not positive-definite.")
                
            else:
                raise ValueError("Need a array of 6 elements to define the inertia: [Ixx, Iyy, Izz, Ixy, Ixz, Iyz].")
            
        else:
            self.J = self.DEFAULTS["inertia"]
            self.invJ = np.linalg.inv(self.J)


        ## TODO default values 
        if "gravity_order" in configuration:
            self._gravity_order = configuration["gravity_order"]
        else:
            self._gravity_order = self.DEFAULTS["gravity_order"]


        if "gravity_degree" in configuration:
            self._gravity_degree = configuration["gravity_degree"]
        else:
            self._gravity_degree = self.DEFAULTS["gravity_degree"]






        ## TODO let the user define 
        ## Found bugs in underlying pysofa in brahe ~
        self.epoch = Epoch(2022,11,26, 12, 0, 5, 0)
        self.epoch_dt = datetime(2022,11,26, 12, 0, 5, 0)

        self.sw = sw_daily()

        if "drag" in configuration:
            self._drag = configuration["drag"]
        else:
            self._drag = self.DEFAULTS["drag"]
        

        # TODO check
        if "crossA" in configuration:
            self._crossA = configuration["crossA"]
        else:
            self._crossA = self.DEFAULTS["crossA"]

        # TODO check
        if "Cd" in configuration:
            self._Cd = configuration["Cd"]
        else:
            self._Cd = self.DEFAULTS["Cd"]



        if "flexible" in configuration:
            self._flexible = configuration["flexible"]
        else:
            self._flexible = self.DEFAULTS["flexible"]



        """
        if "key" in configuration:
            self.key = configuration["key"]
        else:
            self.key = self.DEFAULTS["key"]

        """



    def get_attitude_state(self):
        return self._state[6:13]


    def get_eci_state(self):
        return self._state[0:6]
        


    def get_osculating_state(self):
        """
        1. _a_, Semi-major axis [m]
        2. _e_, Eccentricity [dimensionless]
        3. _i_, Inclination [rad]
        4. _Ω_, Right Ascension of the Ascending Node (RAAN) [rad]
        5. _ω_, Argument of Perigee [ramd]
        6. _M_, Mean anomaly [rad]
        """
        return astro.get_OSC_from_CART(self._state[0:6])


    def get_state(self):
       return self._state


    def orbital_accelerations(self):
        #  x_ecef = R_i2b @ x
        x_eci = self._state
        

        a = np.zeros(3)

        R_i2b = frames.rECItoECEF(self.epoch)
        a += accel_gravity(x_eci[0:6], R_i2b, n_max=self._gravity_degree, m_max=self._gravity_order)
        #a += - (self.µ / np.linalg.norm(x_eci[0:3]) ** 3) * x_eci[0:3]

        if self._drag:
            #r_sun = brahe.sun_position(self.epoch)
            #rho = density_harris_priester(x_eci[0:6], r_sun)
            #print("rho ", rho)
            a += accel_drag(self.epoch_dt, x_eci[0:6], self._mass, self._crossA, self._Cd, R_i2b, self.sw)
            #print("drag", accel_drag(self.epoch_dt, x_eci[0:6], self._mass, self._crossA, self._Cd, R_i2b))


        return a    



    def dynamics(self, x, u):

        # TODO revamp + use of external fct


        xdot = np.zeros(13)
        xdot[0:3] = x[3:6]
        xdot[3:6] = (1 / self._mass) * u[0:3] + self.orbital_accelerations() #- (self.µ / np.linalg.norm(x[0:3]) ** 3) * x[0:3]
        xdot[6:10] = 0.5 * L(x[6:10]) @ self.H @ x[10:13]
        xdot[10:13] = self.invJ.dot(u[3:6] - np.cross(x[10:13], np.dot(self.J, x[10:13])))
        return xdot



    def set_dt(self, new_dt):
        if new_dt <= 0:
            raise ValueError("dt must be positive.")
        self.dt = new_dt
        return True

              
    def advance(self, u=np.zeros(6)):

        # TODO check u size for all cases

        q = self._state[6:10]
        ω = self._state[10:13]

        k1 = self._dt * self.dynamics(self._state, u)
        k2 = self._dt * self.dynamics(self._state + k1 * 0.5, u)
        k3 = self._dt * self.dynamics(self._state + k2 * 0.5, u)
        k4 = self._dt * self.dynamics(self._state + k3, u)

        self._state = self._state + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        self._state[6:10] = np.dot(expm(R(self.H.dot(0.5 * self._dt * ω + (self._dt/6) * (k1[10:13] + k2[10:13] + k3[10:13])))), q)






# Temporary local testing
if __name__ == "__main__":


    config = {
        "mass": 2.0,
        "inertia": [10,20,30,0.0,0.0,0.0],
        "dt": 1.0,
        "flexible": False,
        "initial_attitude":[1.0,0,0,0,0.1,0.1,0.1],
        "initial_orbit_oe":[1.5e6, 0, 0, 0, 0, 0],
        "gravity_order": 5,
        "gravity_degree": 5,
        "drag": True
    }



    spacecraft = Spacecraft(config)
    print(spacecraft.J)
    print(spacecraft.get_state())
    for i in range(10):
        spacecraft.advance()
        print(spacecraft.get_state()[0:6])

    print(spacecraft.epoch)
