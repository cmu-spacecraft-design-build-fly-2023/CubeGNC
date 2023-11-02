import numpy as np 


class Spacecraft:

    REQUIRED_KEYS = ["mass", "mode", "dt"]


    DEFAULTS = {
        "flexible": True,
        "drag": False,
        "third-body": False,
        "solar-radiation": False,
        "gravity-gradient": False,
        "gyroscope": None,
        "sun-sensor": None,
        "gps": None,
    }


    def __init__(self, configuration):
        
        # Check if all required keys are provided in the configuration
        missing_keys = [key for key in self.REQUIRED_KEYS if key not in configuration]
        if missing_keys:
            raise ValueError(f"Missing required keys: {', '.join(missing_keys)}")
        

        # Mode
        if "mode" in configuration:
            if configuration["mode"] in {"attitude", "orbit", "both"}:
                self.mode = configuration["mode"]
            else:
                raise ValueError("Invalid mode provided. 'mode' must be 'attitude', 'orbit', or 'both'.")

        # Mass
        if configuration["mass"] >= 0.0:
            self.mass = configuration["mass"]
        else:
            raise ValueError("Negative mass value")
        
        # Timestep
        if configuration["dt"] >= 0.0:
            self.dt = configuration["dt"]
        else:
            raise ValueError("Negative dt value")
        

        # Default values 
        if "inertia" in configuration:

            # isinstance(configuration["inertia"], np.ndarray) and 
            if len(configuration["inertia"]) == 6:

                Iv = configuration["inertia"]
                inertia = np.array([[Iv[0], Iv[3], Iv[4],],[Iv[3], Iv[1], Iv[5]],[Iv[4], Iv[5], Iv[2]]])

                # Check positive-definiteness
                if np.all(np.linalg.eigvals(inertia) > 0):
                    self.inertia = np.array(inertia)
                else:
                    raise ValueError("Inertia is not positive-definite.")
                
            else:
                raise ValueError("Need a array of 6 elements to define the inertia: [Ixx, Iyy, Izz, Ixy, Ixz, Iyz].")
            
        else:
            self.inertia = self.DEFAULTS["inertia"]




        if "flexible" in configuration:
            self.flexible = configuration["flexible"]
        else:
            self.flexible = self.DEFAULTS["flexible"]


        if "gyroscope" in configuration:
            self.gyroscope = configuration["gyroscope"]
        else:
            self.gyroscope = self.DEFAULTS["gyroscope"]


        """
        if "key" in configuration:
            self.key = configuration["key"]
        else:
            self.key = self.DEFAULTS["key"]

        """

    def get_current_state():
        pass

              
    def advance(self, u, dt):
        # Need three versions? (orbit, attitude, both)
        pass





# Temporary local testing
if __name__ == "__main__":


    config = {
        "mode": "both",
        "mass": 1500.0,
        "inertia": [10,10,10,0.0,0.0,0.0],
        "dt": 1.0,
        "flexible": False,
        #"gyroscope": True
    }


    spacecraft = Spacecraft(config)
    print(spacecraft.inertia)

    """
    std_noise = ... 
    bias = 0.1 
    Gyroscope(std_noise, bias 
    """
