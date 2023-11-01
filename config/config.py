import numpy as np
import yaml

def config_parser(path_to_config_yaml):
    """This function parses the config YAML for the sim"""

    # Load config YAML
    with open(path_to_config_yaml, 'r') as yaml_file:
        data = yaml.safe_load(yaml_file)

    sc_properties = data["sc_properties"]
    gravity = data["gravity"]
    environment = data["environment"]
    init_oe = data["initial_orbital_elements"]
    propagator_parameters = data["propagator_parameters"]
    initial_attitude_conditions = data["initial_attitude_conditions"]
    sensor_properties = data["sensor_properties"]
    actuator_properties = data["actuator_properties"]

    # Inertia matrix
    J = np.array([sc_properties["inertia_rows"][0],
                  sc_properties["inertia_rows"][1],
                  sc_properties["inertia_rows"][2]])
    
    # Inertia Inverse
    invJ = np.linalg.inv(J)

    # Spacecraft properties
    sc_mass = sc_properties["mass"]
    sc_area = sc_properties["area"]
    sc_drag_coeff = sc_properties["drag_coeff"]

    # Spacecraft properties dictionary
    spacecraft = {
        "mass": sc_mass,
        "area": sc_area,
        "cd": sc_drag_coeff,
        "J": J,
        "invJ": invJ
    }

    # Gravity params
    gravity_params = {
            "spherical_harmonic_gravity_bool": gravity["spherical_harmonic_gravity_bool"],
            "grav_deg": gravity["gravity_degree"],
            "grav_order": gravity["gravity_order"]
        }

    # Timing related stuff
    epoch_orbital = propagator_parameters["epoch"]
    dt_orbital = propagator_parameters["dt_orbital"]
    dt_attitude = propagator_parameters["dt_attitude"]
    dt_controller = propagator_parameters["dt_controller"]
    tf = propagator_parameters["time_final"] * 24 * 3600  # seconds

    # Time vectors for simulation
    t_vec_orbital = np.arange(0, tf + dt_orbital, dt_orbital)
    inner_loop_t_vec = np.arange(0, dt_attitude, dt_orbital)
    t_vec_attitude = np.arange(0, t_vec_orbital[-1] + dt_attitude, dt_attitude)

    # Time parameters dictionary
    time_params = {
        "dt_orbital": dt_orbital,
        "dt_attitude": dt_attitude,
        "dt_controller": dt_controller,
        "tf": tf,
        "t_vec_orbital": t_vec_orbital,
        "inner_loop_t_vec": inner_loop_t_vec,
        "t_vec_attitude": t_vec_attitude
    }

    # Initial conditions
    NqB_0 = initial_attitude_conditions["q0"]
    w_0 = initial_attitude_conditions["w0"]

    # Orbital elements
    sma = init_oe["sma"]
    ecc = init_oe["ecc"]
    inc = init_oe["inc"]
    raan = init_oe["raan"]
    argp = init_oe["argp"]
    M = init_oe["M"]

    # Initial orbital elements
    oe_0 = [sma, ecc, inc, raan, argp, M]

    # Convert to ECI r and v
    #eci_rv_0 = sOSCtoCART(oe_0, use_degrees=False) [Need brahe]

    initial_conditions = {
        "epoch_orbital": epoch_orbital,
        "NqB0": NqB_0,
        "w0": w_0,
        "oe0": oe_0
    }

    gyro_properties = sensor_properties["gyro"]
    sun_sensor_properties = sensor_properties["sun_sensor"]
    magnetometer_properties = sensor_properties["magnetometer"]

    gyro = {
        "rotate_std_deg": gyro_properties["rotate_std_deg"],
        "noise_std_degps": gyro_properties["noise_std_degps"],
        "initial_bias_deg": gyro_properties["initial_bias_deg"],
        "bias_noise_std": gyro_properties["bias_noise_std"]
    }

    sun_sensor = {
        "rotate_std_deg": sun_sensor_properties["rotate_std_deg"],
        "noise_std_deg": sun_sensor_properties["noise_std_deg"]
    }

    magnetometer = {
        "rotate_std_deg": magnetometer_properties["rotate_std_deg"],
        "noise_std_deg": magnetometer_properties["noise_std_deg"]
    }

    offsets = {
        "gyro_bias": np.deg2rad(gyro["initial_bias_deg"]) * np.random.randn(3)
    }

    sensors = {
        "gyro": gyro,
        "sun_sensor": sun_sensor,
        "magnetometer": magnetometer,
        "offsets": offsets
    }

    magnetorquer_properties = actuator_properties["magnetorquer"]
    
    magnetorquer = {
        "max_dipole_strength": magnetorquer_properties["max_dipole_strength"]
    }

    actuators = {
        "magnetorquer": magnetorquer
    }

    params = {
        "sc": spacecraft,
        "time_params": time_params,
        "initial_conditions": initial_conditions,
        "sensors": sensors,
        "actuators": actuators,
        "gravity_params": gravity_params
    }

    return params


if __name__ == "__main__":
    params = config_parser('test_config.yml')
    print(params)