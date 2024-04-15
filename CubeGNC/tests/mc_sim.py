import threading
import multiprocessing
import numpy as np
from time import sleep
from CubeGNC.models.spacecraft import Spacecraft
from CubeGNC.visualizer.mc import Visualizer
from CubeGNC.models.sensors import *
from CubeGNC.models.actuators import *
from CubeGNC.controller.bcross import *
from CubeGNC.dynamics.magnetic_field import *
from brahe.orbit_dynamics.gravity import accel_gravity, accel_thirdbody_sun, accel_thirdbody_moon

from CubeGNC.dynamics.drag import *
from CubeGNC.dynamics.magnetic_field import *
from CubeGNC.tests.mc_sim_initialization import *
import numpy as np
import matplotlib.pyplot as plt

def mc_sim_run(spacecraft, bcross_controller, Ntrials = 1, tspan = (0.0, 2 * 60 * 60), integrator_dt=0.1, controller_dt=1.0):
    """
    Performs a Monte Carlo simulation of orbit and attitude dynamics.

    Args:
        get_initial_state: A function that generates a random initial state.
        controllers: A dictionary of controller functions, where keys are controller names and values are the controller functions.
        Ntrials: The number of Monte Carlo trials.
        params: Orbit dynamics parameters.
        tspan: The time span for the simulation.
        integrator_dt: The time step for the integrator.
        controller_dt: The time step for the controller.

    Returns:
        A dictionary of Monte Carlo data, where keys are controller names and values are dictionaries containing simulation results.
    """

    Ntimesteps = int(np.ceil((tspan[1] - tspan[0]) / integrator_dt))
    mc_data = {
            "X": np.zeros((Ntrials, 13, Ntimesteps)),
            "U": np.zeros((Ntrials, 6, Ntimesteps)),
            "T": np.zeros((Ntrials, 1, Ntimesteps))
    }

    # Check for available threads and warn if necessary
    num_threads = multiprocessing.cpu_count()
    print('n_threads', num_threads)
    if threading.active_count() < num_threads:
        print(f"Warning: Consider using more threads for optimal performance. Start Python with --threads {num_threads} to max out your CPU.")

    # Create threads for parallel execution
    threads = []
    for mc_step in range(Ntrials):
        thread = threading.Thread(target=run_trial, args=(spacecraft, bcross_controller, mc_step, tspan, integrator_dt, controller_dt, mc_data))
        threads.append(thread)
        thread.start()

    # Wait for all threads to finish
    for thread in threads:
        thread.join()

    return mc_data

def run_trial(spacecraft, bcross_controller, mc_step, tspan, integrator_dt, controller_dt, mc_data):
    """
    Runs a single Monte Carlo trial.
    """
    x0 = mc_setup_get_initial_state(
    (500e3, 600e3),  # h_range
    (0.0, 0.0),      # e_range
    (deg2rad(20), deg2rad(160)),  # i_range
    (0, 2 * np.pi),   # Ω_range
    (0, 2 * np.pi),   # ω_range
    (0, 2 * np.pi),   # M_range
    (deg2rad(30), deg2rad(30)),  # angular_rate_magnitude_range
    )
    print(f"Thread {threading.get_ident()}, Trial {mc_step + 1}: x0 = {x0}")

    xhist, uhist, thist = simulate_satellite_orbit_attitude_rk4(x0, spacecraft, bcross_controller, tspan = tspan, integrator_dt=integrator_dt, controller_dt=controller_dt)
    mc_data["X"][mc_step, :, :] = xhist
    mc_data["U"][mc_step, :, :] = uhist
    mc_data["T"][mc_step, 0, :] = thist

def simulate_satellite_orbit_attitude_rk4(x0, spacecraft, bcross_controller, tspan, integrator_dt=0.1, controller_dt=1.0):
    """Simulates satellite orbit and attitude using a 4th-order Runge-Kutta integrator.

    Args:
        x0: Initial state vector (numpy.ndarray)
        params: OrbitDynamicsParameters object
        tspan: Time span tuple (start time, end time)
        integrator_dt: Integrator time step (optional, default=0.1)
        controller: Controller function (optional, default=null_controller)
        controller_dt: Controller update time step (optional, default=1.0)

    Returns:
        xhist: State history (numpy.ndarray)
        uhist: Control history (numpy.ndarray)
        thist: Time history (numpy.ndarray)
    """
    #print('x0', x0)
    r0 = x0[:3]  # ECI frame position
    v0 = x0[3:6]  # Inertial frame velocity
    q0 = x0[6:10]  # Body to ECI quaternion
    omega0 = x0[10:13]  # Body frame angular rate

    z0 = np.concatenate([
        r0,
        v0,
        q0,
        omega0
    ])

    Nt = int(np.ceil((tspan[1] - tspan[0]) / integrator_dt))

    thist = np.zeros(Nt)
    zhist = np.zeros((13, Nt))
    uhist = np.zeros((6, Nt))

    zhist[:, 0] = z0
    #print('z0', z0)
    next_controller_update_time = tspan[0]
    spacecraft.set_state(z0)

    for k in range(1, Nt-1):
        #print('k',k)

        t = (integrator_dt * (k - 1) + tspan[0])
        #print('t',t)
        thist[k] = t

        if next_controller_update_time <= t:
            B_eci_nT, B_body = get_magnetic_field(spacecraft.get_state(), spacecraft.get_epoch())
            #print('B_body', B_body)
            uhist[3:, k] = bcross_controller.get_control_dipole(spacecraft.get_state()[10:13], B_body.T)
            #print('u', uhist[:, k].shape)
            next_controller_update_time = t + controller_dt
        else:
            uhist[:, k] = uhist[:, k - 1]
        #print('uhist[:, k]', uhist[:, k])
        znext = spacecraft.advance(uhist[:, k], B_body.T.squeeze())
        #print('znext',znext)
        znext[6:10] /= np.linalg.norm(znext[6:10])  # Normalize quaternion
        zhist[:, k + 1] = znext[:13]
        spacecraft.set_state(znext[:13])

    xhist = np.concatenate([
        zhist[:3, :],
        zhist[3:6, :],
        zhist[6:10, :],
        zhist[10:13, :]
    ], axis=0)

    return xhist, uhist, thist

def get_downsample(Ntrials, max_samples):
    """
    Returns a boolean array indicating which elements to keep for downsampling.

    Args:
        Ntrials (int): Total number of trials.
        max_samples (int): Maximum number of samples to keep.

    Returns:
        numpy.ndarray: Boolean array of length Ntrials, True for elements to keep.
    """

    sample_steps = int(np.ceil(Ntrials / max_samples))
    indices = np.arange(Ntrials)
    downsample = (indices % sample_steps == 0) | (indices == Ntrials - 1)  # Include the last element
    return downsample

def mc_plot_momentum_magnitude_vs_time(mc_results, spacecraft, max_samples=500, file_suffix=""):
    """
    Plots the momentum magnitude vs time for multiple controllers based on Monte Carlo results.

    Args:
        mc_results (dict): Dictionary containing Monte Carlo results, with controller names as keys.
        params (dict): Dictionary containing parameters, including "satellite_model" with "inertia".
        max_samples (int, optional): Maximum number of samples to plot. Defaults to 500.
        file_suffix (str, optional): Suffix to append to saved plot filenames. Defaults to "".
    """

    J = spacecraft.J
    Ntrials = mc_results["T"].shape[0]
    print ('J',J)

    h_average = np.zeros((max_samples + 1,))
    t_plot_average = np.zeros_like(h_average)

    for mc_step in range(Ntrials):
            xhist = mc_results["X"][mc_step, :, :]
            uhist = mc_results["U"][mc_step, :, :]
            thist = mc_results["T"][mc_step, 0, :]
            #downsample = get_downsample(len(thist), max_samples)

            #downsample = np.linspace(0, len(thist) - 1, max_samples + 1).astype(int)  # Even downsamplin]
            omega = xhist[10:13]
            #omega = xhist[10:13, downsample]
            h = J @ omega
            h_mag = np.linalg.norm(h, axis=0)
            #print('omega', omega.shape)
            #print('h', h.shape)
            #print('h_mag', h_mag.shape)
            #t_plot = thist[downsample] / (60 * 60)

            plt.plot(h_mag.T[2:], alpha=0.4, linewidth=2)

            #h_average += h_mag
            #t_plot_average += t_plot

    #h_average /= Ntrials

    # Plot average with a dotted line
    #plt.plot(h_average, color="black", linestyle="dotted", linewidth=2, label="Average")

    plt.xlabel("Time (seconds)")
    plt.ylabel(r"$\|h\|$ (Nms)")
    plt.title("Bcross_control")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig(f"mc_plot_bcross_test10.pdf")
    plt.close()  # Close each plot to avoid figure accumulation


if __name__ == "__main__":
    config = {
        "inertia": [0.0046,0.0046, 0.0046, 0.00002,-0.000003,-0.00002],
        "dt": 1,
        "flexible": False,
        "initial_attitude":[1.0,0,0,0,0.5,-0.02,0.03],
        "initial_orbit_oe":[5.3e6, 0, 0, 0, 0, 0],
        "mass": 1.7,  # kg,
        "gravity_order": 5,
        "gravity_degree": 5,
        "drag": True,
        "third_body": True
    }
    spacecraft = Spacecraft(config)
    magnetometer = Magnetometer(0.0, 5.0)
    magnetorquer = Magnetorquer(10,  0.6e-3, 0.036e-3, 315)
    bcross_controller = Controller(spacecraft, magnetorquer)

    #params = get_mc_params(params)

    mc_results = mc_sim_run(spacecraft, bcross_controller, Ntrials=10, tspan=(0.0, 2*60*60), integrator_dt=0.1, controller_dt=1.0)
    mc_plot_momentum_magnitude_vs_time(mc_results, spacecraft)
    print(mc_results["X"].shape)
    import yaml

    with open("my_data_test10.yaml", "w") as f:
        yaml.dump(mc_results["X"], f)
    import pickle

    with open("my_data_test10.pickle", "wb") as f:
        # Write the dictionary to the file using pickle.dump()
        pickle.dump(mc_results, f)


