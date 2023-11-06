import numpy as np
from CubeGNC.utils.transformations import *


def quat_kinematics(q, ω):
    # TODO - use other formulation to remove unnecesary allocation
    H = np.vstack([np.zeros(3), np.eye(3)])
    q̇ = 0.5 * L(q) @ H @ ω
    return q̇

def euler_rotational_dynamics(params, ω, τ):
    return np.linalg.inv(params.J).dot(τ - np.cross(ω, np.dot(params.J, ω)))

def attitude_dynamics(params, rot_state):
    xdot = np.zeros(7)
    xdot[0:4] = quat_kinematics(rot_state[0:4])
    xdot[4:7] = euler_rotational_dynamics(params, rot_state[4:7])
    return xdot