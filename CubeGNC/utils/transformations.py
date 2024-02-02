import numpy as np
from scipy.linalg import expm

def skew_symmetric(w):
    """
    Returns the skew symmetric form of a numpy array.
    w --> [w]
    """
    
    return np.array([[0, -w[2], w[1]], 
                     [w[2], 0, -w[0]], 
                     [-w[1], w[0], 0]])


### Quaternions 
# https://roboticexplorationlab.org/papers/planning_with_attitude.pdf

def L(q):
    """
    Left-multiply
    """
    L = np.zeros((4,4))
    L[0,0] = q[0]
    L[0,1:] = -q[1:].T
    L[1:,0] = q[1:]
    L[1:,1:] = q[0]*np.eye(3) + skew_symmetric(q[1:])
    return L

def R(q):
    """
    Right-multiply
    """
    R = np.zeros((4,4))
    R[0,0] = q[0]
    R[0,1:] = -q[1:].T
    R[1:,0] = q[1:]
    R[1:,1:] = q[0]*np.eye(3) - skew_symmetric(q[1:])
    return R

def conj(q):
    """
    Inverse of a unit quaternion is its conjugate, i.e. same quaternion with a negated vector part 
    """
    return np.hstack((q[0], -q[1:]))


def rotm_to_quat(r):
    q = np.zeros(4)
    q[0] = 0.5 * np.sqrt(1 + r[0,0] + r[1,1] + r[2,2])
    q[1] = (1/(4*q[0])) * (r[2][1] - r[1][2])
    q[2] = (1/(4*q[0])) * (r[0][2] - r[2][0])
    q[3] = (1/(4*q[0])) * (r[1][0] - r[0][1])
    return np.array(q)


def dcm_from_q(q):
    norm = np.linalg.norm(q)
    q0, q1, q2, q3 = q / norm if norm != 0 else q

    # DCM
    Q = np.array([
        [2*q1**2 + 2*q0**2 - 1,   2*(q1*q2 - q3*q0),   2*(q1*q3 + q2*q0)],
        [2*(q1*q2 + q3*q0),       2*q2**2 + 2*q0**2 - 1,   2*(q2*q3 - q1*q0)],
        [2*(q1*q3 - q2*q0),       2*(q2*q3 + q1*q0),   2*q3**2 + 2*q0**2 - 1]
    ])

    return Q

def quat_to_axisangle(q):
    """
    quat wxyz
    """
    axis = np.zeros(3)
    angle = 2 * np.arccos(q[0])
    axis = q[1:]/np.sqrt(1 - q[0]*q[0])
    return axis*angle


def dcm_from_phi(Φ):
    """Compute DCM from an axis-angle"""
    return expm(skew_symmetric(Φ))

def ecef_Q_ned_mat(longitude, latitude):
    phi = np.radians(latitude)
    lam = np.radians(longitude)

    ecef_Q_ned = np.array([[-np.sin(phi) * np.cos(lam), -np.sin(lam), -np.cos(phi) * np.cos(lam)],
                           [-np.sin(phi) * np.sin(lam), np.cos(lam), -np.cos(phi) * np.sin(lam)],
                           [np.cos(phi), 0.0, -np.sin(phi)]])

    return ecef_Q_ned

