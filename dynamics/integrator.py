import numpy as np
from scipy.linalg import expm

# Runge-Kutta 4
def rk4(dynamics, x, u, dt):
    k1 = dt * dynamics(x, u)
    k2 = dt * dynamics(x + k1 * 0.5, u)
    k3 = dt * dynamics(x + k2 * 0.5, u)
    k4 = dt * dynamics(x + k3, u)
    x = x + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
    return x



# hermite simpson implicit integrator residual
def hermite_simpson(dynamics, x1, x2, u, dt):
    x_k12 = 0.5 * (x1 + x2) + (dt / 8) * (dynamics(x1, u) - dynamics(x2, u))
    return x1 + (dt / 6) * (dynamics(x1, u) + 4 * dynamics(x_k12, u) + dynamics(x2, u)) - x2



# implicit midpoint integrator residual
def implicit_midpoint(dynamics, x1, x2, u, dt):
    return x1 + dt * dynamics(0.5 * (x1 + x2), u) - x2




def rk4mk_attitude(params, dynamics, x1, u, dt):
    H = np.vstack([np.zeros(3), np.eye(3)])
    k1 = dt * dynamics(params, x1, u)
    k2 = dt * dynamics(params, x1 + k1 * 0.5, u)
    k3 = dt * dynamics(params, x1 + k2 * 0.5, u)
    k4 = dt * dynamics(params, x1 + k3, u)

    f = x1 + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)

    q = x1[0:4]
    ω = x1[4:7]

    # Using expm from scipy to exponentiate a matrix
    f[0:4] = np.dot(expm(H.dot(0.5 * dt * ω + (dt/6) * (k1[4:7] + k2[4:7] + k3[4:7]))), q)

    return f