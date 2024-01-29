import pyIGRF
import numpy as np
from brahe.epoch import Epoch
from brahe import frames
from brahe import coordinates
import numpy.linalg as LA
from CubeGNC.utils.transformations import *

def get_magnetic_field(state, epoch):

    # get year
    year = epoch.year()
    r_eci = state[0:3]
    # eci and ecef location
    ecef_Q_eci = frames.rECItoECEF(epoch)
    eci_Q_ecef = np.transpose(ecef_Q_eci)
    # get ecef location
    r_ecef = ecef_Q_eci @ r_eci
    # long lat geod
    longitude, latitude, altitude = coordinates.sECEFtoGEOC(r_ecef, use_degrees='false')
    _, _, _, BN, BE, BD, _ = pyIGRF.igrf_value(latitude, longitude, LA.norm(r_ecef)/1000, year)
    BNED_nT = np.array([[BN], [BE], [-BD]])
    # NED and ECEF DCM
    ecef_Q_ned = ecef_Q_ned_mat(longitude, latitude)

    # convert to eci
    B_eci_nT = eci_Q_ecef @ ecef_Q_ned @ BNED_nT
    # convert from nT to T
    #B_eci_T = B_eci_nT * 1e-9
    return B_eci_nT
