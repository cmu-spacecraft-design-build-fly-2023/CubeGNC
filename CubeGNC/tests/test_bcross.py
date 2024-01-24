
from CubeGNC.controller.bcross import Controller

import numpy as np
import unittest


class TestBcrossGainFunction(unittest.TestCase):

    def test_typical_values_1(self):
        R = 6.371e6 #meters
        altitude = 525*1000 #meters
        x0 = R + altitude
        y0 = 0
        z0 = 0
        semi_major_axis = np.linalg.norm([x0,y0,z0])
        print(semi_major_axis)
        inclination = np.pi/4  # degrees
        min_moment_of_inertia = 0.1 # kg*m^2
        expected_gain = 7515.06219
        actual_gain = Controller._compute_bcross_gain(semi_major_axis, inclination, min_moment_of_inertia)
        self.assertAlmostEqual(actual_gain, expected_gain, places=3)
    
    def test_typical_values_2(self):
        R = 6.371e6 #meters
        altitude = 450*1000 #meters
        x0 = R + altitude
        y0 = 0
        z0 = 0
        semi_major_axis = np.linalg.norm([x0,y0,z0])
        print(semi_major_axis)
        inclination = np.pi/4  # degrees
        min_moment_of_inertia = 0.1 # kg*m^2
        expected_gain = 7639.349585
        actual_gain = Controller._compute_bcross_gain(semi_major_axis, inclination, min_moment_of_inertia)
        self.assertAlmostEqual(actual_gain, expected_gain, places=3)

    def test_typical_values_3(self):
        R = 6.371e6 #meters
        altitude = 700*1000 #meters
        x0 = R + altitude
        y0 = 0
        z0 = 0
        semi_major_axis = np.linalg.norm([x0,y0,z0])
        print(semi_major_axis)
        inclination = np.pi/4  # degrees
        min_moment_of_inertia = 0.1 # kg*m^2
        expected_gain = 7237.810407
        actual_gain = Controller._compute_bcross_gain(semi_major_axis, inclination, min_moment_of_inertia)
        self.assertAlmostEqual(actual_gain, expected_gain, places=3)

class TestComputeMinMomentOfInertia(unittest.TestCase):

    def test_diagonal_matrix(self):
        inertia_matrix = np.diag([10, 5, 2])
        expected_min_moment = 2
        actual_min_moment = Controller._compute_min_moment_of_inertia(inertia_matrix)
        self.assertEqual(actual_min_moment, expected_min_moment)

    def test_symmetric_matrix(self):
        inertia_matrix = np.array([[4, 1, 2],
                                   [1, 6, 3],
                                   [2, 3, 5]])
        expected_min_moment = 1.921
        actual_min_moment = Controller._compute_min_moment_of_inertia(inertia_matrix)
        self.assertAlmostEqual(actual_min_moment, expected_min_moment, places=3)

    def test_typical_values_1(self):
        J = np.array([[4.5e-3,-3.2e-4,0],[-3.2e-4,5.1e-3,0],[0,0,3.7e-3]])
        Jmin = Controller._compute_min_moment_of_inertia(J)
        self.assertGreater(Jmin, 0.0)
        self.assertLess(Jmin, 1.0)

class TestBcrossControl(unittest.TestCase):

    def test_zero_angular_velocity(self):
        angular_velocity = np.zeros(3)  # Zero angular velocity
        B_body = np.array([1, 2, 3])
        k = 0.5
        
        expected_dipole_moment = np.zeros(3)  # Expected zero dipole moment
        actual_dipole_moment = Controller._apply_bcross_control(angular_velocity, B_body, k)
        
        np.testing.assert_allclose(actual_dipole_moment, expected_dipole_moment)

if __name__ == "__main__":
    unittest.main()
