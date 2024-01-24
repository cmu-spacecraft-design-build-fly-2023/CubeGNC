import numpy as np

GM = 3.986004415e14
class Controller:

    def __init__(self,
                 spacecraft,
                 actuator,
                 k = None):
        """
        Param initialization for B-Cross
        """
        semi_major_axis = spacecraft.get_osculating_state()[0]
        inclination = spacecraft.get_osculating_state()[2]
        min_inertia = self._compute_min_moment_of_inertia(spacecraft.J)
        if k is None:
            self.k = self._compute_bcross_gain(
                semi_major_axis, inclination, min_inertia)
        else:
            self.k = k

        self.maximum_dipole_moment = np.array(actuator.maximum_dipole_moment)
        self.output_range = np.array(actuator.output_range)

    @staticmethod
    def _compute_min_moment_of_inertia(inertia_matrix):
        """
        Computes the the minimum eigenvalue of the inertia matrix.

        Args:
            inertia_matrix: The Inertia Matrix for a rigid-body

        Returns:
            The calculated minimum moment of inertia
        """
        return np.min(np.linalg.eigvals(inertia_matrix))

    @staticmethod
    def _compute_bcross_gain(semi_major_axis, inclination, min_moment_of_inertia):
        """
        Calculates the bcross gain value as per the formula given in Avanzini
        and Giulietti's paper for choosing a reasonable gain value[Eqn 30].

        Args:
            semi_major_axis: Semi-major axis of the orbit in meters(m).
            inclination: Inclination of the orbit in radians(rad).
            min_moment_of_inertia: Minimum moment of inertia of the spacecraft in kg*m^2.

        Returns:
            The calculated gain value.
        """

        orbital_rate = GM / np.sqrt(semi_major_axis**3)

        k_gain = 2 * orbital_rate * (1 + np.sin(inclination)) * min_moment_of_inertia
        return k_gain

    @staticmethod
    def _apply_bcross_control(angular_velocity_body, mag_field_body, k):
        """
        Calculates the commanded magnetic dipole moment value based on the equations
        mentioned in Fundamentals of Spacecraft Attitude Determination and Control
        by F. Landis Markley , John L. Crassidis, Page 308

        Args:
            angular_velocity_body: angular velocity in body frame with units rad/s
            mag_field_body: magnetic field in body frame in Tesla (T)
            k: positive scalar gain

        Returns:
            The commanded magnetic dipole moment value
        """
        
        B_norm = np.linalg.norm(mag_field_body)
        b = mag_field_body / B_norm
        commanded_mag_dipole_moment = (k / B_norm) * np.cross(angular_velocity_body, b)
        return commanded_mag_dipole_moment

    @staticmethod
    def _saturate_mag_dipole(commanded_mag_dipole_moment, maximum_dipole_moment):
        """
        Saturate the magnetic dipole to be in the range [-self.maximum_dipole_moment, self.maximum_dipole_moment]
        
        Args:
            commanded_mag_dipole_moment: calculated magnetic dipole moment using the bcross_control_function Am^2

        Returns:
            saturated_mag_dipole_moment: saturated magentic dipole moment value
        """
        saturated_mag_dipole_moment = np.clip(commanded_mag_dipole_moment, -maximum_dipole_moment, maximum_dipole_moment)
        return saturated_mag_dipole_moment

    @staticmethod
    def _scale_mag_dipole(saturated_mag_dipole_moment, maximum_dipole_moment, output_range):
        """
        Scales the saturated magnetic dipole moment value to fit the output range.
        Args:
            commanded_mag_dipole_moment: calculated magnetic dipole moment using the bcross_control_function Am^2

        Returns:
            scaled_mag_dipole_moment: scaled magentic dipole moment value that fits the output range
        """
        scaled_mag_dipole_moment = (saturated_mag_dipole_moment / maximum_dipole_moment) * output_range
        return scaled_mag_dipole_moment

    def get_control_dipole(self, angular_velocity_body, magnetic_field_body):
        """
        Computes the final magnetic dipole moment to be used for control with saturation and scaling applied
        used to detumble the satellite

        Args:
            angular_velocity_body: angular velocity in body frame with units rad/s
            magnetic_field_body: magnetic field in body frame in Tesla (T)
        
        Returns:
            final_mag_dipole: final magnetic dipole moment to be used for control with saturation and scaling applied
        """
        b_cross_control_dipole = self._apply_bcross_control(
            angular_velocity_body, magnetic_field_body, self.k)
        saturated_mag_dipole = self._saturate_mag_dipole(
            b_cross_control_dipole, self.maximum_dipole_moment)
        scaled_mag_dipole = self._scale_mag_dipole(
            saturated_mag_dipole, self.maximum_dipole_moment, self.output_range)
        return scaled_mag_dipole
