import numpy as np
class Magnetorquer():
    def __init__(self, voltage, trace_width, trace_thickness, no_of_turns):
        self.voltage = voltage
        self.trace_width = trace_width
        self.trace_thickness = trace_thickness
        self.no_of_turns = no_of_turns
        self.trace_area = trace_thickness * trace_width
        self.maximum_dipole_moment = self.compute_max_dipole_strength()
        self.output_range = self.maximum_dipole_moment

    def compute_coil_resistance(self, length, width):
        copper_resistivity = 1.77e-8
        perimeter = 2 * (length + width)
        return copper_resistivity * self.no_of_turns * perimeter / self.trace_area

    def compute_max_dipole_strength(self):
        coil_resistance = self.compute_coil_resistance(0.1, 0.1)
        print(coil_resistance)
        coil_current = self.voltage / coil_resistance
        maximum_dipole_moment = self.no_of_turns * coil_current * 0.1 * 0.1
        maximum_dipole_moment = maximum_dipole_moment*np.ones(3)
        return maximum_dipole_moment

# Temporary local testing
if __name__ == "__main__":
    torquer = Magnetorquer(10,  0.6e-3, 0.036e-3, 315)
    maximum_dipole_moment = torquer.compute_max_dipole_strength()
    print(maximum_dipole_moment)
