import os
from time import sleep
import numpy as np

import meshcat
import meshcat.geometry as g
import meshcat.transformations as tf


class Visualizer:

    def __init__(self, mesh_stl_path  = "visualizer/scaled_cubesat_0005x.STL"):
        self.viz = meshcat.Visualizer()
        #stl_filename = "visualizer/scaled_cubesat_0005x.STL"

        self.mesh = meshcat.geometry.StlMeshGeometry.from_file(mesh_stl_path)
        self.viz["spacecraft"].set_object(self.mesh)


    def start_visualization(self):
        self.viz.open()
        self.viz["/Background"].set_property("top_color", [0, 0, 0.1])
        self.viz["/Background"].set_property("bottom_color", [0, 0, 0.1])

    def set_grid(self, on):
        self.viz["/Background"].set_property("grid", on)

    def get_link(self):
        pass


    def apply_position(self, translation_vector):
        # [x, y, z] position
        self.viz["spacecraft"].set_transform(tf.translation_matrix(translation_vector))

    def apply_attitude(self, quat):
        self.viz["spacecraft"].set_transform(tf.quaternion_matrix(quat))



# Temporary local testing
if __name__ == "__main__":


    stl_filename = "visualizer/scaled_cubesat_0005x.STL"
    viz = Visualizer(stl_filename)
    viz.start_visualization()
    viz.set_grid(False)
    sleep(3)
    viz.apply_position([1,0,0])
    q = np.random.rand(4)
    q = q / np.linalg.norm(q)
    q = [0.7071068, 0.7071068, 0, 0]
    viz.apply_attitude(q)
    sleep(3)
    viz.apply_attitude(q)
    sleep(3)


