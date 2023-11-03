import os
from time import sleep
import numpy as np

import meshcat
import meshcat.geometry as g
import meshcat.transformations as tf




class Visualizer:

    def __init__(self, mesh_stl_path):
        self.viz = meshcat.Visualizer()

        #stl_filename = "visualizer/example_cubesat_v7.STL"
        #stl_filename = "visualizer/scaled_cubesat_0005x.STL"

        self.mesh = meshcat.geometry.StlMeshGeometry.from_file(mesh_stl_path)
        self.viz["spacecraft"].set_object(self.mesh)


    def start_visualization(self):
        self.viz.open()


    def get_link(self):
        pass


    def apply_translation(self, translation_vector):
        # [x, y, z] translation
        pass

    def apply_rotation(self, quat):
        pass


viz = meshcat.Visualizer()
viz.open()


stl_filename = "visualizer/scaled_cubesat_0005x.STL"
mesh = meshcat.geometry.StlMeshGeometry.from_file(stl_filename)


translation = [0.0, 0.0, 0.0]  # [x, y, z] translation
rotation = tf.rotation_matrix(np.pi, [0, 0, 0])
pose = np.dot(tf.translation_matrix(translation), rotation)



viz["/Background"].set_property("grid", False)
viz["/Background"].set_property("top_color", [0, 0, 0.1])
viz["/Background"].set_property("bottom_color", [0, 0, 0.1])

meshcat_object = viz["cubesat"]
meshcat_object.set_object(mesh)


sleep(5)


