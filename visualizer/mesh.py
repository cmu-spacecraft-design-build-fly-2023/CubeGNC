import numpy as np

import trimesh
from stl import mesh


mesh = trimesh.load_mesh('visualizer/example_cubesat_v7.STL')


#mesh = trimesh.convex.convex_hull(omesh)
#mesh.export('assets/hulls/chaser_spacecraft_cvx_hull.obj')

scale_factor = 0.005
mesh.apply_scale(scale_factor)



angle = np.pi / 2  # 90 degrees
axis = [1, 0, 0]  
rotation_matrix = trimesh.transformations.rotation_matrix(angle, axis)
mesh.apply_transform(rotation_matrix)

T = trimesh.transformations.translation_matrix([0, 0, 0.35])
mesh.apply_transform(T)



#print(mesh.triangles)
#print("size vertices:", len(mesh.vertices))
#print(mesh.principal_inertia_components)

output_filename = 'scaled_cubesat_' + str(scale_factor) + 'x.STL'
mesh.export(output_filename)

reference_frame = trimesh.creation.axis(axis_length=2)  
combined_mesh = trimesh.util.concatenate([mesh, reference_frame])
combined_mesh.show()

#print(mesh.show())