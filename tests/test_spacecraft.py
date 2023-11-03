from time import sleep
from models.spacecraft import Spacecraft
from visualizer.mc import Visualizer

config = {
    "mode": "all",
    "initial_attitude":[1.0,0,0,0,0.1,-0.23,0.52],
    "initial_orbit_oe":[1.5e6, 0, 0, 0, 0, 0],
    "mass": 1.7, # kg
    "inertia": [0.003,0.003,0.003,0.0,0.0,0.0],
    "dt": 0.01,
    "flexible": False,
    "gravity_order": 5,
    "gravity_degree": 5,
    "drag": True
}


spacecraft = Spacecraft(config)
print(spacecraft.get_state())




stl_filename = "visualizer/scaled_cubesat_0005x.STL"
viz = Visualizer(stl_filename)
viz.start_visualization()
viz.set_grid(False)



dt = config["dt"]
N = 2000

for i in range(N):
    spacecraft.advance()
    q = spacecraft.get_state()[6:10]
    print(q)
    viz.apply_rotation(q)
    sleep(dt)


