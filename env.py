import numpy as np

wall_half_width = 0.05
A = np.array([[-1, 0], [1, 0], [0, -1], [0, 1]])

# Walls 
walls = []

walls.append(np.array([0, 0, 0, 11], dtype = np.float64))
walls.append(np.array([12, 12, 0, 11], dtype = np.float64))
walls.append(np.array([0, 12, 0, 0], dtype = np.float64))
walls.append(np.array([0, 12, 11, 11], dtype = np.float64))

walls.append(np.array([3, 3, 3, 4], dtype = np.float64))
walls.append(np.array([3, 3, 7, 8], dtype = np.float64))
walls.append(np.array([9, 9, 3, 4], dtype = np.float64))
walls.append(np.array([9, 9, 7, 8], dtype = np.float64))

walls.append(np.array([6, 6, 2, 9], dtype = np.float64))
walls.append(np.array([5, 7, 2, 2], dtype = np.float64))
walls.append(np.array([5, 7, 9, 9], dtype = np.float64))


# Obseveration Spots 
b1 = np.array([-0.8, 2.2, -0.8, 2.2], dtype = np.float64)
b2 = np.array([-0.8, 2.2, -8.8, 10.2], dtype = np.float64)
b3 = np.array([-9.8, 11.2, -0.8, 2.2], dtype = np.float64)
b4 = np.array([-9.8, 11.2, -8.8, 10.2], dtype = np.float64)
observation_spots = [(A, b1), (A, b2), (A, b3), (A, b4)]

# Transmitters 
b5 = np.array([-0.8, 2.2, -4.8, 6.2], dtype = np.float64)
b6 = np.array([-9.8, 11.2, -4.8, 6.2], dtype = np.float64)
transmitters = [(A, b5), (A, b6)]

# Charging Statsions 
b7 = np.array([-4, 6, -4, 7], dtype = np.float64)
b8 = np.array([-6, 8, -4, 7], dtype = np.float64)
charging_stations = [(A, b7), (A,b8)]


# obs 
obs = []
for wall in walls:
    if wall[0]==wall[1]:
        wall[0] -= wall_half_width
        wall[1] += wall_half_width
    elif wall[2]==wall[3]:
        wall[2] -= wall_half_width
        wall[3] += wall_half_width
    else:
        raise ValueError('wrong shape for axis-aligned wall')
    wall *= np.array([-1,1,-1,1])
    obs.append((A, wall))