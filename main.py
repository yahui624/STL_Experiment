from viz import * 
from STL import * 
from helper import * 
from node import *
from env import *


def main(): 
    n_rovers = 2
    
    x0s = [[5.5, 6.5], [5.5, 4.5], [7.5, 6.5], [7.5, 4.5]]
    x0s = x0s[:n_rovers]

    tc = 10.
    td = 10.
    tmax = 8.
    vmax = 5.

    finally_visit_observs = []
    specs = []
    for i in range(n_rovers):
        charging = Node('or', deps=[Node('mu', info={'A':A, 'b':b}) for A,b in charging_stations])
        charge_in_tc = Node('F', deps=[charging, ], info={'int':[0, tc]})
        phi_1 = Node('A', deps=[Node('or', deps=[charging, charge_in_tc]),], info={'int':[0,tmax]})

        # Ois = [Node('mu', info={'A':A, 'b':b}) for A, b in observation_spots]
        # phi_2 = Node('and', deps=[Node('F', deps=[Oi,], info={'int':[0,tmax]}) for Oi in Ois])

        transmitting = Node('or', deps=[Node('mu', info={'A':A, 'b':b}) for A, b in transmitters])
        transmitting_in_td = Node('F', deps=[transmitting, ], info={'int':[0, td]})
        notOis = [Node('negmu', info={'A':A, 'b':b}) for A, b in observation_spots]
        notO = Node('and', deps=notOis)
        phi_3 = Node('A', deps=[Node('or', deps=[notO, transmitting_in_td]),], info={'int':[0,tmax]})

        avoid_obs = Node('and', deps=[Node('negmu', info={'A':A, 'b':b}) for A, b in obs])
        phi_4 = Node('A', deps=[avoid_obs,], info={'int':[0,tmax]})

        specs.append(Node('and', deps=[phi_1, phi_3, phi_4]))

        Ois = [Node('mu', info={'A':A, 'b':b}) for A, b in observation_spots]
        finally_visit_observs.append([Node('F', deps=[Oi,], info={'int':[0,tmax]}) for Oi in Ois])

    PWL = plan(x0s, specs, bloat=0.21, tasks=finally_visit_observs, MIPGap = 0.7, num_segs=10, tmax=tmax, vmax=vmax)

    plots = [[transmitters, 'y'], [charging_stations, 'b'], [observation_spots, 'g'], [obs, 'k']]
    return x0s, plots, PWL

if __name__ == '__main__': 
    results = viz(main)