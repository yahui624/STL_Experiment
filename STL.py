from gurobipy import * 
import numpy as np 
import time 
from math import * 

EPS = 1e-2

# Conjunction Node 
class Conjunction(object): 
    def __init__(self, deps=[]) -> None:
        self.deps = deps
        self.constraints = []

# Disjunction Node 
class Disjunction(object): 
    def __init__(self, deps=[]) -> None:
        self.deps = deps
        self.constraints = []

def noIntersection(a, b, c, d):
    # b < c or d < a 
    return Disjunction([c-b, a-d])

def hasIntersection(a, b, c, d): 
    # ~noIntersection 
    return Conjunction([b-c, d-a])

def always(i, a, b, z_phis, piecewise_linear):
    t_i = piecewise_linear[i][1]
    t_i_1 = piecewise_linear[i+1][1]
    conjunctions = [] 
    for j in range(len(piecewise_linear)-1): 
        t_j = piecewise_linear[j][1]
        t_j_1 = piecewise_linear[j+1][1]
        # Recursively check for disjuction 
        conjunctions.append(Disjunction([noIntersection(t_j, t_j_1, t_i + a, t_i_1 + b), z_phis[j]]))
    return Conjunction(conjunctions)

def eventually(i, a, b, z_phis, piecewise_linear): 
    t_i = piecewise_linear[i][1]
    t_i_1 = piecewise_linear[i+1][1]
    z_interval_width = b-a-(t_i_1-t_i)
    disjunctions = []
    for j in range(len(piecewise_linear)-1): 
        t_j = piecewise_linear[j][1]
        t_j_1 = piecewise_linear[j+1][1]
        # Recursively check for conjunction 
        disjunctions.append(Conjunction([hasIntersection(t_j, t_j_1, t_i_1 + a, t_i + b), z_phis[j]]))
    return Conjunction([z_interval_width, Disjunction(disjunctions)])

def bounded_eventually(i, a, b, z_phis, piecewise_linear, tmax): 
    t_i = piecewise_linear[i][1]
    t_i_1 = piecewise_linear[i+1][1]
    z_interval_width = b-a-(t_i_1-t_i)
    disjunctions = []
    for j in range(len(piecewise_linear)-1): 
        t_j = piecewise_linear[j][1]
        t_j_1 = piecewise_linear[j+1][1]
        # Recursively check for conjunction 
        disjunctions.append(Conjunction([hasIntersection(t_j, t_j_1, t_i_1 + a, t_i + b), z_phis[j]]))
    return Disjunction([Conjunction([z_interval_width, Disjunction(disjunctions)]), t_i+b-tmax])

def until(i, a, b, zphi1s, zphi2s, PWL):
    t_i = PWL[i][1]
    t_i_1 = PWL[i+1][1]
    z_intervalWidth = b-a-(t_i_1-t_i)-EPS
    disjunctions = []
    for j in range(len(PWL)-1):
        t_j = PWL[j][1]
        t_j_1 = PWL[j+1][1]
        conjunctions = [hasIntersection(t_j, t_j_1, t_i_1 + a, t_i + b), zphi2s[j]]
        for l in range(j+1):
            t_l = PWL[l][1]
            t_l_1 = PWL[l+1][1]
            conjunctions.append(Disjunction([noIntersection(t_l, t_l_1, t_i, t_i_1 + b), zphi1s[l]]))
        disjunctions.append(Conjunction(conjunctions))
    return Conjunction([z_intervalWidth, Disjunction(disjunctions)])

def release(i, a, b, zphi1s, zphi2s, PWL):
    t_i = PWL[i][1]
    t_i_1 = PWL[i+1][1]
    conjunctions = []
    for j in range(len(PWL)-1):
        t_j = PWL[j][1]
        t_j_1 = PWL[j+1][1]
        disjunctions = [noIntersection(t_j, t_j_1, t_i_1 + a, t_i + b), zphi2s[j]]
        for l in range(j):
            t_l = PWL[l][1]
            t_l_1 = PWL[l+1][1]
            disjunctions.append(Conjunction([hasIntersection(t_l, t_l_1, t_i_1, t_i_1 + b), zphi1s[l]]))
        conjunctions.append(Disjunction(disjunctions))
    return Conjunction(conjunctions)

def mu(i, piecewise_linear, bloat_factor, A, b):
    bloat_factor = np.max([0, bloat_factor])
    b = b.reshape(-1)
    num_edges = len(b)
    conjunctions = []
    for e in range(num_edges): 
        a = A[e, :]
        for j in [i, i+1]: 
            x = piecewise_linear[j][0]
            conjunctions.append(b[e] - np.linalg.norm(a) * bloat_factor - sum([a[k]*x[k] for k in range(len(x))]) - EPS)
    return Conjunction(conjunctions)

def negMu(i, piecewise_linear, bloat_factor, A, b):
    b = b.reshape(-1) 
    num_edges = len(b)
    disjunctions = []
    for e in range(num_edges): 
        a = A[e,:]
        conjunctions = []
        for j in [i, i+1]:
            x = piecewise_linear[j][0]
            conjunctions.append(sum([a[k]*x[k] for k in range(len(x))]) - (b[e] + np.linalg.norm(a) * bloat_factor) - EPS)
        disjunctions.append(Conjunction(conjunctions))
    return Disjunction(disjunctions)


