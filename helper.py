from gurobipy import *
from math import *
from STL import * 
from node import *
import numpy as np
import time

M = 1e3
# a large M causes numerical issues and make the model infeasible to Gurobi
T_MIN_SEP = 1e-2
# see comments in GreaterThanZero
IntFeasTol  = 1e-1 * T_MIN_SEP / M

def setM(v):
    global M, IntFeasTol
    M = v
    IntFeasTol  = 1e-1 * T_MIN_SEP / M

EPS = 1e-2

def _sub(x1, x2):
    return [x1[i] - x2[i] for i in range(len(x1))]

def _add(x1, x2):
    return [x1[i] + x2[i] for i in range(len(x1))]

def L1Norm(model, x):
    xvar = model.addVars(len(x), lb=-GRB.INFINITY)
    abs_x = model.addVars(len(x))
    model.update()
    xvar = [xvar[i] for i in range(len(xvar))]
    abs_x = [abs_x[i] for i in range(len(abs_x))]
    for i in range(len(x)):
        model.addConstr(xvar[i] == x[i])
        model.addConstr(abs_x[i] == abs_(xvar[i]))
    return sum(abs_x)


def add_space_constraints(model, xlist, limits, bloat=0.):
    xlim, ylim = limits
    for x in xlist:
        model.addConstr(x[0] >= (xlim[0] + bloat))
        model.addConstr(x[1] >= (ylim[0] + bloat))
        model.addConstr(x[0] <= (xlim[1] - bloat))
        model.addConstr(x[1] <= (ylim[1] - bloat))

    return None

def add_time_constraints(model, PWL, tmax=None):
    if tmax is not None:
        model.addConstr(PWL[-1][1] <= tmax - T_MIN_SEP)

    for i in range(len(PWL)-1):
        x1, t1 = PWL[i]
        x2, t2 = PWL[i+1]
        model.addConstr(t2 - t1 >= T_MIN_SEP)

def add_velocity_constraints(model, PWL, vmax=3):
    for i in range(len(PWL)-1):
        x1, t1 = PWL[i]
        x2, t2 = PWL[i+1]
        # squared_dist = sum([(x1[j]-x2[j])*(x1[j]-x2[j]) for j in range(len(x1))])
        # model.addConstr(squared_dist <= (vmax**2) * (t2 - t1) * (t2 - t1))
        L1_dist = L1Norm(model, _sub(x1,x2))
        model.addConstr(L1_dist <= vmax * (t2 - t1))

def disjoint_segments(model, seg1, seg2, bloat):
    assert(len(seg1) == 2)
    assert(len(seg2) == 2)
    # assuming that bloat is the error bound in two norm for one agent
    return 0.5 * L1Norm(model, _sub(_add(seg1[0], seg1[1]), _add(seg2[0], seg2[1]))) - 0.5 * (L1Norm(model, _sub(seg1[0], seg1[1])) + L1Norm(model, _sub(seg2[0], seg2[1]))) - 2*bloat*np.sqrt(len(seg1[0])) - EPS

def add_mutual_clearance_constraints(model, PWLs, bloat):
    for i in range(len(PWLs)):
        for j in range(i+1, len(PWLs)):
            PWL1 = PWLs[i]
            PWL2 = PWLs[j]
            for k in range(len(PWL1)-1):
                for l in range(len(PWL2)-1):
                    x11, t11 = PWL1[k]
                    x12, t12 = PWL1[k+1]
                    x21, t21 = PWL2[l]
                    x22, t22 = PWL2[l+1]
                    z_noIntersection = noIntersection(t11, t12, t21, t22)
                    z_disjoint_segments = disjoint_segments(model, [x11, x12], [x21, x22], bloat)
                    z = Disjunction([z_noIntersection, z_disjoint_segments])
                    add_CDTree_Constraints(model, z)


def plan(x0s, specs, bloat, limits=None, num_segs=None, tasks=None, vmax=3., MIPGap=1e-4, max_segs=None, tmax=None, hard_goals=None, size=0.11*4/2):
    if num_segs is None:
        min_segs = 1
        assert max_segs is not None
    else:
        min_segs = num_segs
        max_segs = num_segs

    for num_segs in range(min_segs, max_segs+1):
        for spec in specs:
            clearSpecTree(spec)
        print('----------------------------')
        print('num_segs', num_segs)

        PWLs = []
        m = Model("xref")
        # m.setParam(GRB.Param.OutputFlag, 0)
        m.setParam(GRB.Param.IntFeasTol, IntFeasTol)
        m.setParam(GRB.Param.MIPGap, MIPGap)
        # m.setParam(GRB.Param.NonConvex, 2)
        # m.getEnv().set(GRB_IntParam_OutputFlag, 0)

        for idx_a in range(len(x0s)):
            x0 = x0s[idx_a]
            x0 = np.array(x0).reshape(-1).tolist()
            spec = specs[idx_a]

            dims = len(x0)

            PWL = []
            for i in range(num_segs+1):
                PWL.append([m.addVars(dims, lb=-GRB.INFINITY), m.addVar()])
            PWLs.append(PWL)
            m.update()

            # the initial constriant
            m.addConstrs(PWL[0][0][i] == x0[i] for i in range(dims))
            m.addConstr(PWL[0][1] == 0)

            if hard_goals is not None:
                goal = hard_goals[idx_a]
                m.addConstrs(PWL[-1][0][i] == goal[i] for i in range(dims))

            if limits is not None:
                add_space_constraints(m, [P[0] for P in PWL], limits)

            add_velocity_constraints(m, PWL, vmax=vmax)
            add_time_constraints(m, PWL, tmax)

            handleSpecTree(spec, PWL, bloat, size)
            add_CDTree_Constraints(m, spec.zs[0])

        if tasks is not None:
            for idx_agent in range(len(tasks)):
                for idx_task in range(len(tasks[idx_agent])):
                    handleSpecTree(tasks[idx_agent][idx_task], PWLs[idx_agent], bloat, size)

            conjunctions = []
            for idx_task in range(len(tasks[0])):
                disjunctions = [tasks[idx_agent][idx_task].zs[0] for idx_agent in range(len(tasks))]
                conjunctions.append(Disjunction(disjunctions))
            z = Conjunction(conjunctions)
            add_CDTree_Constraints(m, z)

        add_mutual_clearance_constraints(m, PWLs, bloat)

        # obj = sum([L1Norm(m, _sub(PWL[i][0], PWL[i+1][0])) for PWL in PWLs for i in range(len(PWL)-1)])
        obj = sum([PWL[-1][1] for PWL in PWLs])
        m.setObjective(obj, GRB.MINIMIZE)

        m.write("test.lp")
        print('NumBinVars: %d'%m.getAttr('NumBinVars'))

        # m.computeIIS()
        # import ipdb;ipdb.set_treace()
        try:
            start = time.time()
            m.optimize()
            end = time.time()
            print('sovling it takes %.3f s'%(end - start))
            PWLs_output = []
            for PWL in PWLs:
                PWL_output = []
                for P in PWL:
                    PWL_output.append([[P[0][i].X for i in range(len(P[0]))], P[1].X])
                PWLs_output.append(PWL_output)
            m.dispose()
            return PWLs_output
        except Exception as e:
            m.dispose()
    return [None,]
