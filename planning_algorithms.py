import numpy as np
from numpy.linalg import norm
from scipy.optimize import root
from util import *

def exact_path_following(xp, q0, fk):
    """ Solves the inverse kinematics for every path point """
    if xp.shape[1] != 2:
        raise ValueError("Input path must have shape 2 x N")
        
    Np = xp.shape[0]
    qsol = np.zeros((Np, 2))
    qsol[0] = q0
    for i in range(1, Np):
        res = root(lambda q : xp[i] - fk(q), qsol[i-1])
        if res['success']:
            qsol[i] = res['x']
        else:
            print "IK did not converge for path point: " + str(i)
    return qsol

def taylors_algorithm(x0, xN, q0, qN, delta, fk):
    qsol = [q0, qN]
    xsol = [x0, xN]
    
    i = 0 # current point looked at (i to i+1)
    N = 2 # total number of points in solution
    max_iter = 500

    while(i < (N - 1) and max_iter <= 0):
        max_iter -= 1
        # interpolate in joint space
        q_mid, x_mid, e = mid_point2(qsol, xsol, i, fk)

        # check error in task space
        if e <= delta:
            # add point to solution
            add_point(xsol, qsol, x_mid, q_mid, i+1)
            N += 1
            i += 2
        else:
            # refine grid with ik solver
            x_ref, q_ref = refine_grid(xsol, qsol, i, fk)
            add_point(xsol, qsol, x_ref, q_ref, i+1)
            N +=1
    
    if max_iter <= 0:
        print "Not converged, maximum iterations reached"
        
    qsol = np.array(qsol)
    xsol = np.array(xsol)
    return qsol, xsol

def find_close_point(x1, q1, x2s, Nd, J):
    """ Close in joint space """
    dmin = np.inf
    xmin = 0
    Jk = J(q1)
    Jk_inv = pinv(Jk)
    for j in range(Nd):
        xj = x2s[:, j]
        dj = norm(Jk_inv.dot(xj - x1))

        if dj < dmin:
            dmin = dj
            xmin = xj
        
    return xmin, dmin

def local_optimization(xp, q0, delta, fk, J):
    if xp.shape[1] != 2:
        raise ValueError("Input path must have shape 2 x N")       
    Np = xp.shape[0]
    
    # extend path with discretized ball with radius delta (1D in this case) around each path point
    Nd = 10 # discretization points of ball around path point
    xp_ext = []
    
    for i in range(Np):
        x_ext = np.linspace(xp[i][0], xp[i][0], Nd)
        y_ext = np.linspace(xp[i][1] - delta, xp[i][1] + delta, Nd)
        xp_ext.append(np.vstack((x_ext, y_ext)))
        
    # distance in joint space to the different point?
    qsol = [q0]
    xsol = [xp[0]]

    for i in range(Np-1):
        xmin, dmin = find_close_point(xp[i], qsol[i], xp_ext[i+1], Nd, J)

        # find corresponding joint solution
        sol = root(lambda q : xmin - fk(q), qsol[i])
        qmin = sol['x']

        xsol.append(xmin)
        qsol.append(qmin)
    qsol = np.array(qsol)
    xsol = np.array(xsol)
    return qsol, xsol

def trajectory_shortening(xp, qp, delta, fk, zeta=0.001, angle_max=0.5):
    # check xp shape and length
    if xp.shape[1] != 2:
        raise ValueError("Input path must have shape 2 x N")       
    Np = xp.shape[0]
    
    qsol = [qp[0]]
    xsol = [xp[0]]
    
    for i in range(1, Np-1):
        # save current piont
        qi = qp[i]
        xi = xp[i]
        
        # create joint space vectors
        v1 = qp[i-1] - qp[i]
        v2 = qp[i+1] - qp[i]
        nv1 = norm(v1)
        nv2 = norm(v2)
        
        # vector must have minimum length to allow for sufficient shortening
        if (nv1 < zeta) or (nv2 < zeta):
            qsol.append(qi)
            xsol.append(xi)
            continue
        
        # If angle is big, not much shortening can be done
        angle12 = abs(angle(v1, v2, nv1, nv2))
        # print angle12
        if (angle12 > angle_max):
            qsol.append(qi)
            xsol.append(xi)
            continue
            
        # if we get to this part of the loop, we should try to shorten the trajectory
        # from (i-1) to (i+1) by replacing qi
        w = np.linspace(0.1, 4.0, 10) # try different weights
        for wi in w:
            qmc = mass_center(qp[i-1], qp[i], qp[i+1], wi)
            xmc = fk(qmc)
            if dist(xp[i], xp[i+1], xmc) < delta:
                # print "point replaced in iteration " + str(i)
                # this is a good replacement point!
                qi = qmc
                xi = xmc
                break
            
        qsol.append(qi)
        xsol.append(xi)
    
    qsol.append(qp[-1])
    xsol.append(xp[-1])
    qsol = np.array(qsol)
    xsol = np.array(xsol)
    return qsol, xsol
        
    