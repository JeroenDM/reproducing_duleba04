import numpy as np
from numpy.linalg import norm
from scipy.optimize import root
from util import *

def exact_path_following(xp, q0, fk):
    """ Solves the inverse kinematics for every path point.
    
    Uses root solving algorithm from scipy.optimize to solve the inverse kinematics problem.
    Achieved by solving the equation fk(q) - xp = 0 for ever path point.
    Uses the solution of the previous point as initial value for the next point.
    
    Parameters
    ----------
    xp : ndarray (n, 2)
        Cartesian coordinates for the end effector path, consisting of n points.
    q0 : ndarray (2,)
        Initial guess for joint values when solving the ik problem for the first point
    fk : fun
        Function that takes as input a list or array with joint values and
        returns the cartesian coordinates of the end effector.
        
    Returns
    -------
    dict
        A dictionnary with keys
        success : True when the algorithm converged, False otherwise
        q : (n, 2) ndarray with joint solutions
        
    """
    if xp.shape[1] != 2:
        raise ValueError("Input path must have shape 2 x N")
        
    Np = xp.shape[0]
    qsol = np.zeros((Np, 2))
    qsol[0] = q0
    suc = True
    
    for i in range(1, Np):
        res = root(lambda q : xp[i] - fk(q), qsol[i-1])
        if res['success']:
            qsol[i] = res['x']
        else:
            print "IK did not converge for path point: " + str(i)
            suc = False
            
    #return {'success': suc, 'q': qsol}
    return qsol, xp, suc

def taylors_algorithm(x0, xN, q0, qN, delta, fk):
    """ Linear interpolation in joint space.
    
    This algorithm generates a path close to a straight line segment, given by points x0 and xN.
    It uses linear interpolation in joint space. If this results in a cartesian end effector position
    close enough (delta) to the original line, it is added to the path.
    If not, an extra exact ik solution is added to the path.
    
    Parameters
    ----------
    x0 : ndarray (2,)
        Cartesian starting point on the line segment.
    xN : ndarray (2,)
        Cartesian final point on the line segmennt.
    q0 : ndarray (2,)
        IK solution for x0
    qN : ndarray (2,)
        IK solution for xN
    delta : float
        Maximum end effector deviation from the line segment.
    fk : fun
        Function that takes as input a list or array with joint values and
        returns the cartesian coordinates of the end effector.
        
    Returns
    -------
    dict
        A dictionnary with keys
        success : True when the algorithm converged, False otherwise
        q : (n, 2) ndarray with joint solutions
        x : (n, 2) ndarray with cartesian path of the end effector
    
    """
    qs = [q0, qN]
    xs = [x0, xN]
    
    i = 0 # current point looked at (i to i+1)
    N = 1 # index of the last point in solution vector
    max_iter = 500
    suc = True
    

    while(i < N and max_iter > 0):
        max_iter -= 1
        # interpolate in joint space
        q_mid = (qs[i] + qs[i+1]) / 2
        x_mid = fk(q_mid)
        
        x_goal = (xs[i] + xs[i+1]) / 2

        # check error in task space
        if norm(x_goal - x_mid) <= delta:
            # add point to solution
            add_point(xs, qs, x_mid, q_mid, i+1)
            N += 1
            i += 2
        else:
            # refine grid with ik solver
            x_ref, q_ref = refine_grid(xs, qs, i, fk)
            add_point(xs, qs, x_ref, q_ref, i+1)
            N +=1
    
    if max_iter <= 0:
        print "Maximum iterations reached in taylor's algorithm"
        suc = False
        
    qsol = np.array(qs)
    xsol = np.array(xs)
    
    #return {'success': suc, 'q': qsol, 'x': xsol}
    return qsol, xsol, suc

def local_optimization(xp, q0, delta, fk, J):
    """ Locally optimize path length in joint space
    
    Sequentially plan the motion for a path xp. From every point to the next,
    the path length in joint space is minimized.
    The end effector has a maximum deviation from the path in cartesian space, given by delta.
    
    Parameters
    ----------
    xp : ndarray (n, 2)
        Cartesian coordinates for the end effector path of lenght n.
    q0 : ndarray (2,)
        Initial guess for joint values when solving the ik problem for the first point
    fk : fun
        Takes as input a list or array with joint values and
        returns the cartesian coordinates of the end effector.
    J : fun
        Takes as input a list or arrray with joint values and
        returns the jacobian for the robot at the given points
        (as a numpy.ndarray).
        
    Returns
    -------
    dict
        A dictionnary with keys
        success : True when the algorithm converged, False otherwise
        q : (n, 2) ndarray with joint solutions
        
    """
    
    if xp.shape[1] != 2:
        raise ValueError("Input path must have shape 2 x N")
    
    Np = xp.shape[0]
    suc = True
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
        if sol['success']:
            qmin = sol['x']
        else:
            suc = False
            print "IK did not converge for path point: " + str(i)
            break
        xsol.append(xmin)
        qsol.append(qmin)
        
    qsol = np.array(qsol)
    xsol = np.array(xsol)
    
    #return {'success': suc, 'q': qsol, 'x': xsol}
    return qsol, xsol, suc

def trajectory_shortening(xp, qp, delta, fk, zeta=0.001, angle_max=0.5):
    """ Linear interpolation in joint space
    
    This algorithm tries to generate a path close to a straight line segment, given by points x0 and xN. It used a linear interpolation in joint space when this results in a cartesian end effector position close enough (delta) to the original line. If not, an extra exact ik solution is added to the path.
    
    Parameters
    ----------
    xp : ndarray (n,2)
        Cartesian coordinates for the end effector path of lenght n.
    qp : ndarray (n,2)
        Joint solution path which we will try to make shorter.
    delta : float
        Maximum deviation from cartesian path in joint space.
    fk : fun
        Function that takes as input a list or array with joint values en returns the cartesian coordinates of the end effector.
    zeta : float
        Minimum lenght of a joint space vector to consider trajectory shortening for the current point. (Default 0.001)
    angle_max : float
        Maximum angle between to consecutive joint space vectors to consider trajectory shortening. (Default 0.5)
        
    Returns
    -------
    dict
        A dictionnary with keys
        success : True when the algorithm converged, False otherwise
        q : (n, 2) ndarray with joint solutions
        x : (n, 2) ndarray with cartesian path of the end effector
    
    """
    # check xp shape and length
    if xp.shape[1] != 2:
        raise ValueError("Input path must have shape 2 x N")       
    Np = xp.shape[0]
    suc = True
    
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
        if (angle12 < angle_max):
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
    
    #return {'success': suc, 'q': qsol, 'x': xsol}
    return qsol, xsol, suc
        
    