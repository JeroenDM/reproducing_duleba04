from numpy.linalg import norm, pinv
from numpy import cos, sin, array, pi, dot, cross, arctan2, zeros, linspace, vstack, inf
from scipy.optimize import root
import matplotlib.pyplot as plt

def find_close_point(x1, q1, x2s, Nd, J):
    """ Close in joint space """
    dmin = inf
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

def generate_example_paths():
    """ Generate three vector of 2D points (N x 2) given in the paper """
    Na = 21
    Nb = 17
    Nc = 13
    A0 = [ 1.940, 0.342]
    Af = [-1.860, 0.342]
    B0 = [ 1.700, 0.985]
    Bf = [-1.690, 0.985]
    C0 = [ 1.266, 1.510]
    Cf = [-1.234, 1.510]
    
    x1 = vstack((linspace(A0[0], Af[0], Na),
                 linspace(A0[1], Af[1], Na))).T
    x2 = vstack((linspace(B0[0], Bf[0], Nb),
                 linspace(B0[1], Bf[1], Nb))).T
    x3 = vstack((linspace(C0[0], Cf[0], Nc),
                 linspace(C0[1], Cf[1], Nc))).T
    
    return x1, x2, x3
    
    
def dist(p1, p2, a):
    d = norm(p1 - p2)
    if d <= 1e-12:
        raise ValueError("Line must be given by two different points, d= " + str(d))
    A = abs( (p2[1] - p1[1]) * a[0] - (p2[0] - p1[0]) * a[1] + p2[0] * p1[1] - p2[1] * p1[0] )
    return A / d

# add point in joint and path points array
# WARNING: xsol and qsol are modified by the function
def add_point(xsol, qsol, xin, qin, ind):
    xsol.insert(ind, xin)
    qsol.insert(ind, qin)
    
def refine_grid(xsol, qsol, i, fk):
    x_new = (xsol[i] + xsol[i+1]) / 2
    sol = root(lambda q : x_new - fk(q), qsol[i])
    if sol['success']:
        q_new = sol['x']
    else:
        raise ValueError("IK not converged")
    
    return x_new, q_new

def angle(a, b, na, nb):
    """ At this point it is checked that the norm is not close to zeros"""
    cos_angle = dot(a, b) / (na * nb)
    sin_angle = cross(a, b) / (na * nb)
    
    return arctan2(sin_angle, cos_angle)

def mass_center(q1, q2, q3, w):
    return (q1 + w * q2 + q3) / (2.0 + w)
  
def calc_eta(q_ep, q_met, Np):
    l_ep = zeros(Np-1)
    l_met = zeros(Np-1)
    for k in range(Np-1):
        l_ep[k]  = norm(q_ep[k] - q_ep[k+1])
        l_met[k] = norm(q_met[k] - qs_met[k+1])
    eta = (l_met - l_ep) / l_ep
    return eta
    
def newton(fk, J, x, q0, alpha = 0.1, max_it = 1000, tol = 1e-6):
  it_counter = max_it
  
  qk = q0
  xk = fk(q0)
  Jk = J(q0)
  
  while(it_counter > 0):
    it_counter -= 1
    
    Jk_inv = pinv(Jk)
    q_new = qk + alpha * Jk_inv.dot(x - xk)
    
    xk = fk(q_new)
    if norm(xk - x) <= tol:
      return {'conv': True, 'q': q_new}
    
    qk = q_new
    Jk = J(qk)
    
  return {'conv': False}

def jac_transpose(fk, J, x, q0, alpha = 0.1, max_it = 1000, tol = 1e-6):
  it_counter = max_it
  
  qk = q0
  xk = fk(q0)
  Jk = J(q0)
  
  while(it_counter > 0):
    it_counter -= 1
    
    Jk_t = Jk.T
    q_new = qk + alpha * Jk_t.dot(x - xk)
    
    xk = fk(q_new)
    if norm(xk - x) <= tol:
      return {'conv': True, 'q': q_new}
    
    qk = q_new
    Jk = J(qk)
    
  return {'conv': False}
