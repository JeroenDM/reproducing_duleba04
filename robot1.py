from numpy import cos, sin, array
from scipy.optimize import root
import matplotlib.pyplot as plt

l1 = 1.0
l2 = 1.0

def fk(q, all_links = False):
  
    x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1])
    y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1])
    if all_links:
        x1 = l1 * cos(q[0])
        y1 = l1 * sin(q[0])
        return array([x, y]), array([x1, y1])
    else:
        return array([x, y])

def ik(x, q0):
    solution = root(lambda q : x - fk(q), q0)
    
    return solution

def J(q):
  from numpy import cos, sin, array
  j11 = -l1 * sin(q[0]) - l2 * sin(q[0] + q[1])
  j12 = -l2 * sin(q[0] + q[1])
  j21 =  l2 * cos(q[0] + q[1])
  j22 =  l1 * cos(q[0]) + l2 * cos(q[0] + q[1])
  return array([j11, j12, j21, j22]).reshape(2, 2)
  
def rrplot(q):
  from numpy import cos, sin
  x1 = l1 * cos(q[0])
  y1 = l1 * sin(q[0])
  x2 = l1 * cos(q[0]) + l2 * cos(q[0] + q[1])
  y2 = l1 * sin(q[0]) + l2 * sin(q[0] + q[1])
  x = [ 0.0, x1[0], x2[0] ]
  y = [ 0.0, y1[0], y2[0] ]
  
  plt.figure()
  plt.axis([-3, 3, -3, 3])
  plt.plot(x, y, 'ko-')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.show()