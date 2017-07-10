# Robot manipulator path following

## Description
Python code and IPython notebook files to reproduces results from paper: `Algorithms of trajectory planning with constrained deviation
from a given end-effector path' by Ignacy Dule ̧ba and Iwona Karcz–Dule ̧ba.

More info on the algorithms can be found in the code documentation (documented as doc strings in the functions) or in the paper of course.
In the paper, the algorithms are demostrated on a two-link planar manipulater. The end effector has the follow three different paths within a given accuracy.

## Convergence problems
At the there are convergence problems with Taylor's algorithm for straight line path planning.
They occur when an exact IK solution is added in the middle of a line segment. The robot configuration for this new point can be different than for the other point already calculated. For example, elbow up versus elbow down for a two link manipulator. When interpolating between such two different configurations in joint space, the result deviates to much in cartesian space. Even if the two path points are very close in the cartesian space, they are far appart in the joint space.
