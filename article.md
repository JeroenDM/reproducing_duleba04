# Algorithms for robot manipulator path following

## Introduction

The article `Algorithms of trajectory planning with constrained deviation from a given end-effector path` four different trajectory planning algorithms are described. The task consist of following given path points with the end-effector. Some deviation of the given path points is allowed.

## Algorithms description

In the first algoritm, the inverse kinematic equations are solved for each path point. The solution of the previous path point is used as initial guess when solving the inverse kinematic equations for the current path point. The solution of this problem is used as a refenrence, or even input, for the other algorithms.

## Simulation results for a planar 2-link manipulator

## Simulation results for a planar 3-link manipulator?

## Conclusion

## Convergence problem with Taylor's straight line planning

[pr1](./pr1.png)

Preamble: imagine a two link planar manipulator having an albow up

Conclusion: only use this for short line segments. This reduces the risk of getting stuk between two configurations.
This means that for the two boundary points, the solution of the ik should converge to the same configuration for (the same) initial conditions.
For analytical ik, the same cocnfiguration should be selected.
Should work close to the boundary of two configurations.
