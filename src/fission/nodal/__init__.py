"""Nodal diffusion for the whole core (M4).

``mesh``: the model's core as a 3D grid of homogeneous nodes;
``diffusion``: the coarse-mesh finite-difference operator and its eigenvalue
(the reference on a fine mesh, and the global problem of NEM);
``nem``: the nodal expansion method (two-node problems that correct it).
"""
