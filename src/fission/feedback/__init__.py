"""Thermal-hydraulic feedback for the core (M4b).

``water``: IAPWS-IF97 region 1 (density, enthalpy of compressed water);
``fuel``: steady temperatures of a fuel rod (pellet, gap, clad);
``core``: nodal diffusion iterated with channel enthalpy rise and rod
temperatures until the power and the node states agree.
"""
