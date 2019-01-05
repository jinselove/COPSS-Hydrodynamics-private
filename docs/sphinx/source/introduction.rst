.. _introductions:

Introductions
===============

**COPSS**
---------

**COPSS** (Continuum Particle Simulation Suite) is an open source,
[**LIBMESH**](http://libmesh.github.io/) based, software for continuum simulations.
The package is designed to be easy to use, extendable and scalable. It currently includes
two modules, [**COPSS-Hydrodynamics**](https://bitbucket.org/COPSS/copss-hydrodynamics-public.git)
to solve hydrodynamic interactions in colloidal suspensions and [**COPSS-Polarization**]
(https://bitbucket.org/COPSS/copss-polarization-public) to solve electrostatic interactions
between dielectric particles. The algorithms beneath **COPSS** have been published or
under-review, but the code framework, user-interface, etc., are still rough. We are
working on improving **COPSS** and appreciate your contributions.

**COPSS-Hydrodynamics**
----------------------------

**COPSS-Hydrodynamics** solves the hydrodynamic interactions in colloidal suspensions by
directly solve the Stokes flow.  It is based on an efficient $O(N)$ computational approach
to model the dynamics of hydrodyna- -mically interacting Brownian or micron-sized particles
in arbitrary geometries. A parallel finite element Stokes' solver is the center of the
algorithm.
