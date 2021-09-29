# dune-elastodynamics
!!! PROTOTYPE !!!

Created during my master's thesis for testing partitioned coupling schemes with
the library preCICE.

This DUNE module is used to solve the linear equations of elastodynamics.
Assemblers for the calculation of stiffness, mass and lumped mass operators are provided.
Different explicit and implicit time-stepping methods for the solution of the corresponding
system of second order ordinary differential equations are implemented. 

## Installation

It is assumed the dune core modules and all other necessary libraries are allready installed.

Dependencies:

* dune-common, dune-geometry, dune-uggrid, dune-grid, dune-localfunctions, dune-istl, dune-functions, dune-foamgrid, dune-typetree
* MPI, Message Passing Interface library
* OpenMP

Run the following command to build the module:
`dunecontrol --current all`
