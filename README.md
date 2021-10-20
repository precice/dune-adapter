# dune-adapter
**experimental** preCICE-adapter for DUNE, a modular toolbox for solving partial differential equations

A adapter module enabeling multiphysics simulation inside DUNE with the help of preCICE [1][2].

The dune-adapter module depends on the DUNE core modules and
on dune-elastodynamics, dune-foamgrid providing the necessary tools for structural simulation.

The dune-elastodynamics module can be found here:
https://github.com/maxfirmbach/dune-elastodynamics

The adapter can be build with:
`<path-to-dune-common/bin/dunecontrol> --current all`

For more detailed installation instructions have a look at:
https://www.dune-project.org/doc/installation/

## tutorial

The executable for the perpendicular-flap tutorial case can be
found in:
`dune-precice/build-cmake/src/`
after the module was build and needs to be copied to the solid-dune folder.

<a id="1">[1]</a> 
Firmbach M. (2021).
[Aeroelastic simulation of slender wings for electric aircraft - A partitioned approach with DUNE and preCICE](https://mediatum.ub.tum.de/node?id=1609293), Master's Thesis

<a id="2">[2]</a> 
Firmbach M., Callies R. (2021).
[Aeroelastic simulation of slender wings for electric aircraft - A partitioned approach with DUNE and preCICE](https://athene-forschung.unibw.de/138607), Conference Contribution
