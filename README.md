# dune-adapter
**experimental** preCICE-adapter for DUNE, a modular toolbox for solving partial differential equations

## dune-precice
A adapter module enabeling multiphysics simulation inside DUNE with the help of preCICE [1][2].

The adapter can be built with:
`<path-to-dune-common/bin/dunecontrol> --current all`

For more detailed installation instructions have a look at:
https://www.dune-project.org/doc/installation/

## dune-precice-howto
This module contains tutorial cases to get familiar with preCICE and the respective dune-adapter.

The dune-adapter module depends on the DUNE core modules (version >= 2.7), dune-precice (precice version >= 2.0.0) and
on dune-elastodynamics, dune-foamgrid (version >= 2.7) providing the necessary tools for structural simulation.

The dune-elastodynamics module can be found here:
https://github.com/maxfirmbach/dune-elastodynamics

The executable `dune-perpendicular-flap` for the perpendicular-flap tutorial case can be
found in:
`dune-precice-howto/build-cmake/examples/`
after the module was built and needs to be copied to the `solid-dune` folder inside the respective tutorial case.

<a id="1">[1]</a> 
Firmbach M. (2021).
[Aeroelastic simulation of slender wings for electric aircraft - A partitioned approach with DUNE and preCICE](https://mediatum.ub.tum.de/node?id=1609293), Master's Thesis

<a id="2">[2]</a> 
Firmbach M., Callies R. (2021).
[Aeroelastic simulation of slender wings for electric aircraft - A partitioned approach with DUNE and preCICE](https://athene-forschung.unibw.de/138607), Conference Contribution
