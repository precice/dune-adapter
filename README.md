# DUNE-preCICE adapter

**experimental** preCICE-adapter for DUNE, a modular toolbox for solving partial differential equations.

## Installing the adapter

Build the adapter in the following way:

```bash
cd dune-precice
<path-to-dune-common/bin/dunecontrol> --current all
```

For more detailed information about how `dunecontrol` is used for installing DUNE modules, look at the *Building DUNE Modules* in the [DUNE documentation](https://www.dune-project.org/doc/installation/installation-buildsrc/).

## dune-precice-howto

This module contains a DUNE based solver for the solid participant of the [perpendicular-flap](https://github.com/precice/tutorials/tree/master/perpendicular-flap) tutorial case. The [dune-elastodynamics](https://github.com/maxfirmbach/dune-elastodynamics) module is an external requirement which needs to be separately downloaded.

After the module is built, the executable `dune-perpendicular-flap` case can be found in `dune-precice-howto/build-cmake/examples/`. Copy the executable to the [solid-dune](https://github.com/precice/tutorials/tree/master/perpendicular-flap/solid-dune) folder and then run the case.

## Citing

- Firmbach M. (2021). Aeroelastic simulation of slender wings for electric aircraft - A partitioned approach with DUNE and preCICE, [Master Thesis](https://mediatum.ub.tum.de/node?id=1609293)
- Firmbach M., Callies R. (2021). Aeroelastic simulation of slender wings for electric aircraft - A partitioned approach with DUNE and preCICE, [Conference Contribution](https://athene-forschung.unibw.de/138607)
