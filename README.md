---
title: The DUNE adapter
permalink: adapter-dune.html
keywords: DUNE, FSI, C++
summary: "A DUNE module-type adapter to couple to other codes using preCICE"
---

# DUNE-preCICE adapter

**experimental** preCICE-adapter for [DUNE](https://www.dune-project.org/), a modular toolbox for solving partial differential equations.

The DUNE-preCICE adapter is a DUNE module. An example solver that uses the adapter is available in `dune-precice-howto/`. This is another DUNE module, based on the [dune-elastodynamics](https://github.com/maxfirmbach/dune-elastodynamics) module. To build all modules together, the easiest way is to clone all repositories in the same directory, and build all modules present in that directory. See how we install DUNE in the [preCICE Demo VM](https://github.com/precice/vm/blob/develop/provisioning/install-dune.sh).

## Installing the adapter

We need to build the adapter module with `dunecontrol`. For more detailed information about how `dunecontrol` is used for installing DUNE modules, look at the *Building DUNE Modules* in the [DUNE documentation](https://www.dune-project.org/doc/installation/installation-buildsrc/).

Clone the [adapter repository](https://github.com/precice/dune-adapter), navigate into it, and build the module found in the `dune-precice` directory:

```bash
cd dune-precice
<path-to-dune-common/bin/dunecontrol> --current all
```

## dune-precice-howto

This module contains a DUNE based solver for the solid participant of the [perpendicular-flap](https://precice.org/tutorials-perpendicular-flap.html) tutorial case. The [dune-elastodynamics](https://github.com/maxfirmbach/dune-elastodynamics) module is an external requirement, which needs to be separately downloaded.

This module is built in the same way as the adapter. The executable `dune-perpendicular-flap` is in `dune-precice-howto/build-cmake/examples/`. Copy the executable to the [solid-dune](https://github.com/precice/tutorials/tree/master/perpendicular-flap/solid-dune) folder and then run the case.

## Citing

- Firmbach M. (2021). Aeroelastic simulation of slender wings for electric aircraft - A partitioned approach with DUNE and preCICE, [Master Thesis](https://mediatum.ub.tum.de/node?id=1609293)
- Firmbach M., Callies R. (2021). Aeroelastic simulation of slender wings for electric aircraft - A partitioned approach with DUNE and preCICE, [Conference Contribution](https://athene-forschung.unibw.de/138607)
