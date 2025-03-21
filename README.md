# OpenMM ASE Calculator

This repository provides an implementation of an [OpenMM](https://github.com/openmm/openmm) calculator for the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/) package. It enables seamless integration of OpenMM's Contexts with ASE's tools.

## Installation

This calculator can be installed via pip. Within this repository, simply run:

```bash
pip install -e .
```

## Implemented Properties

The calculator currently supports the following properties:
- **Energies**
- **Forces**

## Quick Start

To use the OpenMM ASE calculator, you'll need an `Atoms` object and an OpenMM `Context`. A utility function, `openmm_topology_to_ase_atoms`, is available to convert an OpenMM `Topology` into an ASE `Atoms` object.

Here's a simple example of setting up the calculator:

```python
from openmm_ase import OpenMMCalculator, openmm_topology_to_ase_atoms

# Convert OpenMM Topology to ASE Atoms
atoms = openmm_topology_to_ase_atoms(pdb.topology, pdb.positions)

# Create and assign the OpenMM calculator
calc = OpenMMCalculator(atoms, context=simulation.context)
atoms.calc = calc
```

## References

1. OpenMM documentation: http://docs.openmm.org/latest/userguide/
1. OpenMM repository: https://github.com/openmm/openmm
1. ASE documentation: https://wiki.fysik.dtu.dk/ase/
1. ASE repository: https://gitlab.com/ase/