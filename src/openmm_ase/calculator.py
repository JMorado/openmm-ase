"""This module defines an ASE interface to OpenMM."""

import ase.units as units
import numpy as np
from ase.calculators.calculator import Calculator

import openmm as mm


class OpenMMCalculator(Calculator):
    """
    Interface to OpenMM using an ASE calculator.

    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object to calculate.
    label : str
        The label of the calculator.
    context : openmm.Context
        The OpenMM context to use for the calculation.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, atoms=None, label=None, context=None):
        # if not have_openmm:
        #     raise RuntimeError("OpenMM is not installed.")
        Calculator.__init__(self, label, atoms)
        self._context = context

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        if system_changes:
            if "energy" in self.results:
                del self.results["energy"]
            if "forces" in self.results:
                del self.results["forces"]

        if "energy" not in self.results:
            # Set the positions of the atoms in the OpenMM context
            self._context.setPositions(atoms.get_positions() * mm.unit.angstrom)

            # Compute the potential energy
            energy = self._context.getState(getEnergy=True).getPotentialEnergy()

            # Store the energy in the results dictionary
            self.results["energy"] = energy * units.kJ / units.mol
        if "forces" not in self.results:
            # Set the positions of the atoms in the OpenMM context
            self._context.setPositions(atoms.get_positions() * mm.unit.angstrom)

            # Compute the forces
            forces = self._context.getState(getForces=True, getEnergy=True).getForces(
                asNumpy=True
            )
            energy = self._context.getState(getEnergy=True).getPotentialEnergy()

            # Store the forces and energy in the results dictionary
            self.results["forces"] = forces * units.kJ / units.mol / units.nm
            self.results["energy"] = energy * units.kJ / units.mol
