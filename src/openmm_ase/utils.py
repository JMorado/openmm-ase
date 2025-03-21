"""Module containing utility functions for OpenMM and ASE interoperability."""

def openmm_topology_to_ase_atoms(openmm_topology, positions=None):
    """
    Convert an OpenMM topology to an ASE Atoms object.

    Parameters
    ----------
    openmm_topology : openmm.Topology
        The OpenMM topology object.
    positions : np.ndarray, optional
        The positions of the atoms. If not provided, the positions are set to None.
    """
    from ase import Atoms     
    import ase.units as units
    import numpy as np
    import openmm as mm 

    # Extract information from the OpenMM topology
    symbols = [atom.element.symbol for atom in openmm_topology.atoms()]
    atomic_numbers = [atom.element.atomic_number for atom in openmm_topology.atoms()]
    masses = [
        atom.element.mass.value_in_unit(mm.unit.dalton)
        for atom in openmm_topology.atoms()
    ]
    cell = (
        openmm_topology.getUnitCellDimensions().value_in_unit(mm.unit.nanometer)
        * units.nm
    )
    pbc = cell is not None

    # Create and set the information in the ASE Atoms object
    atoms = Atoms(symbols)
    atoms.set_chemical_symbols(symbols)
    atoms.set_atomic_numbers(atomic_numbers)
    if positions is not None:
        positions = np.array(positions.value_in_unit(mm.unit.nanometer)) * units.nm
        atoms.set_positions(positions)
    atoms.set_masses(masses)
    atoms.set_cell(cell)
    atoms.set_pbc(pbc)

    return atoms
