import ase.optimize
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from ase.constraints import FixInternals
from openmm_ase import OpenMMCalculator, openmm_topology_to_ase_atoms
from openmmml import MLPotential

# Load PDB file and set the FFs
pdb = app.PDBFile("alanine-dipeptide-explicit.pdb")
ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# Create the OpenMM MM System and ML potential
mmSystem = ff.createSystem(pdb.topology, nonbondedMethod=app.PME)
potential = MLPotential("ani2x")

# Choose the ML atoms
mlAtoms = [a.index for a in next(pdb.topology.chains()).atoms()]

# Create the mixed ML/MM system (we're using the nnpops implementation for performance)
mixedSystem = potential.createMixedSystem(
    pdb.topology, mmSystem, mlAtoms, interpolate=False, implementation="nnpops"
)

# Choose to run on a GPU (CUDA)
platform = mm.Platform.getPlatformByName("CUDA")

# Setup the simulation
simulation = app.Simulation(
    pdb.topology, mixedSystem, mm.VerletIntegrator(1.0 * unit.femtosecond), platform
)

# Create the ASE calculator
atoms = openmm_topology_to_ase_atoms(pdb.topology, pdb.positions)

# Set the calculator
calc = OpenMMCalculator(atoms, context=simulation.context)
atoms.calc = calc

# ------------------------------------------------------------------ #
#                           2D Relaxed PES                           #
# ------------------------------------------------------------------ #

dihedral1 = [4, 6, 8, 14]
dihedral2 = [8, 14, 16, 18]

# Set the dihedral angles
grid = np.linspace(-180, 180, 10)

# Set the optimization parameters
basename = "alanine-dipeptide"
ase_options = {}
fmax = 0.01
steps = 1000

for dih1_val in grid:
    for dih2_val in grid:
        # Remove any previous constraints
        atoms._del_constraints()

        # Set the dihedral angles
        atoms.set_dihedral(*dihedral1, dih1_val)
        atoms.set_dihedral(*dihedral2, dih2_val)
        x = atoms.get_positions() * ase.units.Angstrom

        # Create the constraints
        c = FixInternals(dihedrals_deg=[[dih1_val, dihedral1], [dih2_val, dihedral2]])
        atoms.set_constraint(c)

        # Run the optimization
        opt = ase.optimize.FIRE(
            atoms=atoms,
            trajectory=basename + f"_opt_{dih1_val:.2f}_{dih2_val:.2f}.traj",
            restart=basename + "_opt.pkl",
            **ase_options,
        )
        opt.run(fmax=fmax, steps=steps)

        # Get the energy
        energy = atoms.get_potential_energy()

        # Write the PDB file
        app.PDBFile.writeFile(
            pdb.topology,
            atoms.get_positions() * unit.angstrom,
            open(basename + f"_opt_{dih1_val:.2f}_{dih2_val:.2f}.pdb", "w"),
        )
        print(
            "Dihedral1:",
            atoms.get_dihedral(*dihedral1),
            "Dihedral2:",
            atoms.get_dihedral(*dihedral2),
            "Energy:",
            energy,
        )
