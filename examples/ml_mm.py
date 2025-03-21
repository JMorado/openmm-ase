import openmm as mm
import openmm.app as app
import openmm.unit as unit
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

# Get the energy
EV_TO_KCAL_MOL = 23.0605419453293
energy = atoms.get_potential_energy() * EV_TO_KCAL_MOL
print("Energy:", energy)

# Get the OpenMm Energy
simulation.context.setPositions(pdb.positions)
state = simulation.context.getState(getEnergy=True)
print(
    "OpenMM Energy:",
    state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole),
)
