"""fes_ml base package."""

__version__ = "0.0.1"
__author__ = "Joao Morado"

from .calculator import OpenMMCalculator
from .utils import openmm_topology_to_ase_atoms