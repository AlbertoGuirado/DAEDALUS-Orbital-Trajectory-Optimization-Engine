"""Orbital mechanics and simulation in SI units (Phase 1)."""

from .dynamics.two_body import EARTH_MU, two_body_dynamics
from .orbital.elements import OrbitalElements, elements_to_state, state_to_elements
from .orbital.mechanics import (
    STANDARD_GRAVITY, delta_v, eccentricity_vector, node_vector, orbital_period,
    propellant_mass, specific_angular_momentum, specific_energy, total_delta_v,
)
from .orbital.transfers import HohmannTransfer, hohmann_transfer
from .simulation.impulsive import propagate_with_impulses
from .simulation.propagator import PropagationResult, propagate

__all__ = [
    "EARTH_MU", "two_body_dynamics", "OrbitalElements", "elements_to_state",
    "state_to_elements", "delta_v", "eccentricity_vector", "node_vector",
    "orbital_period", "specific_angular_momentum", "specific_energy",
    "total_delta_v", "PropagationResult", "propagate", "STANDARD_GRAVITY",
    "propellant_mass", "HohmannTransfer", "hohmann_transfer",
    "propagate_with_impulses",
]
