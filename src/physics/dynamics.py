"""Compatibility imports; implementations live in orbital_optimizer."""
from orbital_optimizer.dynamics.two_body import EARTH_MU, two_body_dynamics, validate_mu, validate_state
from orbital_optimizer.orbital.mechanics import delta_v, specific_angular_momentum, specific_energy, total_delta_v
