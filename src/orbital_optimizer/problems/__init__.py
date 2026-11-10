"""Adapters that express orbital problems through the generic optimization interfaces.

This layer may import both the physics (dynamics, orbital, simulation) and
orbital_optimizer.optimization; neither of them imports it.
"""

from .coplanar_transfer import (
    TANGENCY_FACTOR, CoplanarTwoImpulseProblem, PlanVerification, TwoImpulsePlan,
    hohmann_two_impulse_plan, reach_margin, verify_two_impulse_plan,
)

__all__ = [
    "CoplanarTwoImpulseProblem", "TwoImpulsePlan", "PlanVerification",
    "hohmann_two_impulse_plan", "verify_two_impulse_plan", "reach_margin", "TANGENCY_FACTOR",
]
