"""Adapters that express orbital problems through the generic optimization interfaces.

This layer may import both the physics (dynamics, orbital, simulation) and
orbital_optimizer.optimization; neither of them imports it.

coplanar_transfer: ascending two-impulse reference problem (Phase 3 baseline).
coplanar_multi_impulse: N impulses, ascending or descending, with the impulse
    instants in the decision vector; N-impulse plans and their verification.
objectives: weighted delta-v / time / energy / risk objective.
constrained: mission constraints (Phase 4) applied to a transfer problem.
references: independent reference solutions for constrained transfers.
"""

from .constrained import ConstrainedTransferProblem
from .coplanar_multi_impulse import CoplanarMultiImpulseProblem, ImpulsivePlan, verify_impulsive_plan
from .coplanar_transfer import (
    TANGENCY_FACTOR, CoplanarTwoImpulseProblem, PlanVerification, TwoImpulsePlan,
    hohmann_two_impulse_plan, reach_margin, verify_two_impulse_plan,
)
from .objectives import TERM_METRICS, transfer_objective
from .references import TransferReference, reach_delta_v, time_limited_reference, two_impulse_closed_form

__all__ = [
    "ConstrainedTransferProblem", "TransferReference", "time_limited_reference", "two_impulse_closed_form",
    "reach_delta_v",
    "CoplanarTwoImpulseProblem", "TwoImpulsePlan", "PlanVerification",
    "hohmann_two_impulse_plan", "verify_two_impulse_plan", "reach_margin", "TANGENCY_FACTOR",
    "CoplanarMultiImpulseProblem", "ImpulsivePlan", "verify_impulsive_plan",
    "transfer_objective", "TERM_METRICS",
]
