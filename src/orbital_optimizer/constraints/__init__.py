"""Mission constraints and feasibility (Phase 4), SI units.

Physical layer: it imports the constants, the spacecraft model and the
simulation, and nothing from the optimization, problem or mission layers.
How an optimizer treats an infeasible candidate is a separate matter
(optimization.constraint_handling); this package only states, for one
transfer, which limits hold and by how much each one is exceeded.

model: ConstraintValue and ConstraintReport (g <= 0 convention, tolerances,
    normalized violations, hard and soft constraints).
limits: MissionLimits, TransferSummary and evaluate_constraints, the
    catalogue of mass, thrust, delta-v, duration, altitude and propellant.
verification: the same quantities recomputed from a numerical propagation.
"""

from .limits import (
    CONSTRAINT_NAMES, DEFAULT_RELATIVE_TOLERANCE, MissionLimits, TransferSummary, evaluate_constraints,
    limits_from_mapping,
)
from .model import ConstraintReport, ConstraintValue
from .verification import arc_minimum_radius, numerical_transfer_summary, trajectory_minimum_radius

__all__ = [
    "ConstraintValue", "ConstraintReport", "MissionLimits", "TransferSummary", "evaluate_constraints",
    "limits_from_mapping", "CONSTRAINT_NAMES", "DEFAULT_RELATIVE_TOLERANCE",
    "arc_minimum_radius", "trajectory_minimum_radius", "numerical_transfer_summary",
]
