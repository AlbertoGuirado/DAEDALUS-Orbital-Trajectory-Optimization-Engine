"""Mission constraints applied to a transfer problem (Phase 4).

``ConstrainedTransferProblem`` wraps a transfer problem and evaluates, for
every valid candidate, the mission constraints of constraints.limits on the
quantities that the base problem reports:

    burn magnitudes        metadata["burn_delta_v_m_s"]
    duration               metadata["transfer_time_s"]
    trajectory minimum     metadata["minimum_radius_m"]   (exact, whole path)

The objective of the evaluation stays the raw objective of the base problem.
The constraints are returned beside it: ``constraint_violation`` (total
normalized violation of the hard constraints), ``constraints`` (per
constraint) and the full report in ``metadata["constraint_report"]``. Nothing
is added to the objective here; how an optimizer ranks infeasible candidates
is the business of ``constraint_handling`` (default: Deb's feasibility rules).

An evaluation that the base problem declares invalid stays invalid: its
constraints are not evaluated, because there is no trajectory to evaluate
them on. A candidate is therefore in one of three explicit states:
invalid (no transfer), valid and infeasible (a transfer that breaks a limit,
with the amounts recorded) or valid and feasible (RF-4.07, RF-4.08).

Thrust: the impulsive model has no thrust level. ``commanded_thrust_n`` is
the thrust at which the mission intends to execute the burns; with it the
thrust constraint is evaluated and ideal burn durations are reported
(Spacecraft.burn_durations_s). Without it the constraint is reported as not
applicable. The durations are diagnostics of how far the burns are from
instantaneous; they are not constrained and no finite-burn dynamics is
simulated (Phase 5).

Mass and thrust limits do not depend on the decision vector: when one of
them is violated every candidate is infeasible, and ``configuration_report``
says so before any optimization.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from ..constraints.limits import MissionLimits, TransferSummary, evaluate_constraints
from ..constraints.model import ConstraintReport
from ..optimization.constraint_handling import ConstraintHandling, FeasibilityRules
from ..optimization.problem import Bounds, Evaluation, EvaluationStatus
from ..spacecraft import Spacecraft

CONFIGURATION_CONSTRAINTS: tuple[str, ...] = ("max_wet_mass", "max_thrust")


class ConstrainedTransferProblem:
    """A transfer problem with mission constraints (see the module docstring).

    problem: a Problem whose valid evaluations carry the three metadata keys
    above (e.g. CoplanarMultiImpulseProblem). spacecraft, limits: the vehicle
    and the mission limits. handling: constraint-handling strategy used by the
    optimizers (FeasibilityRules by default). commanded_thrust_n: thrust level
    of the burns in newtons, or None. Attributes that this class does not
    define (``plan``, ``reference``, counters...) are those of the base problem.
    """

    def __init__(self, problem: Any, spacecraft: Spacecraft, limits: MissionLimits, *,
                 handling: ConstraintHandling | None = None, commanded_thrust_n: float | None = None) -> None:
        if not isinstance(getattr(problem, "bounds", None), Bounds) or not callable(getattr(problem, "evaluate", None)):
            raise ValueError("problem must implement the Problem protocol")
        if not isinstance(spacecraft, Spacecraft):
            raise ValueError("spacecraft must be a Spacecraft")
        if not isinstance(limits, MissionLimits):
            raise ValueError("limits must be MissionLimits")
        if handling is None:
            handling = FeasibilityRules()
        if not callable(getattr(handling, "key", None)) or not callable(getattr(handling, "describe", None)):
            raise ValueError("handling must be a constraint-handling strategy")
        if commanded_thrust_n is not None:
            commanded_thrust_n = TransferSummary((0.0,), 0.0, 1.0, commanded_thrust_n).thrust_n
        self.problem = problem
        self.spacecraft = spacecraft
        self.limits = limits
        self.constraint_handling = handling
        self.commanded_thrust_n = commanded_thrust_n
        self.bounds = problem.bounds
        self.name = f"{problem.name}:constrained"

    def __getattr__(self, name: str) -> Any:
        # Only called for attributes this object lacks: delegate them to the base problem.
        if name == "problem":
            raise AttributeError(name)
        return getattr(self.problem, name)

    def summary(self, metadata: Any) -> TransferSummary:
        """TransferSummary of the metadata of a valid base evaluation."""
        return TransferSummary.from_metadata(metadata, self.commanded_thrust_n)

    def evaluate(self, x: ArrayLike) -> Evaluation:
        """Evaluate the base problem at x and its mission constraints (see the module docstring)."""
        base = self.problem.evaluate(x)
        if not base.is_valid:
            return base
        summary = self.summary(base.metadata)
        report = evaluate_constraints(summary, self.spacecraft, self.limits)
        burns = np.array(summary.burn_delta_v_m_s)
        metadata = dict(base.metadata)
        finite = bool(np.all(np.isfinite(burns)))
        metadata["constraint_report"] = report.to_dict()
        metadata["feasible"] = report.feasible
        metadata["propellant_per_burn_kg"] = self.spacecraft.propellant_per_burn_kg(burns) if finite else None
        metadata["propellant_required_kg"] = report.by_name("propellant").value
        metadata["final_mass_kg"] = self.spacecraft.wet_mass_kg - report.by_name("propellant").value
        metadata["burn_thrust_n"] = self.commanded_thrust_n
        metadata["burn_durations_s"] = (self.spacecraft.burn_durations_s(burns, self.commanded_thrust_n)
                                        if finite and self.commanded_thrust_n is not None else None)
        return Evaluation(base.x, base.objective, EvaluationStatus.VALID, base.message, metadata,
                          report.total_violation, report.normalized())

    def report(self, x: ArrayLike) -> ConstraintReport:
        """ConstraintReport of a decision vector; an invalid one raises ValueError."""
        base = self.problem.evaluate(x)
        if not base.is_valid:
            raise ValueError(f"decision vector is not valid: {base.message}")
        return evaluate_constraints(self.summary(base.metadata), self.spacecraft, self.limits)

    def configuration_report(self) -> ConstraintReport:
        """The constraints that do not depend on the trajectory: wet mass and thrust.

        Evaluated on a transfer without burns; only ``max_wet_mass`` and
        ``max_thrust`` are kept. If it is not feasible, no decision vector can be.
        """
        report = evaluate_constraints(TransferSummary((), 0.0, self.limits.body_radius_m, self.commanded_thrust_n),
                                      self.spacecraft, self.limits)
        kept = tuple(value for value in report.values if value.name in CONFIGURATION_CONSTRAINTS)
        skipped = {name: reason for name, reason in report.not_applicable.items() if name in CONFIGURATION_CONSTRAINTS}
        return ConstraintReport(kept, skipped)

    def describe(self) -> dict[str, Any]:
        """JSON-compatible record of the vehicle, the limits and the handling strategy."""
        return {
            "base_problem": self.problem.name,
            "spacecraft": self.spacecraft.describe(),
            "limits": self.limits.describe(),
            "commanded_thrust_n": self.commanded_thrust_n,
            "constraint_handling": self.constraint_handling.describe(),
        }
