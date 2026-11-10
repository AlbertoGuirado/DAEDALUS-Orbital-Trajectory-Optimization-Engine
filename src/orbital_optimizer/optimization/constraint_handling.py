"""Constraint-handling strategies for comparison-based optimizers (RF-4.07, RF-4.09, RF-4.10).

Every optimizer of this package only *compares* candidates, through a
ranking key (lexicographic, lower is better). A strategy is therefore just a
rule that turns an Evaluation -- raw objective f, validity, total normalized
violation V and per-constraint normalized violations v_j -- into such a key.
The raw objective is never modified, and the point reported as the best is
always chosen by feasibility first, whatever the strategy (see
optimizer.EvaluationBudget).

    FeasibilityRules   feasible (0, f) < infeasible (1, V) < invalid (2, 0)
    DeathPenalty       feasible (0, f) < everything else (1, 0)
    StaticPenalty      valid (0, f + rho * sum_j v_j^2) < invalid (1, 0)

FeasibilityRules is Deb's method (2000): a feasible point beats an
infeasible one, two feasible points compare by objective and two infeasible
points by violation. It has no parameter, which is why it is the default.

DeathPenalty discards every infeasible point. It keeps no information about
how far a point is from feasibility, so it cannot make progress while the
population contains no feasible point.

StaticPenalty is the quadratic exterior penalty of docs/03 section 15.2,
written with normalized violations so that rho has the unit of the objective.
It is an *inexact* penalty: for any finite rho the minimizer of the
penalized function violates an active constraint by an amount of order
(objective sensitivity) / rho. The coefficient must be stated with its
rationale (RF-4.09); there is no default.

An evaluation of a problem without constraints (``constraint_violation`` is
None) is feasible whenever it is valid, so every strategy reduces to
"valid points first, then lower objective" there.
"""

from __future__ import annotations

import math
from typing import Any, Protocol, runtime_checkable

from .problem import Evaluation

Key = tuple[int, float]


@runtime_checkable
class ConstraintHandling(Protocol):
    """Rule that ranks evaluations of a constrained problem."""

    name: str

    def key(self, evaluation: Evaluation) -> Key: ...

    def describe(self) -> dict[str, Any]: ...


class FeasibilityRules:
    """Deb's feasibility rules: no parameter (see the module docstring)."""

    name: str = "feasibility-rules"

    def key(self, evaluation: Evaluation) -> Key:
        if not evaluation.is_valid:
            return (2, 0.0)
        if evaluation.is_feasible:
            return (0, evaluation.objective)
        return (1, float(evaluation.constraint_violation))  # type: ignore[arg-type]

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "reference": "Deb (2000), An efficient constraint handling method for genetic algorithms",
                "ranking": "feasible by objective < infeasible by total normalized violation < invalid"}


class DeathPenalty:
    """Reject every infeasible or invalid point: no parameter, no gradient towards feasibility."""

    name: str = "death-penalty"

    def key(self, evaluation: Evaluation) -> Key:
        return (0, evaluation.objective) if evaluation.is_feasible else (1, 0.0)

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "ranking": "feasible by objective < infeasible or invalid (all tied)"}


class StaticPenalty:
    """Quadratic exterior penalty f + rho * sum_j v_j^2 on valid points.

    coefficient: rho > 0, in the unit of the objective (the v_j are
    dimensionless). rationale: why this value; required and recorded.
    If an evaluation reports only the total violation V, V^2 is used.
    """

    name: str = "static-penalty"

    def __init__(self, coefficient: float, rationale: str) -> None:
        if isinstance(coefficient, bool) or not isinstance(coefficient, (int, float)) \
                or not math.isfinite(coefficient) or coefficient <= 0:
            raise ValueError("coefficient must be a positive finite number in the unit of the objective")
        if not isinstance(rationale, str) or not rationale.strip():
            raise ValueError("rationale must be a nonempty string justifying the coefficient")
        self.coefficient = float(coefficient)
        self.rationale = rationale.strip()

    def penalized(self, evaluation: Evaluation) -> float:
        """f + rho * sum_j v_j^2 of a valid evaluation; the raw objective when it has no violation."""
        if not evaluation.is_constrained:
            return evaluation.objective
        if evaluation.constraints:
            squares = sum(amount * amount for amount in evaluation.constraints.values())
        else:
            squares = float(evaluation.constraint_violation) ** 2  # type: ignore[arg-type]
        return evaluation.objective + self.coefficient * squares

    def key(self, evaluation: Evaluation) -> Key:
        return (0, self.penalized(evaluation)) if evaluation.is_valid else (1, 0.0)

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "coefficient": self.coefficient, "rationale": self.rationale,
                "form": "objective + coefficient * sum_j normalized_violation_j^2",
                "ranking": "valid by penalized objective < invalid",
                "note": "inexact penalty: its minimizer may violate active constraints"}


def feasibility_key(evaluation: Evaluation) -> Key:
    """Ordering used to *report* the best point of a run, independent of the search strategy."""
    return FeasibilityRules().key(evaluation)


def handling_from_mapping(config: Any) -> ConstraintHandling:
    """Strategy from {"type": "feasibility_rules" | "death_penalty" | "static_penalty", ...}.

    ``static_penalty`` requires ``coefficient`` and ``rationale``. Unknown
    types or keys raise ValueError.
    """
    if not hasattr(config, "keys") or "type" not in config:
        raise ValueError("constraint handling must be a mapping with a 'type'")
    kind = config["type"]
    extra = set(config) - {"type"}
    if kind == "feasibility_rules" and not extra:
        return FeasibilityRules()
    if kind == "death_penalty" and not extra:
        return DeathPenalty()
    if kind == "static_penalty" and extra == {"coefficient", "rationale"}:
        return StaticPenalty(config["coefficient"], config["rationale"])
    raise ValueError(f"invalid constraint handling configuration {dict(config)!r}")
