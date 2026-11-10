"""Constraint quantities recomputed from a numerical propagation of a plan.

The optimizer evaluates the path constraint with the exact minimum radius of
each Keplerian arc. This module obtains the same quantities by another
route -- numerical integration of the burns as applied -- so that a reported
feasibility can be checked independently before it is trusted (RF-6B.04).

Minimum radius of a numerically propagated arc:
  1. integrate the arc with RK45 and take the accepted step with the lowest
     |r| (steps are at most ``max_step`` apart, 60 s by default);
  2. integrate again the interval between its two neighbours with
     ``refine_samples`` equally spaced outputs;
  3. fit a parabola through the lowest sample and its neighbours and take
     its vertex.
With a spacing h between refined samples the vertex differs from the true
minimum by O(|d^3 r/dt^3| h^3); for h below one second this is below a
millimetre in low Earth orbit. The dominant error is the global error of the
integrator itself (of the order of 0.1 m over a LEO-GEO transfer with the
default tolerances; docs/phase1-phase2-audit.md), not this search. The
result is a numerical estimate with that accuracy, not a bound.

A minimum at an end point of an arc is the radius of that end point. The
delta-v and the duration are read from the plan itself.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..constants import EARTH_MU
from ..simulation.impulsive import ImpulsiveTrajectory, propagate_impulsive_arcs
from ..simulation.propagator import PropagationResult, propagate
from .limits import TransferSummary


def _parabola_minimum(values: NDArray[np.float64]) -> float:
    """Minimum of the parabola through three equally spaced samples whose middle one is the lowest."""
    left, middle, right = (float(value) for value in values)
    curvature = left - 2 * middle + right
    if curvature <= 0:
        return middle
    offset = 0.5 * (left - right) / curvature
    return middle - 0.25 * (left - right) * offset


def arc_minimum_radius(arc: PropagationResult, *, refine_samples: int = 201,
                       method: Literal["RK45", "RK4"] = "RK45", max_step: float = 60.0,
                       rtol: float = 1e-10, atol: ArrayLike | None = None) -> float:
    """Minimum |r| in metres over a propagated arc (see the module docstring).

    ``arc`` holds the accepted steps of a numerical propagation. The interval
    around its lowest step is integrated again with the same options and
    ``refine_samples`` >= 3 outputs.
    """
    if isinstance(refine_samples, bool) or not isinstance(refine_samples, (int, np.integer)) or refine_samples < 3:
        raise ValueError("refine_samples must be an integer >= 3")
    radii = np.linalg.norm(arc.states[:, :3], axis=1)
    lowest = int(np.argmin(radii))
    if lowest == 0 or lowest == len(radii) - 1:
        # The radius is monotonic next to this end point up to the step size: check the adjacent step.
        first, last = (0, 1) if lowest == 0 else (len(radii) - 2, len(radii) - 1)
    else:
        first, last = lowest - 1, lowest + 1
    span = float(arc.times[last] - arc.times[first])
    refined = propagate(arc.states[first], span, arc.mu, method=method,
                        sample_times=np.linspace(0.0, span, refine_samples), max_step=max_step, rtol=rtol, atol=atol)
    fine = np.linalg.norm(refined.states[:, :3], axis=1)
    index = int(np.argmin(fine))
    if index == 0 or index == len(fine) - 1:
        return float(fine[index])
    return _parabola_minimum(fine[index - 1:index + 2])


def trajectory_minimum_radius(trajectory: ImpulsiveTrajectory, **options) -> float:
    """Minimum |r| in metres over every arc of a numerically propagated trajectory."""
    if not trajectory.arcs:
        raise ValueError("the trajectory has no arc")
    return min(arc_minimum_radius(arc, **options) for arc in trajectory.arcs)


def numerical_transfer_summary(
    initial_state: ArrayLike,
    impulses: Sequence[tuple[float, ArrayLike]],
    mu: float = EARTH_MU,
    *,
    thrust_n: float | None = None,
    refine_samples: int = 201,
    method: Literal["RK45", "RK4"] = "RK45",
    max_step: float = 60.0,
    rtol: float = 1e-10,
    atol: ArrayLike | None = None,
) -> TransferSummary:
    """TransferSummary of a plan from a numerical propagation up to its last burn.

    ``initial_state`` (6,) in m and m/s; ``impulses`` as for
    simulation.impulsive.propagate_with_impulses, with at least one burn
    after t = 0 (the flight must have a duration). The burn magnitudes and
    the duration come from ``impulses``; the minimum radius from the
    integrated arcs. Bad input raises ValueError; integration failure raises
    RuntimeError.
    """
    burns = list(impulses)
    if not burns:
        raise ValueError("at least one impulse is required")
    duration = float(burns[-1][0])
    trajectory = propagate_impulsive_arcs(initial_state, burns, duration, mu, method=method,
                                          max_step=max_step, rtol=rtol, atol=atol)
    minimum = trajectory_minimum_radius(trajectory, refine_samples=refine_samples, method=method,
                                        max_step=max_step, rtol=rtol, atol=atol)
    magnitudes = tuple(float(np.linalg.norm(vector)) for vector in trajectory.burn_vectors)
    return TransferSummary(magnitudes, duration, minimum, thrust_n)
