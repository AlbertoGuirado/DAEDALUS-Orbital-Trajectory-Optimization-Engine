"""Optimization problem contracts on a physics-independent decision space.

A problem declares physical decision bounds, with one unit label per variable.
Optimizers search the unit hypercube instead, through the affine map

    u = (x - lower) / (upper - lower),   x = lower + u * (upper - lower),

so variables of very different scales are treated alike (RF-3.08). This module
imports only NumPy and the standard library; it knows nothing about orbits.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType
from typing import Any, Protocol, runtime_checkable

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _real_array(value: ArrayLike, name: str) -> NDArray[np.float64]:
    """Return a float64 copy of real numerical data; reject other dtypes."""
    try:
        array = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain real numbers") from exc
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numbers")
    return np.array(array, dtype=float, copy=True)


def _read_only(array: NDArray) -> NDArray:
    array.setflags(write=False)
    return array


def _labels(values: Sequence[str] | None, default: Sequence[str], name: str) -> tuple[str, ...]:
    if values is None:
        return tuple(default)
    if isinstance(values, str):
        raise ValueError(f"{name} must be a sequence of strings, not a single string")
    labels = tuple(values)
    if len(labels) != len(default):
        raise ValueError(f"{name} must have one entry per decision variable")
    if not all(isinstance(label, str) and label for label in labels):
        raise ValueError(f"{name} must contain nonempty strings")
    return labels


@dataclass(frozen=True, eq=False)
class Bounds:
    """Box bounds of a physical decision vector of dimension D >= 1.

    ``lower`` and ``upper`` are finite (D,) arrays in the units named by
    ``units``, with lower < upper componentwise. ``periodic`` marks variables
    whose domain is the half-open period [lower, upper), e.g. an angle in
    radians; the optimizer wraps those instead of reflecting them, and
    ``from_unit``, ``contains`` and ``validate`` never accept or return
    ``upper`` for them. Defaults: names ``x0..x{D-1}``, units
    ``"dimensionless"``, no periodic variable. The arrays are independent
    copies of the inputs and are read-only. Invalid bounds raise ValueError.
    """

    lower: NDArray[np.float64]
    upper: NDArray[np.float64]
    names: tuple[str, ...] | None = None
    units: tuple[str, ...] | None = None
    periodic: NDArray[np.bool_] | None = None

    def __post_init__(self) -> None:
        lower = _real_array(self.lower, "lower")
        upper = _real_array(self.upper, "upper")
        if lower.ndim != 1 or lower.size < 1 or lower.shape != upper.shape:
            raise ValueError("lower and upper must be nonempty one-dimensional arrays of equal shape")
        if not (np.all(np.isfinite(lower)) and np.all(np.isfinite(upper))):
            raise ValueError("bounds must be finite")
        if np.any(lower >= upper):
            raise ValueError("every lower bound must be strictly less than its upper bound")
        with np.errstate(over="ignore"):
            width = upper - lower
        if not np.all(np.isfinite(width)):
            raise ValueError("upper - lower must be finite")
        dimension = lower.size
        names = _labels(self.names, [f"x{j}" for j in range(dimension)], "names")
        if len(set(names)) != dimension:
            raise ValueError("names must be unique")
        units = _labels(self.units, ["dimensionless"] * dimension, "units")
        if self.periodic is None:
            periodic = np.zeros(dimension, dtype=bool)
        else:
            raw = np.asarray(self.periodic)
            if raw.dtype != np.bool_ or raw.shape != (dimension,):
                raise ValueError("periodic must be a boolean array of shape (D,)")
            periodic = np.array(raw, dtype=bool, copy=True)
        object.__setattr__(self, "lower", _read_only(lower))
        object.__setattr__(self, "upper", _read_only(upper))
        object.__setattr__(self, "names", names)
        object.__setattr__(self, "units", units)
        object.__setattr__(self, "periodic", _read_only(periodic))

    @property
    def dimension(self) -> int:
        """Number D of decision variables."""
        return int(self.lower.size)

    @property
    def width(self) -> NDArray[np.float64]:
        """upper - lower, shape (D,), in each variable's unit."""
        return self.upper - self.lower

    def to_unit(self, x: ArrayLike) -> NDArray[np.float64]:
        """Map physical x, shape (D,) or (N,D), to u = (x - lower)/(upper - lower).

        No bound check is made: values outside the box map outside [0,1].
        """
        values = _real_array(x, "x")
        if values.shape[-1:] != (self.dimension,) or values.ndim not in (1, 2):
            raise ValueError("x must have shape (D,) or (N,D)")
        return (values - self.lower) / self.width

    def from_unit(self, u: ArrayLike) -> NDArray[np.float64]:
        """Map u in [0,1], shape (D,) or (N,D), to x = lower + u*(upper - lower).

        Non-periodic variables are clipped to [lower, upper] only to absorb
        floating-point rounding of the affine map. Periodic variables always
        land in the half-open period [lower, upper): whenever the rounded map
        gives x >= upper, lower is returned instead. Because lower + 1*(upper
        - lower) is itself rounded, u = 1 may give either lower (when the sum
        rounds to upper) or a value just below upper (when it rounds down);
        u slightly below 1 may likewise round up to upper and then give lower.
        lower and upper are the same point of the period, so the replacement
        keeps the physical point and picks its canonical representative;
        stepping to nextafter(upper, lower) would instead move the point by
        one ulp. u outside [0,1] raises ValueError.
        """
        values = _real_array(u, "u")
        if values.shape[-1:] != (self.dimension,) or values.ndim not in (1, 2):
            raise ValueError("u must have shape (D,) or (N,D)")
        if not np.all(np.isfinite(values)) or np.any(values < 0.0) or np.any(values > 1.0):
            raise ValueError("u must be finite and lie in [0,1]")
        x = np.clip(self.lower + values * self.width, self.lower, self.upper)
        return np.where(self.periodic & (x >= self.upper), self.lower, x)

    def _inside(self, values: NDArray[np.float64]) -> bool:
        above = np.where(self.periodic, values >= self.upper, values > self.upper)
        return bool(np.all(values >= self.lower) and not np.any(above))

    def contains(self, x: ArrayLike) -> bool:
        """True if x, shape (D,), is finite and in the box.

        The box is [lower, upper] for non-periodic and [lower, upper) for
        periodic variables.
        """
        values = _real_array(x, "x")
        return bool(values.shape == (self.dimension,) and np.all(np.isfinite(values)) and self._inside(values))

    def validate(self, x: ArrayLike) -> NDArray[np.float64]:
        """Return an independent read-only copy of x, shape (D,), inside the box.

        Wrong shape, non-finite values or values outside the box ([lower, upper],
        or [lower, upper) for periodic variables) raise ValueError: they are
        caller errors, not invalid evaluations.
        """
        values = _real_array(x, "x")
        if values.shape != (self.dimension,):
            raise ValueError(f"x must have shape ({self.dimension},)")
        if not np.all(np.isfinite(values)):
            raise ValueError("x must be finite")
        if not self._inside(values):
            raise ValueError("x must lie within the bounds")
        return _read_only(values)

    def describe(self) -> dict[str, Any]:
        """JSON-compatible description of the bounds."""
        return {
            "dimension": self.dimension,
            "names": list(self.names),
            "units": list(self.units),
            "lower": self.lower.tolist(),
            "upper": self.upper.tolist(),
            "periodic": self.periodic.tolist(),
        }


class EvaluationStatus(Enum):
    """Outcome of one objective evaluation.

    VALID: finite objective. INVALID_DOMAIN: the candidate lies outside the
    domain where the problem's evaluation is defined. NON_FINITE: the objective
    was NaN or +/-inf. NUMERICAL_FAILURE: a declared computational exception.
    Only VALID evaluations can become the optimizer's best solution.
    """

    VALID = "valid"
    INVALID_DOMAIN = "invalid_domain"
    NON_FINITE = "non_finite"
    NUMERICAL_FAILURE = "numerical_failure"

    @property
    def is_valid(self) -> bool:
        return self is EvaluationStatus.VALID


@dataclass(frozen=True, eq=False)
class Evaluation:
    """Result of evaluating one physical decision vector x, shape (D,).

    ``objective`` is in the problem's objective unit. It must be finite when
    ``status`` is VALID; otherwise it records the raw value (NaN by default).
    ``metadata`` is a read-only mapping of problem-specific diagnostics.

    Constraints (Phase 4). ``constraint_violation`` is None for a problem
    without constraints. A constrained problem sets it to the total
    normalized violation V >= 0 of its hard constraints (dimensionless; 0
    means feasible, +inf a constraint that could not be evaluated) and
    ``constraints`` to the normalized violation of each one by name. The
    objective is always the raw objective: a violation is never folded into
    it. How violations affect the search is decided by the optimizer's
    constraint handling, not here.
    """

    x: NDArray[np.float64]
    objective: float
    status: EvaluationStatus
    message: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)
    constraint_violation: float | None = None
    constraints: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.status, EvaluationStatus):
            raise ValueError("status must be an EvaluationStatus")
        x = _real_array(self.x, "x")
        if x.ndim != 1:
            raise ValueError("x must be one-dimensional")
        objective = _real_array(self.objective, "objective")
        if objective.shape != ():
            raise ValueError("objective must be a scalar")
        value = float(objective)
        if self.status.is_valid and not math.isfinite(value):
            raise ValueError("a VALID evaluation requires a finite objective")
        if not isinstance(self.message, str):
            raise ValueError("message must be a string")
        violation = self.constraint_violation
        if violation is not None:
            if isinstance(violation, bool) or not isinstance(violation, (int, float, np.integer, np.floating)) \
                    or math.isnan(violation) or violation < 0:
                raise ValueError("constraint_violation must be None or a real number >= 0")
            violation = float(violation)
        constraints = {}
        for name, amount in dict(self.constraints).items():
            if not isinstance(name, str) or isinstance(amount, bool) \
                    or not isinstance(amount, (int, float, np.integer, np.floating)) \
                    or math.isnan(amount) or amount < 0:
                raise ValueError("constraints must map names to normalized violations >= 0")
            constraints[name] = float(amount)
        if constraints and violation is None:
            raise ValueError("constraints require a constraint_violation")
        object.__setattr__(self, "x", _read_only(x))
        object.__setattr__(self, "objective", value)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))
        object.__setattr__(self, "constraint_violation", violation)
        object.__setattr__(self, "constraints", MappingProxyType(constraints))

    @property
    def is_valid(self) -> bool:
        return self.status.is_valid

    @property
    def is_constrained(self) -> bool:
        """True if the problem reported a constraint violation (possibly zero)."""
        return self.constraint_violation is not None

    @property
    def is_feasible(self) -> bool:
        """VALID and without constraint violation; for an unconstrained problem, the same as valid."""
        return self.status.is_valid and (self.constraint_violation is None or self.constraint_violation == 0.0)


@runtime_checkable
class Problem(Protocol):
    """Single-objective minimization problem over box bounds.

    ``evaluate`` receives a validated, read-only physical vector of shape (D,)
    inside ``bounds`` and returns an Evaluation. Domain or numerical failures
    that the problem anticipates are reported through the status; any other
    exception propagates and aborts the optimization run.
    """

    name: str
    bounds: Bounds

    def evaluate(self, x: NDArray[np.float64]) -> Evaluation: ...


class FunctionProblem:
    """Wrap a scalar function f(x) -> float as a Problem.

    Non-finite values (NaN, +inf, -inf) yield NON_FINITE evaluations that keep
    the raw value. Exceptions propagate unless their type is listed in
    ``invalid_exceptions``; those become NUMERICAL_FAILURE evaluations. A
    return value that is not a real scalar raises ValueError.
    """

    def __init__(
        self,
        function: Callable[[NDArray[np.float64]], float],
        bounds: Bounds,
        name: str = "function",
        *,
        invalid_exceptions: tuple[type[Exception], ...] = (),
    ) -> None:
        if not callable(function):
            raise ValueError("function must be callable")
        if not isinstance(bounds, Bounds):
            raise ValueError("bounds must be a Bounds instance")
        if not isinstance(name, str) or not name:
            raise ValueError("name must be a nonempty string")
        if not isinstance(invalid_exceptions, tuple) or not all(
            isinstance(kind, type) and issubclass(kind, Exception) for kind in invalid_exceptions
        ):
            raise ValueError("invalid_exceptions must be a tuple of Exception subclasses")
        self.function = function
        self.bounds = bounds
        self.name = name
        self.invalid_exceptions = invalid_exceptions

    def evaluate(self, x: ArrayLike) -> Evaluation:
        """Evaluate f at x, shape (D,), inside the bounds (else ValueError)."""
        point = self.bounds.validate(x)
        try:
            raw = self.function(point)
        except self.invalid_exceptions as exc:
            return Evaluation(point, math.nan, EvaluationStatus.NUMERICAL_FAILURE, f"{type(exc).__name__}: {exc}")
        value = np.asarray(raw)
        if value.shape != () or value.dtype.kind not in "iuf":
            raise ValueError("the objective function must return a real scalar")
        objective = float(value)
        if not math.isfinite(objective):
            return Evaluation(point, objective, EvaluationStatus.NON_FINITE, f"objective is {objective}")
        return Evaluation(point, objective, EvaluationStatus.VALID)
