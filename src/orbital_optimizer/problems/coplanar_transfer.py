"""Coplanar two-impulse transfer between circular orbits as an optimization Problem.

Adapter layer (Phase 3): it imports the Phase 1-2 physics and the generic
optimization interfaces; neither of those layers imports this module.

Model and decision vector
-------------------------
Two-body point mass with constant mu; concentric coplanar circular orbits of
radii r1 < r2 in the xy plane, prograde; two instantaneous impulses; free
transfer time; the target is the orbit, not a point on it. By the rotational
symmetry of the initial circle the departure point is fixed at
r0 = (r1, 0, 0), v0 = (0, v_c1, 0), t1 = 0, with v_c = sqrt(mu/r). Decision
(SI): x = (dv1, alpha1),

    dv1    in [0, (sqrt(2) - 1) v_c1]  m/s  magnitude of the first impulse
    alpha1 in [-pi, pi)  rad, periodic  angle from the local velocity
                                         direction, positive radially outward
    first impulse  dv1_vec = dv1 (sin(alpha1), cos(alpha1), 0)

Below the upper bound |v0 + dv1_vec| < sqrt(2) v_c1 for every alpha1, so all
candidates are bound, and v_t = v_c1 + dv1 cos(alpha1) > 0 keeps them
prograde. The second impulse is not a decision variable: at the first
outward crossing of r2 it is the unique burn that circularizes there,
dv2_vec = sqrt(mu/r2) (z_hat x r_hat) - v(tau). The transfer time tau is derived.

Objective (m/s): J(x) = |dv1_vec| + |dv2_vec| (orbital.mechanics.total_delta_v).
Within this restricted class (two impulses, coplanar, circular, free time,
no phasing) the classical Hohmann transfer is the minimum; it is the reference
J_H, t_H from orbital.transfers.hohmann_transfer. No claim is made about
other numbers of impulses, constraints or perturbations.

Evaluation chain
----------------
1. Analytic domain predicates, before any physics call: bound and not
   numerically parabolic (the threshold of state_to_elements), not radial,
   and the transfer conic reaches r2 (tangency test below). Failing
   candidates are INVALID_DOMAIN.
2. state_to_elements(y1) gives a, e and the departure true anomaly nu1.
   cos(nu2) = (p/r2 - 1)/e with p = a(1-e^2); nu2 = arccos(clip(cos(nu2),
   -1, 1)) is the outward crossing, which is always the first one reached
   because the departure radius r1 is below r2. Clipping happens only at
   accepted tangency and its amount is recorded.
3. tau = elliptic_time_of_flight(nu1, nu2, a, e, mu) and
   propagate_kepler(y1, tau) give the arrival state (analytic f/g
   propagation, no RK45 in the objective).
4. dv2_vec as above and J. A ValueError or RuntimeError raised by the physics
   calls is counted (``captured_value_errors``, ``captured_runtime_errors``)
   and reported as INVALID_DOMAIN or NUMERICAL_FAILURE; in normal use the
   predicates make both counts zero. Other exceptions propagate.

Tangency tolerance
------------------
The exact Hohmann point lies on the feasibility boundary r_a = r2, so a
strict floating-point test can reject it (measured: r_a - r2 = -4.5e-8 m at
the float Hohmann dv1 for LEO->GEO via state_to_elements). The conic through
the departure point (r_p <= r1 < r2) reaches r2 iff the radial speed there is
real, i.e. iff

    Q = 2 eps r2^2 + 2 mu r2 - h^2 >= 0,  eps = v^2/2 - mu/r1,  h = r1 v_t.

Q is a plain sum of three terms with no square root, so its forward error
has a scale-free bound in terms of S = 2 r2^2 (v^2/2 + mu/r1) + 2 mu r2 + h^2,
the sum of the magnitudes before cancellation, with unit roundoff
u = 2^-53:
  (a) arithmetic: the longest path from the inputs to Q has about 10
      roundings (v_c1 ~1.5, v_t 1, v^2 2, eps 1, times r2^2 2, h^2 2,
      two additions), so |Q_computed - Q| <~ 10 u S;
  (b) representation of the decision nearest the exact tangent: the float
      Hohmann dv1 from hohmann_transfer carries an absolute error of about
      3u v_p + 1.5u v_c1 <= 4.5u v, and |dQ/dv| dv <= 2 S dv/v <= 9u S.
The sum is about 19u S; the tolerance takes the next power of two for the
neglected second-order terms:

    accept iff Q_computed >= -K u S,  K = 32.

K depends on neither r2/r1 nor scale; it was fixed by this argument, not
fitted to tests. Confirmation sweep (2026-11-03, 3 values of mu x 3 scales x 40
ratios r2/r1 in [1.0001, 100]): at the float Hohmann dv1 the measured
|Q|/(u S) was at most 1.56, and Q was negative in 153 of 360 cases, so a strict
test would have rejected them.

Admitted infeasible band (measured by bisection, alpha1 = 0): for LEO->GEO,
dv1 down to 93 ulps (4.2e-11 m/s) below the float Hohmann value is accepted,
equivalent to an apoapsis deficit of 6.3e-14 relative (2.7e-6 m); J at the
band edge is J_H - 1.4e-10 m/s. Because Q barely depends on dv1 when r2 ~ r1,
the band in dv1 widens for near-equal radii (r1 = 6 678 km, mu = EARTH_MU:
relative 2.2e-10 at r2/r1 = 1.01, 2.1e-6 at 1.0001), while the
apoapsis-equivalent deficit stayed below 3e-10 at the 12 ratios measured.
Thresholds such as "J >= J_H - eps" must account for this when r2/r1 -> 1.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..dynamics.two_body import EARTH_MU, real_array, validate_mu, validate_state
from ..optimization.problem import Bounds, Evaluation, EvaluationStatus
from ..orbital.elements import SINGULAR_TOL, state_to_elements
from ..orbital.kepler import elliptic_time_of_flight
from ..orbital.mechanics import orbital_period, total_delta_v
from ..orbital.transfers import hohmann_transfer
from ..simulation.impulsive import propagate_with_impulses
from ..simulation.kepler_propagator import propagate_kepler
from ..simulation.propagator import propagate

UNIT_ROUNDOFF: float = float(np.finfo(float).eps) / 2
TANGENCY_FACTOR: float = 32.0
# Phase 2 arrival criteria (tests/test_impulses.py), reused for verification.
VERIFICATION_RADIUS_TOLERANCE_M: float = 1.0
VERIFICATION_ECCENTRICITY_TOLERANCE: float = 1e-6


def _radius(value: float, name: str) -> float:
    radius = real_array(value, name)
    if radius.shape != () or radius <= 0:
        raise ValueError(f"{name} must be a positive finite scalar in metres")
    return float(radius)


def _frozen_vector(value: ArrayLike, name: str, size: int) -> NDArray[np.float64]:
    array = real_array(value, name)
    if array.shape != (size,):
        raise ValueError(f"{name} must have shape ({size},)")
    array.setflags(write=False)
    return array


def _hohmann_reference(r1: float, r2: float, mu: float):
    """hohmann_transfer without NumPy warnings; callers check finiteness explicitly.

    hohmann_transfer (Phase 2) evaluates abs(v_perigee - v_circular) outside
    its errstate block, so at extreme scales inf - inf emits a RuntimeWarning
    and yields NaN. That function is not modified here; the warning is
    suppressed for this call and the NaN is rejected by the caller.
    """
    with np.errstate(over="ignore", invalid="ignore"):
        return hohmann_transfer(r1, r2, mu)


def reach_margin(initial_radius: float, final_radius: float, mu: float,
                 tangential_speed: float, radial_speed: float) -> tuple[float, float]:
    """Return (Q, S) of the tangency test; units m^4/s^2.

    Q = 2 eps r2^2 + 2 mu r2 - h^2 and S = 2 r2^2 (v^2/2 + mu/r1) + 2 mu r2 + h^2
    for the state at r1 with velocity components v_t (along the circular
    velocity) and v_r (radial), in m/s. Q >= 0 iff the conic reaches r2.
    """
    r1, r2 = initial_radius, final_radius
    speed2 = tangential_speed * tangential_speed + radial_speed * radial_speed
    energy = speed2 / 2 - mu / r1
    h2 = (r1 * tangential_speed) ** 2
    q = 2 * energy * r2 * r2 + 2 * mu * r2 - h2
    s = 2 * r2 * r2 * (speed2 / 2 + mu / r1) + 2 * mu * r2 + h2
    return float(q), float(s)


@dataclass(frozen=True, eq=False)
class TwoImpulsePlan:
    """Oriented two-impulse plan in the inertial frame, SI units.

    initial_state (6,) [m, m/s] at t=0; first_burn (3,) m/s applied at t=0;
    second_burn (3,) m/s applied at transfer_time (s, > 0); target circular
    orbit of final_radius (m) about mu (m^3/s^2). Arrays are read-only.
    """

    initial_state: NDArray[np.float64]
    first_burn: NDArray[np.float64]
    second_burn: NDArray[np.float64]
    transfer_time: float
    final_radius: float
    mu: float = EARTH_MU

    def __post_init__(self) -> None:
        state = validate_state(self.initial_state)
        state.setflags(write=False)
        object.__setattr__(self, "initial_state", state)
        object.__setattr__(self, "first_burn", _frozen_vector(self.first_burn, "first_burn", 3))
        object.__setattr__(self, "second_burn", _frozen_vector(self.second_burn, "second_burn", 3))
        time = real_array(self.transfer_time, "transfer_time")
        if time.shape != () or time <= 0:
            raise ValueError("transfer_time must be a positive finite scalar in seconds")
        object.__setattr__(self, "transfer_time", float(time))
        object.__setattr__(self, "final_radius", _radius(self.final_radius, "final_radius"))
        object.__setattr__(self, "mu", validate_mu(self.mu))

    @property
    def total_delta_v(self) -> float:
        """|first_burn| + |second_burn|, m/s."""
        return total_delta_v([self.first_burn, self.second_burn])

    def impulses(self) -> list[tuple[float, NDArray[np.float64]]]:
        """[(0, first_burn), (transfer_time, second_burn)] for propagate_with_impulses."""
        return [(0.0, self.first_burn.copy()), (self.transfer_time, self.second_burn.copy())]


def hohmann_two_impulse_plan(initial_radius: float, final_radius: float, mu: float = EARTH_MU) -> TwoImpulsePlan:
    """Oriented Hohmann plan, ascending or descending, from hohmann_transfer magnitudes.

    Departure at (r1, 0, 0) with prograde velocity (0, v_c1, 0). With
    s = sign(r2 - r1) the burns are (0, s dv1, 0) at t = 0 and (0, -s dv2, 0)
    at t_H, where the vehicle is at (-r2, 0, 0) moving along -y. Equal radii
    raise ValueError (no transfer); a non-finite reference raises RuntimeError.
    """
    r1, r2 = _radius(initial_radius, "initial_radius"), _radius(final_radius, "final_radius")
    gm = validate_mu(mu)
    if r1 == r2:
        raise ValueError("initial and final radii must differ")
    reference = _hohmann_reference(r1, r2, gm)
    values = (reference.first_delta_v, reference.second_delta_v, reference.transfer_time)
    if not all(np.isfinite(values)):
        raise RuntimeError("Hohmann reference is not finite at this numerical scale")
    sign = 1.0 if r2 > r1 else -1.0
    state = np.array([r1, 0.0, 0.0, 0.0, np.sqrt(gm / r1), 0.0])
    return TwoImpulsePlan(state, [0.0, sign * reference.first_delta_v, 0.0],
                          [0.0, -sign * reference.second_delta_v, 0.0], reference.transfer_time, r2, gm)


@dataclass(frozen=True, eq=False)
class PlanVerification:
    """Outcome of verify_two_impulse_plan. ``passed`` is True iff every check holds.

    checks: name -> bool. metrics (SI): final_semi_major_axis_error_m,
    final_eccentricity, final_radius_error_m, coast_max_radius_error_m,
    coast_final_eccentricity (absent if not reached). final_state: (6,)
    after the second burn, or None if propagation failed.
    """

    passed: bool
    checks: Mapping[str, bool]
    metrics: Mapping[str, float]
    final_state: NDArray[np.float64] | None
    message: str
    tolerances: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in ("checks", "metrics", "tolerances"):
            object.__setattr__(self, name, MappingProxyType(dict(getattr(self, name))))
        if self.final_state is not None:
            state = np.array(self.final_state, dtype=float, copy=True)
            state.setflags(write=False)
            object.__setattr__(self, "final_state", state)


def verify_two_impulse_plan(plan: TwoImpulsePlan) -> PlanVerification:
    """Independently propagate a plan with RK45 and check arrival on the target circle.

    1. propagate_with_impulses(initial_state, [(0, b1), (tau, b2)], tau, mu)
       with its RK45 defaults (rtol 1e-10, atol 1e-3 m / 1e-6 m/s, max_step 60 s).
    2. Explicit finiteness and state_to_elements on the final state, because
       that wrapper does not revalidate a burn applied at the final instant.
    3. |a - r2| <= 1 m, e < 1e-6, | |r| - r2 | <= 1 m: the Phase 2 arrival
       criteria of tests/test_impulses.py.
    4. Coast of one target period with propagate (RK45 defaults): max over the
       accepted steps of | |r| - r2 | <= 1 m and final e < 1e-6.
    Propagation or element errors (ValueError/RuntimeError) make the
    verification fail with a message instead of raising.
    """
    if not isinstance(plan, TwoImpulsePlan):
        raise ValueError("plan must be a TwoImpulsePlan")
    r2, gm = plan.final_radius, plan.mu
    radius_tol, ecc_tol = VERIFICATION_RADIUS_TOLERANCE_M, VERIFICATION_ECCENTRICITY_TOLERANCE
    tolerances = {"radius_m": radius_tol, "eccentricity": ecc_tol}
    checks: dict[str, bool] = {}
    metrics: dict[str, float] = {}

    def failed(message: str, final: NDArray | None = None) -> PlanVerification:
        return PlanVerification(False, checks, metrics, final, message, tolerances)

    try:
        final = propagate_with_impulses(plan.initial_state, plan.impulses(), plan.transfer_time, gm)
    except (ValueError, RuntimeError) as exc:
        return failed(f"transfer propagation failed: {type(exc).__name__}: {exc}")
    checks["final_state_finite"] = bool(np.all(np.isfinite(final)))
    if not checks["final_state_finite"]:
        return failed("final state is not finite", final)
    try:
        elements = state_to_elements(final, gm)
    except (ValueError, RuntimeError) as exc:
        return failed(f"final state is not a supported ellipse: {exc}", final)
    metrics["final_semi_major_axis_error_m"] = elements.semi_major_axis - r2
    metrics["final_eccentricity"] = elements.eccentricity
    metrics["final_radius_error_m"] = float(np.linalg.norm(final[:3]) - r2)
    checks["final_semi_major_axis"] = abs(metrics["final_semi_major_axis_error_m"]) <= radius_tol
    checks["final_eccentricity"] = elements.eccentricity < ecc_tol
    checks["final_radius"] = abs(metrics["final_radius_error_m"]) <= radius_tol
    try:
        coast = propagate(final, orbital_period(r2, gm), gm)
        coast_elements = state_to_elements(coast.final_state, gm)
    except (ValueError, RuntimeError) as exc:
        return failed(f"coast propagation failed: {type(exc).__name__}: {exc}", final)
    radii = np.linalg.norm(coast.states[:, :3], axis=1)
    metrics["coast_max_radius_error_m"] = float(np.max(np.abs(radii - r2)))
    metrics["coast_final_eccentricity"] = coast_elements.eccentricity
    checks["coast_radius"] = metrics["coast_max_radius_error_m"] <= radius_tol
    checks["coast_eccentricity"] = coast_elements.eccentricity < ecc_tol
    passed = all(checks.values())
    failing = [name for name, ok in checks.items() if not ok]
    message = "all checks passed" if passed else "failed checks: " + ", ".join(failing)
    return PlanVerification(passed, checks, metrics, final, message, tolerances)


class CoplanarTwoImpulseProblem:
    """Ascending coplanar circular-to-circular two-impulse transfer (see module docstring).

    Decision x = (dv1 [m/s], alpha1 [rad, periodic]); objective J in m/s.
    Requires 0 < r1 < r2 (descending transfers are verification-only in
    Phase 3: use hohmann_two_impulse_plan); r1 == r2 or r2 < r1 raise
    ValueError. A non-finite Hohmann reference or bound raises RuntimeError.
    The counters ``captured_value_errors`` and ``captured_runtime_errors``
    accumulate over the lifetime of the instance.
    """

    def __init__(self, initial_radius: float, final_radius: float, mu: float = EARTH_MU) -> None:
        r1, r2 = _radius(initial_radius, "initial_radius"), _radius(final_radius, "final_radius")
        gm = validate_mu(mu)
        if r1 == r2:
            raise ValueError("initial and final radii must differ")
        if r2 < r1:
            raise ValueError("only ascending transfers (r2 > r1) are optimization problems in Phase 3")
        reference = _hohmann_reference(r1, r2, gm)
        with np.errstate(over="ignore", invalid="ignore"):
            circular_speed = float(np.sqrt(gm / r1))
            upper = float((np.sqrt(2.0) - 1.0) * circular_speed)
        values = (reference.first_delta_v, reference.second_delta_v, reference.total_delta_v,
                  reference.transfer_time, circular_speed, upper)
        if not all(np.isfinite(values)) or upper <= 0:
            raise RuntimeError("Hohmann reference or decision bounds are not finite at this numerical scale")
        self.name = "coplanar-two-impulse"
        self.initial_radius, self.final_radius, self.mu = r1, r2, gm
        self.reference = reference
        self.bounds = Bounds([0.0, -np.pi], [upper, np.pi], names=("delta_v1", "alpha1"),
                             units=("m/s", "rad"), periodic=np.array([False, True]))
        self._circular_speed = circular_speed
        self._initial_state = np.array([r1, 0.0, 0.0, 0.0, circular_speed, 0.0])
        self._initial_state.setflags(write=False)
        self.captured_value_errors = 0
        self.captured_runtime_errors = 0

    @property
    def initial_state(self) -> NDArray[np.float64]:
        """Departure state (6,) before the first burn, a copy."""
        return self._initial_state.copy()

    @property
    def hohmann_decision(self) -> NDArray[np.float64]:
        """Decision vector of the Hohmann transfer: (first_delta_v, 0)."""
        return np.array([self.reference.first_delta_v, 0.0])

    def _invalid(self, x: NDArray, status: EvaluationStatus, message: str, **metadata: Any) -> Evaluation:
        return Evaluation(x, np.nan, status, message, metadata)

    def evaluate(self, x: ArrayLike) -> Evaluation:
        """Evaluate J(x) in m/s; x = (dv1 m/s, alpha1 rad) inside the bounds (else ValueError)."""
        point = self.bounds.validate(x)
        r1, r2, gm = self.initial_radius, self.final_radius, self.mu
        dv1, alpha = float(point[0]), float(point[1])
        burn1 = np.array([dv1 * np.sin(alpha), dv1 * np.cos(alpha), 0.0])
        tangential = self._circular_speed + burn1[1]
        radial = burn1[0]

        # 1. Analytic predicates; they mirror the domain of state_to_elements.
        speed2 = tangential * tangential + radial * radial
        energy = speed2 / 2 - gm / r1
        if not energy < -SINGULAR_TOL * gm / r1:
            return self._invalid(point, EvaluationStatus.INVALID_DOMAIN,
                                 "transfer orbit is unbound or numerically parabolic",
                                 specific_energy_j_kg=energy)
        if not tangential > SINGULAR_TOL * np.sqrt(speed2):
            return self._invalid(point, EvaluationStatus.INVALID_DOMAIN,
                                 "transfer orbit is radial or retrograde", tangential_speed_m_s=tangential)
        q, s = reach_margin(r1, r2, gm, tangential, radial)
        if q < -TANGENCY_FACTOR * UNIT_ROUNDOFF * s:
            return self._invalid(point, EvaluationStatus.INVALID_DOMAIN,
                                 "transfer orbit does not reach the target radius",
                                 reach_margin_relative=q / s)

        # 2-4. Physics: elements, crossing, analytic propagation, second burn.
        state = self._initial_state + np.concatenate(([0.0, 0.0, 0.0], burn1))
        try:
            elements = state_to_elements(state, gm)
            a, e, nu1 = elements.semi_major_axis, elements.eccentricity, elements.true_anomaly
            if not e > 0:
                raise RuntimeError("transfer orbit is numerically circular and cannot cross r2")
            cos_nu2 = (a * (1 - e) * (1 + e) / r2 - 1) / e
            nu2 = float(np.arccos(np.clip(cos_nu2, -1.0, 1.0)))
            tau = elliptic_time_of_flight(nu1, nu2, a, e, gm)
            if not tau > 0:
                raise RuntimeError("non-positive time of flight")
            arrival = propagate_kepler(state, tau, gm).final_state
            radius = float(np.linalg.norm(arrival[:3]))
            target = np.sqrt(gm / r2) * np.cross([0.0, 0.0, 1.0], arrival[:3] / radius)
            burn2 = target - arrival[3:]
            objective = total_delta_v([burn1, burn2])
        except ValueError as exc:
            self.captured_value_errors += 1
            return self._invalid(point, EvaluationStatus.INVALID_DOMAIN,
                                 f"captured ValueError: {exc}", captured_exception="ValueError")
        except RuntimeError as exc:
            self.captured_runtime_errors += 1
            return self._invalid(point, EvaluationStatus.NUMERICAL_FAILURE,
                                 f"captured RuntimeError: {exc}", captured_exception="RuntimeError")
        metadata = {
            "delta_v_total_m_s": objective,
            "delta_v1_m_s": float(np.linalg.norm(burn1)),
            "delta_v2_m_s": float(np.linalg.norm(burn2)),
            "first_burn_m_s": burn1,
            "second_burn_m_s": burn2,
            "transfer_time_s": tau,
            "semi_major_axis_m": a,
            "eccentricity": e,
            "periapsis_radius_m": a * (1 - e),
            "apoapsis_radius_m": a * (1 + e),
            "departure_true_anomaly_rad": nu1,
            "arrival_true_anomaly_rad": nu2,
            "arrival_radius_residual_m": radius - r2,
            "tangency_clip": float(max(0.0, -1.0 - cos_nu2)),
            "reach_margin_relative": q / s,
            "reference_delta_v_m_s": self.reference.total_delta_v,
            "reference_transfer_time_s": self.reference.transfer_time,
        }
        numbers = [value for value in metadata.values() if not isinstance(value, np.ndarray)]
        if not (np.all(np.isfinite(numbers)) and np.all(np.isfinite(burn2))):
            self.captured_runtime_errors += 1
            return self._invalid(point, EvaluationStatus.NUMERICAL_FAILURE, "non-finite transfer quantities")
        return Evaluation(point, objective, EvaluationStatus.VALID, "", metadata)

    def plan(self, x: ArrayLike) -> TwoImpulsePlan:
        """TwoImpulsePlan of a VALID decision vector; an invalid one raises ValueError."""
        evaluation = self.evaluate(x)
        if not evaluation.is_valid:
            raise ValueError(f"decision vector is not valid: {evaluation.message}")
        meta = evaluation.metadata
        return TwoImpulsePlan(self._initial_state, meta["first_burn_m_s"], meta["second_burn_m_s"],
                              meta["transfer_time_s"], self.final_radius, self.mu)
