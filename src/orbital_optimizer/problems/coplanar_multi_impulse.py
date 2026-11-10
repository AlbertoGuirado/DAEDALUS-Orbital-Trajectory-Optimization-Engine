"""Coplanar circular-to-circular transfer with N impulses as an optimization Problem.

Adapter layer. It generalizes problems/coplanar_transfer.py (ascending, two
impulses), which stays as the validated Phase 3 reference: for N = 2 and
r2 > r1 both evaluators define the same objective, and a test compares them.
Here the transfer may ascend or descend and may use N >= 2 impulses, so the
instants of the impulses are part of the decision vector (RF-3.05, RF-3.07).

Model
-----
Two-body point mass with constant mu; concentric coplanar circular orbits of
radii r1 != r2 in the xy plane, prograde; instantaneous impulses; free
transfer time; the target is the orbit, not a point on it. By the rotational
symmetry of the initial circle the departure is fixed at r0 = (r1, 0, 0),
v0 = (0, v_c1, 0), t1 = 0, with v_c = sqrt(mu/r).

Local frame at a burn: r_hat (radial, outward) and theta_hat = z_hat x r_hat
(transverse, prograde). A free burn k of magnitude dv_k and angle alpha_k is

    dv_k_vec = dv_k (sin(alpha_k) r_hat + cos(alpha_k) theta_hat),

so alpha = 0 is prograde, alpha = +-pi retrograde and 0 < alpha < pi points outward.
At the departure point theta_hat is the velocity direction, which is the
convention of the two-impulse reference problem.

Decision vector (SI), D = 2 + 3 (N - 2):

    x = (dv_1, alpha_1, [M_k, dv_k, alpha_k  for k = 2 .. N-1])
    dv_k     in [0, dv_max]   m/s
    alpha_k  in [-pi, pi) rad when ascending, [0, 2 pi) rad when descending;
             periodic in both cases
    M_k      in [0, 2 pi)     rad, periodic: mean-anomaly sweep of the coast
             between burns k-1 and k, on the orbit left by burn k-1.

The two angle intervals describe the same directions; only the branch cut
moves. It is placed opposite to the burn direction of the Hohmann-like
solutions (prograde when ascending, retrograde when descending) so that the
expected optimum is not on the seam of the periodic variable. Reason: the
mutation of DE/rand/1/bin@1 wraps periodic variables but does not take
shortest-arc differences, and with the cut at the optimum (alpha = +-pi for a
descending Hohmann transfer) its final gaps over five seeds at 2000
evaluations were 0.04 to 12.9 m/s (GEO->LEO, 2026-11-10). A descending
bi-elliptic transfer, whose first burn is prograde, still has that burn on
the seam; this only matters for r1/r2 > 11.94, where it is the better strategy.

The impulse instants follow from the decision: t_1 = 0 and
t_k = t_{k-1} + M_k / n_{k-1}, n = sqrt(mu/a^3) of the orbit being coasted.
The last impulse N is not a decision variable: at the first crossing of r2
after burn N-1 (outward if the vehicle is below r2, inward if above) it is
the unique burn that circularizes there, dv_N_vec = v_c2 theta_hat - v. This
removes the equality constraint "arrive on the target orbit" by construction.

dv_max = (sqrt(2) - 1) sqrt(mu / min(r1, r2)): escape excess of the lower
orbit. For an ascending transfer it is the bound of the reference problem
(every first burn stays bound); it exceeds every Hohmann and bi-elliptic burn
in both directions. It does not by itself keep later burns bound, so each
post-burn state is checked.

Objective (m/s): J(x) = sum_k |dv_k_vec| (orbital.mechanics.total_delta_v).
Analytic references: Hohmann (N = 2) and bi-elliptic (N = 3) from
orbital.transfers. No claim of optimality is made for this problem class
beyond what those closed forms state.

Evaluation chain
----------------
For each free burn: (coast M_k analytically with propagate_kepler when
k >= 2) -> add the burn -> predicates on the new state: bound and not
numerically parabolic, not radial or retrograde (thresholds of
state_to_elements). Then the reach test of the reference problem for r2
(``reach_margin`` with the tangency tolerance K = 32), the crossing true
anomaly, elliptic_time_of_flight and propagate_kepler to the crossing, and
the circularizing burn. Candidates failing a predicate are INVALID_DOMAIN. A
zero coast between two burns (coincident impulses) is INVALID_DOMAIN: the
impulse instants must be strictly increasing (DT-13). ValueError or
RuntimeError raised by the physics calls are counted and reported as
INVALID_DOMAIN or NUMERICAL_FAILURE; other exceptions propagate.

Tangency tolerance after coast arcs
-----------------------------------
The reference problem accepts Q >= -K u S with K = 32 because its only burn
is applied to an exactly known departure state. Here a free burn can follow
coast arcs, and the bi-elliptic transfer -- the analytic reference for
N = 3 -- is tangent to the target after such an arc, so the propagated state
carries the roundoff of the analytic propagation into Q. Error model of
simulation.kepler_propagator: the semi-major axis recomputed from a
propagated state has relative error up to 4 eps (3+e)/(1-e) = 8 u kappa,
kappa = (3+e)/(1-e) of the arc. In Q = 2 eps_orb r2^2 + 2 mu r2 - h^2:
  - energy term: |2 r2^2 d(eps_orb)| = 2 r2^2 |eps_orb| (da/a) <= 8 u kappa S;
  - angular momentum term: the f and g coefficients lose the same factor to
    cancellation (g_dot = 1 - 2/(1+e) at apoapsis), so |2 h dh| <= 8 u kappa S.
Each preceding coast arc therefore adds 16 kappa to the factor:

    accept iff Q_computed >= -(32 + 16 sum_arcs kappa_arc) u S.

The constant 16 follows from this argument, not from a fit. Confirmation
(2026-11-10, bi-elliptic decisions with r2/r1 from 0.05 to 20 and apoapsis up
to 12 000 r1, kappa up to 2.4e4): |Q| / (u S kappa) was at most 0.61, and Q
was negative in 4 of 7 cases, which a tolerance without the arc term rejected.
The admitted infeasible band grows accordingly: a relative deficit of about
16 kappa u in the apsis radius (2e-14 for a coast arc with e = 0.7, 4e-11 for
e = 0.9998). ``tangency_factor`` in the metadata records the factor used.

Conditioning at tangency: the crossing anomaly is arccos(c) with c = +-1
there, so an error d in c moves it by sqrt(2 d). The transfer time of a
tangent arrival therefore carries an error of up to about 4 sqrt(F u) radians
of true anomaly at the apsis rate (F the tangency factor): of order 1e-3 s
for the transfers tested, against 1e-6 s for a single-arc Hohmann transfer.
The objective is not affected to first order, because the circularizing burn
is stationary with respect to the arrival point at an apsis.

Trajectory extrema: ``minimum_radius_m`` and ``maximum_radius_m`` are the
exact extrema of |r| over the whole transfer, from
orbital.kepler.radius_extrema_on_arc on every coast arc (not from samples).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..dynamics.two_body import EARTH_MU, real_array, validate_mu, validate_state
from ..optimization.problem import Bounds, Evaluation, EvaluationStatus
from ..orbital.elements import SINGULAR_TOL, state_to_elements
from ..orbital.kepler import elliptic_time_of_flight, radius_extrema_on_arc
from ..orbital.mechanics import orbital_period, total_delta_v
from ..orbital.transfers import bielliptic_transfer
from ..simulation.impulsive import propagate_with_impulses
from ..simulation.kepler_propagator import propagate_kepler
from ..simulation.propagator import propagate
from .coplanar_transfer import (
    TANGENCY_FACTOR, UNIT_ROUNDOFF, VERIFICATION_ECCENTRICITY_TOLERANCE, VERIFICATION_RADIUS_TOLERANCE_M,
    PlanVerification, _hohmann_reference, _radius, reach_margin,
)


# Added to TANGENCY_FACTOR once per preceding coast arc, times (3+e)/(1-e) of that arc.
ARC_TANGENCY_FACTOR: float = 16.0


@dataclass(frozen=True, eq=False)
class ImpulsivePlan:
    """Oriented N-impulse plan in the inertial frame, SI units.

    initial_state (6,) [m, m/s] at t = 0; burn_times (N,) s, strictly
    increasing and >= 0; burn_vectors (N, 3) m/s, burn k applied at
    burn_times[k]; target circular orbit of final_radius (m) about mu
    (m^3/s^2). N >= 1. Arrays are read-only copies.
    """

    initial_state: NDArray[np.float64]
    burn_times: NDArray[np.float64]
    burn_vectors: NDArray[np.float64]
    final_radius: float
    mu: float = EARTH_MU

    def __post_init__(self) -> None:
        state = validate_state(self.initial_state)
        state.setflags(write=False)
        times = real_array(self.burn_times, "burn_times")
        vectors = real_array(self.burn_vectors, "burn_vectors")
        if times.ndim != 1 or times.size < 1:
            raise ValueError("burn_times must be a nonempty one-dimensional array in seconds")
        if times[0] < 0 or np.any(np.diff(times) <= 0):
            raise ValueError("burn_times must be nonnegative and strictly increasing")
        if vectors.shape != (times.size, 3):
            raise ValueError("burn_vectors must have shape (N, 3) in m/s")
        times.setflags(write=False)
        vectors.setflags(write=False)
        object.__setattr__(self, "initial_state", state)
        object.__setattr__(self, "burn_times", times)
        object.__setattr__(self, "burn_vectors", vectors)
        object.__setattr__(self, "final_radius", _radius(self.final_radius, "final_radius"))
        object.__setattr__(self, "mu", validate_mu(self.mu))

    @property
    def total_delta_v(self) -> float:
        """sum_k |burn_k|, m/s."""
        return total_delta_v(self.burn_vectors)

    @property
    def transfer_time(self) -> float:
        """Instant of the last burn, s."""
        return float(self.burn_times[-1])

    def impulses(self) -> list[tuple[float, NDArray[np.float64]]]:
        """[(t_k, burn_k)] for simulation.impulsive.propagate_with_impulses."""
        return [(float(time), vector.copy()) for time, vector in zip(self.burn_times, self.burn_vectors)]


def verify_impulsive_plan(
    plan: ImpulsivePlan,
    *,
    radius_tolerance_m: float = VERIFICATION_RADIUS_TOLERANCE_M,
    eccentricity_tolerance: float = VERIFICATION_ECCENTRICITY_TOLERANCE,
    method: Literal["RK45", "RK4"] = "RK45",
    max_step: float = 60.0,
    rtol: float = 1e-10,
    atol: ArrayLike | None = None,
) -> PlanVerification:
    """Independently propagate a plan numerically and check arrival on the target circle (RV-3.04).

    1. propagate_with_impulses(initial_state, impulses, t_N, mu) with the
       given integrator options (defaults: those of ``propagate``).
    2. Explicit finiteness and state_to_elements on the final state.
    3. |a - r2| <= radius_tolerance_m, e < eccentricity_tolerance,
       | |r| - r2 | <= radius_tolerance_m. Defaults 1 m and 1e-6: the Phase 2
       arrival criteria.
    4. Coast of one target period with ``propagate``: the same radius bound
       over the accepted steps and the same eccentricity bound at the end.
    A plan whose last burn is at t = 0 has no transfer arc and fails with a
    message. Propagation or element errors make the verification fail with a
    message instead of raising. Local integrator tolerances do not bound the
    global error: for transfers much longer than LEO-GEO the caller must
    tighten ``rtol``/``atol`` rather than loosen the arrival tolerances.
    """
    if not isinstance(plan, ImpulsivePlan):
        raise ValueError("plan must be an ImpulsivePlan")
    limits = real_array([radius_tolerance_m, eccentricity_tolerance], "tolerances")
    if np.any(limits <= 0):
        raise ValueError("verification tolerances must be positive")
    radius_tol, ecc_tol = float(limits[0]), float(limits[1])
    r2, gm = plan.final_radius, plan.mu
    options = {"method": method, "max_step": max_step, "rtol": rtol, "atol": atol}
    tolerances = {"radius_m": radius_tol, "eccentricity": ecc_tol}
    checks: dict[str, bool] = {}
    metrics: dict[str, float] = {}

    def failed(message: str, final: NDArray | None = None) -> PlanVerification:
        return PlanVerification(False, checks, metrics, final, message, tolerances)

    if plan.transfer_time <= 0:
        return failed("the plan has no transfer arc: its last burn is at t = 0")
    try:
        final = propagate_with_impulses(plan.initial_state, plan.impulses(), plan.transfer_time, gm, **options)
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
        coast = propagate(final, orbital_period(r2, gm), gm, **options)
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


@dataclass(frozen=True, eq=False)
class _Burn:
    time: float
    vector: NDArray[np.float64]


@dataclass(frozen=True, eq=False)
class _Failure:
    status: EvaluationStatus
    message: str
    metadata: Mapping[str, Any] = field(default_factory=dict)


def _local_frame(state: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """Radial and transverse unit vectors (xy plane, prograde) and the radius, m."""
    radius = float(np.hypot(state[0], state[1]))
    radial = np.array([state[0] / radius, state[1] / radius, 0.0])
    return radial, np.array([-radial[1], radial[0], 0.0]), radius


class CoplanarMultiImpulseProblem:
    """Coplanar circular-to-circular transfer with N >= 2 impulses (see module docstring).

    Decision x = (dv_1, alpha_1, [M_k, dv_k, alpha_k]); objective J in m/s.
    r1 != r2, both > 0, in metres; mu in m^3/s^2. Equal radii raise
    ValueError; a non-finite reference or bound raises RuntimeError. The
    counters ``captured_value_errors`` and ``captured_runtime_errors``
    accumulate over the lifetime of the instance.
    """

    def __init__(self, initial_radius: float, final_radius: float, mu: float = EARTH_MU, *,
                 impulses: int = 2) -> None:
        r1, r2 = _radius(initial_radius, "initial_radius"), _radius(final_radius, "final_radius")
        gm = validate_mu(mu)
        if r1 == r2:
            raise ValueError("initial and final radii must differ")
        if isinstance(impulses, bool) or not isinstance(impulses, (int, np.integer)) or impulses < 2:
            raise ValueError("impulses must be an integer >= 2")
        reference = _hohmann_reference(r1, r2, gm)
        with np.errstate(over="ignore", invalid="ignore"):
            circular_speed = float(np.sqrt(gm / r1))
            upper = float((np.sqrt(2.0) - 1.0) * np.sqrt(gm / min(r1, r2)))
        values = (reference.first_delta_v, reference.second_delta_v, reference.total_delta_v,
                  reference.transfer_time, circular_speed, upper)
        if not all(np.isfinite(values)) or upper <= 0:
            raise RuntimeError("Hohmann reference or decision bounds are not finite at this numerical scale")
        self.impulses = int(impulses)
        self.name = f"coplanar-{self.impulses}-impulse"
        self.initial_radius, self.final_radius, self.mu = r1, r2, gm
        self.ascending = r2 > r1
        self.reference = reference
        angle_low = -np.pi if self.ascending else 0.0
        angle_high = angle_low + 2 * np.pi
        lower, high = [0.0, angle_low], [upper, angle_high]
        names, units, periodic = ["delta_v1", "alpha1"], ["m/s", "rad"], [False, True]
        for k in range(2, self.impulses):
            lower += [0.0, 0.0, angle_low]
            high += [2 * np.pi, upper, angle_high]
            names += [f"coast{k}", f"delta_v{k}", f"alpha{k}"]
            units += ["rad", "m/s", "rad"]
            periodic += [True, False, True]
        self.bounds = Bounds(lower, high, names=tuple(names), units=tuple(units), periodic=np.array(periodic))
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
        """Decision of the Hohmann transfer; defined for N = 2 only (else ValueError).

        (dv1_H, 0) for an ascending transfer and (dv1_H, pi) for a descending
        one (retrograde first burn).
        """
        if self.impulses != 2:
            raise ValueError("the Hohmann decision is defined for two impulses")
        return np.array([self.reference.first_delta_v, 0.0 if self.ascending else np.pi])

    def bielliptic_decision(self, intermediate_radius: float) -> NDArray[np.float64]:
        """Decision of the bi-elliptic transfer through rb; defined for N = 3 only (else ValueError).

        (dv1, 0, pi, dv2, alpha2) with the magnitudes of bielliptic_transfer
        and a half-revolution coast to the apoapsis rb. The first burn raises
        the apoapsis to rb, so it is prograde in both directions; the second
        moves the periapsis from r1 to r2: prograde (alpha2 = 0) when
        r2 > r1, retrograde (alpha2 = pi) otherwise.
        """
        if self.impulses != 3:
            raise ValueError("the bi-elliptic decision is defined for three impulses")
        transfer = bielliptic_transfer(self.initial_radius, self.final_radius, intermediate_radius, self.mu)
        return np.array([transfer.first_delta_v, 0.0, np.pi, transfer.second_delta_v,
                         0.0 if self.ascending else np.pi])

    def _trajectory(self, point: NDArray[np.float64]) -> _Failure | dict[str, Any]:
        """Build the transfer of a validated decision vector; return its record or a failure."""
        r2, gm = self.final_radius, self.mu
        state = self._initial_state.copy()
        clock = 0.0
        burns: list[_Burn] = []
        low = high = self.initial_radius
        tangency_factor = TANGENCY_FACTOR
        index = 0
        for k in range(1, self.impulses):
            if k > 1:
                sweep = float(point[index])
                index += 1
                elements = state_to_elements(state, gm)
                a, e = elements.semi_major_axis, elements.eccentricity
                coast = float(sweep * np.sqrt(a / gm) * a)
                if not coast > 0:
                    return _Failure(EvaluationStatus.INVALID_DOMAIN,
                                    "coincident impulses: the coast between two burns is zero")
                arc_low, arc_high = radius_extrema_on_arc(elements.true_anomaly, coast, a, e, gm)
                low, high = min(low, arc_low), max(high, arc_high)
                state = propagate_kepler(state, coast, gm).final_state
                clock += coast
                tangency_factor += ARC_TANGENCY_FACTOR * (3 + e) / (1 - e)
            magnitude, alpha = float(point[index]), float(point[index + 1])
            index += 2
            radial_hat, transverse_hat, radius = _local_frame(state)
            burn = magnitude * (np.sin(alpha) * radial_hat + np.cos(alpha) * transverse_hat)
            state[3:] = state[3:] + burn
            burns.append(_Burn(clock, burn))
            radial = float(np.dot(state[3:], radial_hat))
            tangential = float(np.dot(state[3:], transverse_hat))
            speed2 = tangential * tangential + radial * radial
            energy = speed2 / 2 - gm / radius
            if not energy < -SINGULAR_TOL * gm / radius:
                return _Failure(EvaluationStatus.INVALID_DOMAIN,
                                f"orbit after burn {k} is unbound or numerically parabolic",
                                {"specific_energy_j_kg": energy})
            if not tangential > SINGULAR_TOL * np.sqrt(speed2):
                return _Failure(EvaluationStatus.INVALID_DOMAIN,
                                f"orbit after burn {k} is radial or retrograde",
                                {"tangential_speed_m_s": tangential})

        q, s = reach_margin(radius, r2, gm, tangential, radial)
        if q < -tangency_factor * UNIT_ROUNDOFF * s:
            return _Failure(EvaluationStatus.INVALID_DOMAIN, "transfer orbit does not reach the target radius",
                            {"reach_margin_relative": q / s, "tangency_factor": tangency_factor})
        if radius == r2:
            return _Failure(EvaluationStatus.INVALID_DOMAIN,
                            "the last free burn is on the target radius: the crossing is undefined")
        elements = state_to_elements(state, gm)
        a, e, nu = elements.semi_major_axis, elements.eccentricity, elements.true_anomaly
        if not e > 0:
            raise RuntimeError("transfer orbit is numerically circular and cannot cross the target radius")
        cos_crossing = (a * (1 - e) * (1 + e) / r2 - 1) / e
        outward = float(np.arccos(np.clip(cos_crossing, -1.0, 1.0)))
        crossing = outward if radius < r2 else 2 * np.pi - outward
        tau = elliptic_time_of_flight(nu, crossing, a, e, gm)
        if not tau > 0:
            raise RuntimeError("non-positive time of flight to the target radius")
        arc_low, arc_high = radius_extrema_on_arc(nu, tau, a, e, gm)
        low, high = min(low, arc_low), max(high, arc_high)
        arrival = propagate_kepler(state, tau, gm).final_state
        radial_hat, transverse_hat, arrival_radius = _local_frame(arrival)
        final_burn = np.sqrt(gm / r2) * transverse_hat - arrival[3:]
        burns.append(_Burn(clock + tau, final_burn))
        vectors = np.array([burn.vector for burn in burns])
        magnitudes = np.hypot.reduce(vectors, axis=1)
        return {
            "delta_v_total_m_s": total_delta_v(vectors),
            "burn_delta_v_m_s": magnitudes,
            "burn_times_s": np.array([burn.time for burn in burns]),
            "burn_vectors_m_s": vectors,
            "transfer_time_s": clock + tau,
            "minimum_radius_m": low,
            "maximum_radius_m": high,
            "impulse_kinetic_energy_j_kg": float(np.sum(magnitudes**2) / 2),
            "closest_approach_index": min(self.initial_radius, r2) / low,
            "final_arc_semi_major_axis_m": a,
            "final_arc_eccentricity": e,
            "arrival_radius_residual_m": arrival_radius - r2,
            "tangency_clip": float(max(0.0, abs(cos_crossing) - 1.0)),
            "reach_margin_relative": q / s,
            "tangency_factor": tangency_factor,
            "reference_delta_v_m_s": self.reference.total_delta_v,
            "reference_transfer_time_s": self.reference.transfer_time,
            "impulses": self.impulses,
        }

    def evaluate(self, x: ArrayLike) -> Evaluation:
        """Evaluate J(x) in m/s for x inside the bounds (else ValueError)."""
        point = self.bounds.validate(x)
        try:
            outcome = self._trajectory(point)
        except ValueError as exc:
            self.captured_value_errors += 1
            return Evaluation(point, np.nan, EvaluationStatus.INVALID_DOMAIN,
                              f"captured ValueError: {exc}", {"captured_exception": "ValueError"})
        except RuntimeError as exc:
            self.captured_runtime_errors += 1
            return Evaluation(point, np.nan, EvaluationStatus.NUMERICAL_FAILURE,
                              f"captured RuntimeError: {exc}", {"captured_exception": "RuntimeError"})
        if isinstance(outcome, _Failure):
            return Evaluation(point, np.nan, outcome.status, outcome.message, outcome.metadata)
        finite = all(np.all(np.isfinite(value)) for value in outcome.values())
        if not finite:
            self.captured_runtime_errors += 1
            return Evaluation(point, np.nan, EvaluationStatus.NUMERICAL_FAILURE, "non-finite transfer quantities")
        return Evaluation(point, outcome["delta_v_total_m_s"], EvaluationStatus.VALID, "", outcome)

    def plan(self, x: ArrayLike) -> ImpulsivePlan:
        """ImpulsivePlan of a VALID decision vector; an invalid one raises ValueError."""
        evaluation = self.evaluate(x)
        if not evaluation.is_valid:
            raise ValueError(f"decision vector is not valid: {evaluation.message}")
        meta = evaluation.metadata
        return ImpulsivePlan(self._initial_state, meta["burn_times_s"], meta["burn_vectors_m_s"],
                             self.final_radius, self.mu)
