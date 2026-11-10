"""Independent reference solutions for constrained two-impulse transfers.

The Hohmann transfer is the reference of the unconstrained problem. With a
limit on the transfer time shorter than the Hohmann time there is no closed
form for the minimum, so the reference is computed here by a method that
shares nothing with the optimizers or with the evaluation chain of the
problems (state_to_elements, time of flight, propagate_kepler):

  1. ``two_impulse_closed_form`` evaluates cost and time of an ascending
     two-impulse transfer directly from vis-viva, the conservation of
     angular momentum and Kepler's equation in eccentric anomaly;
  2. ``time_limited_reference`` minimizes that cost along the boundary of
     the feasible set and cross-checks the result against an exhaustive
     two-dimensional grid.

Closed form, for a departure from the circular orbit r1 with first impulse
dv at angle alpha from the velocity (positive outward), r2 > r1:

    v_t = v_c1 + dv cos(alpha),  v_r = dv sin(alpha),  h = r1 v_t
    eps = (v_t^2 + v_r^2)/2 - mu/r1,  a = -mu/(2 eps),  e^2 = 1 - h^2/(mu a)
    speed at r2:    v2^2 = 2 (eps + mu/r2),  transverse h/r2,
                    radial sqrt(v2^2 - (h/r2)^2)
    second impulse: dv2 = hypot(v_c2 - h/r2, radial)        (circularization)
    e sin(E1) = r1 v_r / sqrt(mu a),  e cos(E1) = 1 - r1/a
    cos(E2) = (1 - r2/a)/e, E2 in [0, pi]                   (outward crossing)
    tau = ((E2 - e sin E2) - (E1 - e sin E1)) mod 2 pi / sqrt(mu/a^3)
    minimum radius: a(1 - e) if v_r < 0 (the vehicle first descends), else r1

Feasible set for a direction alpha. The transfer reaches r2 iff
dv >= dv_reach(alpha), the positive root of the quadratic obtained from
Q = 2 eps r2^2 + 2 mu r2 - h^2 = 0:

    A dv^2 + B dv + C = 0,   A = r2^2 - r1^2 cos^2(alpha),
    B = 2 v_c1 cos(alpha) (r2^2 - r1^2),   C = -mu (r2 - r1)^2 / r1.

Above it the transfer time decreases as dv grows (a faster departure arrives
sooner), so the candidates with time <= t_max are dv >= dv_min(alpha), where
dv_min is dv_reach if the tangent arrival already meets the limit and the
root of tau(dv, alpha) = t_max otherwise. Since the cost grows with dv along
that ray in the region of interest, the reference is

    min over alpha of J(dv_min(alpha), alpha).

Both assumptions (time decreasing and cost increasing with dv) are *checked*,
not trusted: an exhaustive grid over the whole decision box is evaluated as
well, and a feasible node cheaper than the boundary minimum raises
RuntimeError instead of returning a wrong reference. The returned point is
feasible by construction, so it is an upper bound of the true minimum; the
grid check bounds how much lower the minimum could be by the variation of
the cost over one grid cell.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..constants import EARTH_MU
from ..dynamics.two_body import real_array, validate_mu
from ..orbital.transfers import hohmann_transfer


def two_impulse_closed_form(initial_radius: float, final_radius: float, mu: float,
                            delta_v: ArrayLike, alpha: ArrayLike) -> dict[str, NDArray[np.float64]]:
    """Cost, time and minimum radius of ascending two-impulse transfers, vectorized.

    ``delta_v`` (m/s) and ``alpha`` (rad) broadcast together. Returns arrays
    ``valid`` (bound, prograde and reaching r2), ``delta_v_total_m_s``,
    ``delta_v2_m_s``, ``transfer_time_s`` and ``minimum_radius_m``; entries
    that are not valid hold NaN. Requires 0 < r1 < r2. See the module
    docstring for the equations.
    """
    r1, r2, gm = float(initial_radius), float(final_radius), validate_mu(mu)
    if not 0 < r1 < r2:
        raise ValueError("the closed form requires 0 < initial_radius < final_radius")
    dv, angle = np.broadcast_arrays(np.asarray(delta_v, dtype=float), np.asarray(alpha, dtype=float))
    circular_1, circular_2 = np.sqrt(gm / r1), np.sqrt(gm / r2)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        tangential = circular_1 + dv * np.cos(angle)
        radial = dv * np.sin(angle)
        energy = (tangential**2 + radial**2) / 2 - gm / r1
        momentum = r1 * tangential
        a = -gm / (2 * energy)
        e = np.sqrt(np.maximum(0.0, 1 - momentum**2 / (gm * a)))
        valid = (energy < 0) & (tangential > 0) & (a * (1 + e) >= r2) & (e > 0)
        speed2 = 2 * (energy + gm / r2)
        transverse_2 = momentum / r2
        radial_2 = np.sqrt(np.maximum(0.0, speed2 - transverse_2**2))
        second = np.hypot(circular_2 - transverse_2, radial_2)
        start = np.arctan2(r1 * radial / np.sqrt(gm * a), 1 - r1 / a)
        arrival = np.arccos(np.clip((1 - r2 / a) / e, -1.0, 1.0))
        sweep = np.mod((arrival - e * np.sin(arrival)) - (start - e * np.sin(start)), 2 * np.pi)
        time = sweep * np.sqrt(a / gm) * a
        lowest = np.where(radial < 0, a * (1 - e), r1)
    blank = np.where(valid, 1.0, np.nan)
    return {
        "valid": valid,
        "delta_v_total_m_s": (dv + second) * blank,
        "delta_v2_m_s": second * blank,
        "transfer_time_s": time * blank,
        "minimum_radius_m": lowest * blank,
    }


@dataclass(frozen=True)
class TransferReference:
    """Reference solution of a two-impulse transfer under a time limit, SI units.

    delta_v_total_m_s, delta_v1_m_s, alpha1_rad, transfer_time_s and
    minimum_radius_m describe the reference point, which is feasible.
    active: whether the time limit binds (False: the Hohmann transfer is the
    solution). check_grid_best_m_s: cost of the cheapest feasible node of the
    exhaustive cross-check grid (NaN when the limit does not bind); it is
    never below the reference by more than the check tolerance. method: how
    the reference was obtained.
    """

    delta_v_total_m_s: float
    delta_v1_m_s: float
    alpha1_rad: float
    transfer_time_s: float
    minimum_radius_m: float
    active: bool
    check_grid_best_m_s: float
    method: str

    def to_dict(self) -> dict[str, Any]:
        return dict(self.__dict__)


# A feasible grid node cheaper than the boundary minimum by more than this (m/s) is an error.
CHECK_TOLERANCE_M_S: float = 1e-6


def reach_delta_v(initial_radius: float, final_radius: float, mu: float, alpha: ArrayLike) -> NDArray[np.float64]:
    """Smallest first impulse (m/s) at angle alpha (rad) whose transfer orbit reaches r2 > r1.

    Positive root of A dv^2 + B dv + C = 0 (module docstring). Vectorized in alpha.
    """
    r1, r2 = float(initial_radius), float(final_radius)
    cosine = np.cos(np.asarray(alpha, dtype=float))
    circular = np.sqrt(mu / r1)
    a = r2 * r2 - r1 * r1 * cosine * cosine
    b = 2 * circular * cosine * (r2 * r2 - r1 * r1)
    c = -mu * (r2 - r1) ** 2 / r1
    return (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)


def _boundary_cost(r1: float, r2: float, mu: float, limit: float, upper: float,
                   alpha: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """(dv_min, cost) for every direction: the cheapest feasible first impulse and its total cost.

    Entries are (NaN, +inf) for directions in which no first impulse below
    ``upper`` meets the time limit.
    """
    def time_of(speed: NDArray[np.float64]) -> NDArray[np.float64]:
        return two_impulse_closed_form(r1, r2, mu, speed, alpha)["transfer_time_s"]

    # Just above tangency and just below the bound, so that both ends are inside the valid domain.
    low = reach_delta_v(r1, r2, mu, alpha) * (1 + 1e-12)
    high = np.full_like(alpha, upper * (1 - 1e-12))
    with np.errstate(invalid="ignore"):
        tangent_meets_limit = (low < high) & (time_of(low) <= limit)
        needs_root = (low < high) & ~tangent_meets_limit & (time_of(np.where(low < high, low, high)) > limit) \
            & (time_of(high) <= limit)
    left, right = low.copy(), high.copy()
    for _ in range(100):                       # bisection keeps `right` on the feasible side
        middle = 0.5 * (left + right)
        with np.errstate(invalid="ignore"):
            meets = time_of(middle) <= limit
        right = np.where(needs_root & meets, middle, right)
        left = np.where(needs_root & ~meets, middle, left)
    speed = np.where(tangent_meets_limit, low, np.where(needs_root, right, np.nan))
    cost = two_impulse_closed_form(r1, r2, mu, speed, alpha)["delta_v_total_m_s"]
    return speed, np.where(np.isfinite(cost), cost, np.inf)


def time_limited_reference(initial_radius: float, final_radius: float, max_duration: float,
                           mu: float = EARTH_MU, *, directions: int = 4001, check_nodes: int = 401) -> TransferReference:
    """Minimum total delta-v of an ascending two-impulse transfer with transfer time <= max_duration.

    Radii in metres with 0 < r1 < r2, ``max_duration`` > 0 in seconds, mu in
    m^3/s^2; the first impulse is below (sqrt(2)-1) v_c1, the bound of the
    transfer problems. If the limit is not shorter than the Hohmann time the
    Hohmann transfer is returned. Otherwise (module docstring):

      1. for ``directions`` angles in (-pi, pi), the cheapest feasible first
         impulse is found (tangent arrival, or bisection of tau = t_max);
      2. the best direction is refined by a golden-section search between its
         two neighbours;
      3. an exhaustive ``check_nodes`` x ``check_nodes`` grid over the whole
         decision box is evaluated as a cross-check.

    No feasible direction raises ValueError (the limit cannot be met within
    the bounds). A feasible check node cheaper than the result by more than
    CHECK_TOLERANCE_M_S raises RuntimeError.
    """
    r1, r2, gm = float(initial_radius), float(final_radius), validate_mu(mu)
    limit = real_array(max_duration, "max_duration")
    if limit.shape != () or limit <= 0:
        raise ValueError("max_duration must be a positive scalar in seconds")
    limit = float(limit)
    if not 0 < r1 < r2:
        raise ValueError("the reference requires 0 < initial_radius < final_radius")
    for value, name in ((directions, "directions"), (check_nodes, "check_nodes")):
        if isinstance(value, bool) or not isinstance(value, int) or value < 5:
            raise ValueError(f"{name} must be an integer >= 5")
    hohmann = hohmann_transfer(r1, r2, gm)
    if limit >= hohmann.transfer_time:
        return TransferReference(hohmann.total_delta_v, hohmann.first_delta_v, 0.0, hohmann.transfer_time,
                                 r1, False, float("nan"), "Hohmann closed form: the time limit does not bind")
    upper = float((np.sqrt(2.0) - 1.0) * np.sqrt(gm / r1))

    angles = np.linspace(-np.pi, np.pi, directions + 2)[1:-1]
    speeds, costs = _boundary_cost(r1, r2, gm, limit, upper, angles)
    index = int(np.argmin(costs))
    if not np.isfinite(costs[index]):
        raise ValueError("no two-impulse transfer meets the time limit within the decision bounds")
    best = (float(costs[index]), float(speeds[index]), float(angles[index]))

    def at(angle: float) -> tuple[float, float, float]:
        speed, cost = _boundary_cost(r1, r2, gm, limit, upper, np.array([angle]))
        return float(cost[0]), float(speed[0]), angle

    step = angles[1] - angles[0]
    left, right = best[2] - step, best[2] + step
    ratio = (np.sqrt(5.0) - 1.0) / 2.0
    inner_left, inner_right = at(right - ratio * (right - left)), at(left + ratio * (right - left))
    for _ in range(80):
        if inner_left[0] <= inner_right[0]:
            right, inner_right = inner_right[2], inner_left
            inner_left = at(right - ratio * (right - left))
        else:
            left, inner_left = inner_left[2], inner_right
            inner_right = at(left + ratio * (right - left))
        best = min(best, inner_left, inner_right)

    grid_speeds = np.linspace(0.0, upper, check_nodes + 1)[1:-1]
    grid_angles = np.linspace(-np.pi, np.pi, check_nodes)
    grid = two_impulse_closed_form(r1, r2, gm, grid_speeds[:, None], grid_angles[None, :])
    with np.errstate(invalid="ignore"):
        feasible = grid["valid"] & (grid["transfer_time_s"] <= limit)
    grid_best = float(np.min(np.where(feasible, grid["delta_v_total_m_s"], np.inf)))
    if grid_best < best[0] - CHECK_TOLERANCE_M_S:
        raise RuntimeError(f"cross-check failed: a feasible grid node costs {grid_best:.9f} m/s, "
                           f"below the boundary minimum {best[0]:.9f} m/s")
    final = two_impulse_closed_form(r1, r2, gm, best[1], best[2])
    return TransferReference(
        delta_v_total_m_s=best[0], delta_v1_m_s=best[1], alpha1_rad=best[2],
        transfer_time_s=float(final["transfer_time_s"]), minimum_radius_m=float(final["minimum_radius_m"]),
        active=True, check_grid_best_m_s=grid_best,
        method=f"closed-form transfer; minimum over {directions} directions of the cheapest first impulse "
               f"meeting the limit, golden-section refinement, {check_nodes}x{check_nodes} exhaustive check grid",
    )
