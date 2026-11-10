"""Elliptic Kepler equation in eccentric anomaly; radians, dimensionless e.

Mean anomaly M=E-e*sin(E). The propagation form uses the anomaly change
x=E-E0 from a state with c=e*cos(E0), s=e*sin(E0):

    tau = x - c*sin(x) + s*(1-cos(x)),  dtau/dx = 1 - c*cos(x) + s*sin(x) = r/a.

The derivative lies in [1-e, 1+e] > 0, so the root is unique, and since
|x-tau| = |e*sin(E0+x) - e*sin(E0)| <= 2e it lies in [tau-2e, tau+2e].
M=E-e*sin(E) is the case c=e, s=0. References: Vallado, Fundamentals of
Astrodynamics and Applications, 4th ed., Sec. 2.2-2.3; Battin (1999), Sec. 4.3.
"""

import numpy as np

from ..dynamics.two_body import EARTH_MU, real_array, validate_mu

_EPS: float = float(np.finfo(float).eps)
# Bisection alone shrinks the width-4 bracket below float64 spacing in <60 steps.
_MAX_ITERATIONS: int = 100


def _reduce(angle: float) -> tuple[float, int]:
    """Return (angle - 2*pi*k, k) with the remainder in [-pi, pi]."""
    turns = round(angle / (2 * np.pi))
    return angle - 2 * np.pi * turns, turns


def _solve_reduced(target: float, c: float, s: float) -> float:
    """Solve x - c*sin(x) + s*(1-cos(x)) = target for |target| <= pi.

    Newton-Raphson, replaced by bisection whenever a step leaves the
    sign-maintained bracket. Requires e=hypot(c, s) < 1.
    """
    e = float(np.hypot(c, s))
    low, high = target - 2 * e, target + 2 * e
    x = target
    for _ in range(_MAX_ITERATIONS):
        residual = x - c * np.sin(x) + s * (1 - np.cos(x)) - target
        if residual == 0:
            return x
        if residual < 0:
            low = x
        else:
            high = x
        scale = _EPS * max(1.0, abs(x))
        if high - low <= 4 * scale:
            return (low + high) / 2
        step = residual / (1 - c * np.cos(x) + s * np.sin(x))
        candidate = x - step
        if not low < candidate < high:
            candidate = (low + high) / 2
        elif abs(step) <= 4 * scale:
            return candidate
        x = candidate
    raise RuntimeError("Kepler equation did not converge")


def _validate_eccentricity(eccentricity: float) -> float:
    e = real_array(eccentricity, "eccentricity")
    if e.shape != () or not 0 <= e < 1:
        raise ValueError("eccentricity must be a scalar in [0,1)")
    return float(e)


def solve_kepler(mean_anomaly: float, eccentricity: float) -> float:
    """Return eccentric anomaly E (rad) with E - e*sin(E) = M, for 0<=e<1.

    E is continuous in M and keeps its revolution count: no wrapping to
    [0,2*pi). M is reduced to [-pi,pi] internally, so absolute roundoff grows
    with |M|. Conditioning: dE/dM = 1/(1-e*cos(E)), up to 1/(1-e) near periapsis
    as e->1. Non-finite M or e outside [0,1) raises ValueError.
    """
    mean = real_array(mean_anomaly, "mean_anomaly")
    if mean.shape != ():
        raise ValueError("mean_anomaly must be a scalar in radians")
    e = _validate_eccentricity(eccentricity)
    reduced, turns = _reduce(float(mean))
    return _solve_reduced(reduced, e, 0.0) + 2 * np.pi * turns


def _scalar(value: float, name: str) -> float:
    array = real_array(value, name)
    if array.shape != ():
        raise ValueError(f"{name} must be a scalar")
    return float(array)


def true_to_eccentric_anomaly(true_anomaly: float, eccentricity: float) -> float:
    """Eccentric anomaly E (rad) from true anomaly nu (rad), for 0<=e<1.

    tan(E/2) = sqrt((1-e)/(1+e)) tan(nu/2), evaluated without the tangent as
    E = 2*atan2(sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2)) on nu reduced to
    [-pi,pi], so E is in [-pi,pi] there. The revolution count of nu is kept:
    E(nu + 2*pi*k) = E(nu) + 2*pi*k (Vallado, 4th ed., Sec. 2.2). Non-finite
    nu or e outside [0,1) raises ValueError.
    """
    nu = _scalar(true_anomaly, "true_anomaly")
    e = _validate_eccentricity(eccentricity)
    reduced, turns = _reduce(nu)
    half = reduced / 2
    eccentric = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(half), np.sqrt(1 + e) * np.cos(half))
    return float(eccentric) + 2 * np.pi * turns


def true_to_mean_anomaly(true_anomaly: float, eccentricity: float) -> float:
    """Mean anomaly M = E - e*sin(E) (rad) from true anomaly nu (rad), 0<=e<1.

    E comes from ``true_to_eccentric_anomaly``; M keeps the revolution count
    of nu, so it is the inverse of solve_kepler followed by E -> nu.
    """
    nu = _scalar(true_anomaly, "true_anomaly")
    e = _validate_eccentricity(eccentricity)
    reduced, turns = _reduce(nu)
    eccentric = true_to_eccentric_anomaly(reduced, e)
    return float(eccentric - e * np.sin(eccentric)) + 2 * np.pi * turns


def elliptic_time_of_flight(
    true_anomaly_start: float,
    true_anomaly_end: float,
    semi_major_axis: float,
    eccentricity: float,
    mu: float = EARTH_MU,
) -> float:
    """Time (s) to move forward on an ellipse from nu_start to the next nu_end.

    dt = ((M(nu_end) - M(nu_start)) mod 2*pi) / n, n = sqrt(mu/a^3), with M
    from ``true_to_mean_anomaly``. The result is the first arrival in
    [0, period): whole revolutions are not counted and equal angles (mod
    2*pi) give 0. Angles in radians, a > 0 in metres, 0<=e<1, mu in m^3/s^2.
    Two-body Keplerian motion only. Invalid input raises ValueError; a
    non-finite result raises RuntimeError.
    """
    start = _scalar(true_anomaly_start, "true_anomaly_start")
    end = _scalar(true_anomaly_end, "true_anomaly_end")
    a = _scalar(semi_major_axis, "semi_major_axis")
    if a <= 0:
        raise ValueError("semi_major_axis must be positive, in metres")
    e = _validate_eccentricity(eccentricity)
    gm = validate_mu(mu)
    sweep = (true_to_mean_anomaly(_reduce(end)[0], e) - true_to_mean_anomaly(_reduce(start)[0], e)) % (2 * np.pi)
    if sweep >= 2 * np.pi:  # a tiny negative difference can round to exactly 2*pi
        sweep = 0.0
    with np.errstate(over="ignore", invalid="ignore"):
        time = float(sweep * np.sqrt(a / gm) * a)
    if not np.isfinite(time):
        raise RuntimeError("time of flight exceeds floating-point range")
    return time


def radius_extrema_on_arc(
    true_anomaly_start: float,
    duration: float,
    semi_major_axis: float,
    eccentricity: float,
    mu: float = EARTH_MU,
) -> tuple[float, float]:
    """Minimum and maximum radius (m) on the Keplerian arc flown for ``duration`` seconds.

    r(nu) = a(1-e^2)/(1+e*cos(nu)) decreases from apoapsis to periapsis and
    increases back, so over an arc its extrema are at the end points or at an
    apsis inside the arc. With M0 = M(nu_start) reduced to [0, 2*pi) and
    sweep = n*duration, n = sqrt(mu/a^3):

        periapsis (M = 2*pi*k) is inside iff M0 + sweep >= 2*pi;
        apoapsis (M = pi + 2*pi*k) is inside iff M0 <= pi <= M0 + sweep
            or M0 + sweep >= 3*pi;
        both are inside when sweep >= 2*pi.

    The end radius is a(1 - e*cos(E)) with E = solve_kepler(M0 + sweep, e);
    an apsis that coincides with an end point is covered by that end point.
    The result is exact for two-body motion up to roundoff: it does not
    sample the arc, so a periapsis between two samples cannot be missed
    (RF-4.05). Angles in radians, a > 0 in metres, 0 <= e < 1, duration >= 0
    in seconds, mu in m^3/s^2. Invalid input raises ValueError; a non-finite
    result raises RuntimeError.
    """
    start = _scalar(true_anomaly_start, "true_anomaly_start")
    span = _scalar(duration, "duration")
    if span < 0:
        raise ValueError("duration must be nonnegative, in seconds")
    a = _scalar(semi_major_axis, "semi_major_axis")
    if a <= 0:
        raise ValueError("semi_major_axis must be positive, in metres")
    e = _validate_eccentricity(eccentricity)
    gm = validate_mu(mu)
    periapsis, apoapsis = a * (1 - e), a * (1 + e)
    with np.errstate(over="ignore", invalid="ignore"):
        sweep = float(span * np.sqrt(gm / a) / a)
        start_radius = float(a * (1 - e) * (1 + e) / (1 + e * np.cos(start)))
    if not (np.isfinite(sweep) and np.isfinite(start_radius) and np.isfinite(apoapsis)):
        raise RuntimeError("arc extrema exceed floating-point range")
    if sweep >= 2 * np.pi:
        return float(periapsis), float(apoapsis)
    mean_start = true_to_mean_anomaly(_reduce(start)[0], e) % (2 * np.pi)
    mean_end = mean_start + sweep
    end_radius = float(a * (1 - e * np.cos(solve_kepler(mean_end, e))))
    low, high = min(start_radius, end_radius), max(start_radius, end_radius)
    if mean_end >= 2 * np.pi:
        low = periapsis
    if mean_start <= np.pi <= mean_end or mean_end >= 3 * np.pi:
        high = apoapsis
    # Roundoff can leave an end-point radius an ulp outside [r_p, r_a].
    return float(max(periapsis, low)), float(min(apoapsis, high))
