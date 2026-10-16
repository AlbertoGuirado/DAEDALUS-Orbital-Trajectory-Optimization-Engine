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

from ..dynamics.two_body import real_array

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
