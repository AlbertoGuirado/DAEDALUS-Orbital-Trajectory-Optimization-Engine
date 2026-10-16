"""Analytic two-body propagation with Lagrange f and g coefficients, SI units.

Exact for the Newtonian point-mass model up to floating-point roundoff; it is
a reference for numerical propagation, not an extensible force model.
"""

import numpy as np
from numpy.typing import ArrayLike

from ..dynamics.two_body import EARTH_MU, validate_mu, validate_state
from ..orbital.elements import state_to_elements
from ..orbital.kepler import _reduce, _solve_reduced
from ..orbital.mechanics import specific_energy
from .propagator import PropagationResult, _positive_scalar, _sample_times


def propagate_kepler(
    initial_state: ArrayLike,
    duration: float,
    mu: float = EARTH_MU,
    *,
    sample_times: ArrayLike | None = None,
) -> PropagationResult:
    """Propagate ``[x, y, z, vx, vy, vz]`` analytically for ``duration`` seconds.

    Positions are metres, velocities m/s, mu m^3/s^2. From a=-mu/(2*epsilon),
    n=sqrt(mu/a^3), c=1-r0/a and s=(r0.v0)/sqrt(mu*a), solve for the change in
    eccentric anomaly x at each time t (see orbital/kepler.py):

        n*t = x - c*sin(x) + s*(1-cos(x)),   r = a*(1 - c*cos(x) + s*sin(x)),
        f = 1 - a/r0*(1-cos(x)),             g = t - (x - sin(x))/n,
        fdot = -sqrt(mu*a)*sin(x)/(r*r0),    gdot = 1 - a/r*(1-cos(x)),
        r(t) = f*r0 + g*v0,                  v(t) = fdot*r0 + gdot*v0.

    No orbital angles are formed, so circular and equatorial orbits are not
    singular. Whole revolutions are removed from n*t and t before evaluating
    f and g. Roundoff grows linearly with elapsed time: a is derived from the
    rounded state with relative error up to ~4*eps*(3+e)/(1-e), so the period,
    and hence the along-track timing, drifts by 1.5*t*delta_a/a. This is the
    conditioning of the initial-value problem, not an integration error.

    The domain matches ``propagate``: bound, nonradial elliptic states. There
    is no body surface or collision model. Sampling follows ``propagate``:
    None returns only [0, duration]; otherwise strictly increasing times in
    [0, duration], with missing endpoints added. ``nfev`` is 0 because no
    derivative is evaluated; method is "Kepler". Bad input raises ValueError.
    """
    state = validate_state(initial_state)
    gm = validate_mu(mu)
    state_to_elements(state, gm)
    duration = _positive_scalar(duration, "duration")
    samples = _sample_times(sample_times, duration)
    times = np.array([0.0, duration]) if samples is None else samples

    r0, v0 = state[:3], state[3:]
    radius0 = float(np.linalg.norm(r0))
    a = -gm / (2 * specific_energy(state, gm))
    n = float(np.sqrt(gm / a**3))
    period = 2 * np.pi / n
    c = 1 - radius0 / a
    s = float(np.dot(r0, v0)) / np.sqrt(gm * a)

    states = np.empty((len(times), 6))
    for row, t in enumerate(times):
        reduced, turns = _reduce(n * t)
        x = _solve_reduced(reduced, c, s)
        sin_x, one_minus_cos = np.sin(x), 1 - np.cos(x)
        radius = a * (1 - c * np.cos(x) + s * sin_x)
        f = 1 - a / radius0 * one_minus_cos
        g = (t - turns * period) - (x - sin_x) / n
        fdot = -np.sqrt(gm * a) * sin_x / (radius * radius0)
        gdot = 1 - a / radius * one_minus_cos
        states[row] = np.concatenate((f * r0 + g * v0, fdot * r0 + gdot * v0))
    if not np.all(np.isfinite(states)):
        raise RuntimeError("Kepler propagation produced non-finite states")
    return PropagationResult(times, states, "Kepler", 0, gm)
