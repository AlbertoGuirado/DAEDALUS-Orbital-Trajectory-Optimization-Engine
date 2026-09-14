"""Classical elliptic elements; equations and singular conventions: docs/phase1.md."""

from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..dynamics.two_body import EARTH_MU, real_array, validate_mu, validate_state
from .mechanics import eccentricity_vector, specific_angular_momentum, specific_energy

# Dimensionless roundoff threshold, not a physical eccentricity/altitude limit.
SINGULAR_TOL: float = 32 * np.finfo(float).eps


@dataclass(frozen=True)
class OrbitalElements:
    """Elliptic (a,e,i,Omega,omega,nu): a in m, e dimensionless, angles radians.

    a>0, 0<=e<1, 0<=i<=pi; other angles may wrap. Near-parabolic and near-radial
    geometry is ill-conditioned. Hyperbolic/parabolic elements are unsupported.
    """

    semi_major_axis: float
    eccentricity: float
    inclination: float
    raan: float = 0.0
    argument_of_periapsis: float = 0.0
    true_anomaly: float = 0.0

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            value = real_array(getattr(self, name), name)
            if value.shape != ():
                raise ValueError(f"{name} must be a scalar")
            object.__setattr__(self, name, float(value))
        if self.semi_major_axis <= 0 or not 0 <= self.eccentricity < 1:
            raise ValueError("elliptic elements require a>0 and 0<=e<1")
        if not 0 <= self.inclination <= np.pi:
            raise ValueError("inclination must be in [0,pi] radians")


def elements_to_state(elements: OrbitalElements, mu: float = EARTH_MU) -> NDArray[np.float64]:
    """Return (6,) state in m and m/s from elliptic elements (m, radians).

    p=a(1-e^2), r_PQW=p/(1+e*cos(nu))*[cos(nu),sin(nu),0],
    v_PQW=sqrt(mu/p)*[-sin(nu),e+cos(nu),0]. Rotate by Rz(Omega)Rx(i)Rz(omega).
    """
    gm = validate_mu(mu)
    if not isinstance(elements, OrbitalElements):
        raise ValueError("elements must be OrbitalElements")
    a, e, inc = elements.semi_major_axis, elements.eccentricity, elements.inclination
    node, arg, nu = elements.raan, elements.argument_of_periapsis, elements.true_anomaly
    cO, sO, ci, si, cw, sw = np.cos(node), np.sin(node), np.cos(inc), np.sin(inc), np.cos(arg), np.sin(arg)
    rotation = np.array([
        [cO*cw-sO*sw*ci, -cO*sw-sO*cw*ci],
        [sO*cw+cO*sw*ci, -sO*sw+cO*cw*ci],
        [sw*si, cw*si],
    ])
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        try:
            p = a * (1-e) * (1+e)
            position = p / (1+e*np.cos(nu)) * np.array([np.cos(nu), np.sin(nu)])
            velocity = np.sqrt(gm/p) * np.array([-np.sin(nu), e+np.cos(nu)])
            result = np.concatenate((rotation @ position, rotation @ velocity))
            return validate_state(result)
        except (FloatingPointError, ValueError) as exc:
            raise RuntimeError("elements cannot be represented at this numerical scale") from exc


def _angle(first: NDArray, second: NDArray, normal: NDArray) -> float:
    first, second = first/np.linalg.norm(first), second/np.linalg.norm(second)
    return float(np.arctan2(np.dot(np.cross(first, second), normal), np.dot(first, second)) % (2*np.pi))


def state_to_elements(state: ArrayLike, mu: float = EARTH_MU) -> OrbitalElements:
    """Convert (6,) SI inertial state to elliptic elements in metres/radians.

    a=-mu/(2*epsilon), e=|v cross h / mu-r/|r||, i=atan2(|h_xy|,h_z).
    Angles are oriented about h and returned in [0,2*pi). Circular: omega=0;
    equatorial: Omega=0. Remaining angles encode latitude/longitude in the
    direction of motion, including retrograde equatorial orbits.
    Zero nodes/eccentricities use 32*machine-epsilon dimensionless thresholds.
    Radial, unbound and numerically parabolic states raise ValueError.
    """
    y, gm = validate_state(state), validate_mu(mu)
    r, v = y[:3], y[3:]
    radius, speed = np.linalg.norm(r), np.linalg.norm(v)
    h = specific_angular_momentum(y)
    hnorm = np.linalg.norm(h)
    if speed == 0 or hnorm <= SINGULAR_TOL * radius * speed:
        raise ValueError("radial states have no supported classical orbital elements")
    energy = specific_energy(y, gm)
    if energy >= -SINGULAR_TOL * gm/radius:
        raise ValueError("only bound elliptic states away from the parabolic limit are supported")
    ev = eccentricity_vector(y, gm)
    e = float(np.linalg.norm(ev))
    if e >= 1:
        raise ValueError("eccentricity must be less than one")
    hhat = h/hnorm
    n = np.array([-h[1], h[0], 0.])
    equatorial = np.linalg.norm(n)/hnorm <= SINGULAR_TOL
    circular = e <= SINGULAR_TOL
    inc = float(np.arctan2(np.hypot(h[0], h[1]), h[2]))
    raan = 0.0 if equatorial else float(np.arctan2(n[1], n[0]) % (2*np.pi))
    reference = np.array([1., 0., 0.]) if equatorial else n
    arg = 0.0 if circular else _angle(reference, ev, hhat)
    nu = _angle(reference if circular else ev, r, hhat)
    return OrbitalElements(-gm/(2*energy), 0.0 if circular else e, inc, raan, arg, nu)
