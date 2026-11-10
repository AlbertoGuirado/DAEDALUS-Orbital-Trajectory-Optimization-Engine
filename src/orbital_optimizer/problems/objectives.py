"""Weighted transfer objective with delta-v, time, energy and risk terms (RF-3.19 to RF-3.21).

    J = w_dv dv/dv_ref + w_t t/t_ref + w_E E/E_ref + w_R R/R_ref

Every term is dimensionless: the metric divided by its value on the Hohmann
transfer between the same two orbits, so each term equals 1 at the analytic
baseline (docs/03 section 10.3). The metrics are those reported by
problems.coplanar_multi_impulse.CoplanarMultiImpulseProblem:

    delta_v   delta_v_total_m_s            sum_k |dv_k|, m/s
    time      transfer_time_s              instant of the last impulse, s
    energy    impulse_kinetic_energy_j_kg  E = sum_k |dv_k|^2 / 2, J/kg
    risk      closest_approach_index       R = min(r1, r2) / min_t r(t), dimensionless

Meaning and limits of the two proxies (they are modelling choices, not
validated models):
  - E is the kinetic energy per unit mass that each impulse imparts in the
    frame moving with the vehicle just before the burn. For a fixed total
    delta-v it is lowest when the burns are equal, so it penalizes
    concentrating the velocity change in one burn. It is not the chemical
    energy of the propellant, which follows the delta-v (Tsiolkovsky).
  - R compares the closest approach to the central body with the lower of the
    two mission orbits: R = 1 if the transfer never descends below it (as in
    the Hohmann transfer) and R > 1 if it does. It is a geometric proxy for
    atmospheric and collision exposure; it is not a probability.

Weights are dimensionless, nonnegative mission preferences. Their rationale
must be supplied by the caller and is recorded with the objective (RF-3.20).
With only the delta-v weight positive the minimizer is that of J = delta-v.
A weighted sum reaches only the convex part of a Pareto front.
"""

from __future__ import annotations

from ..optimization.objective import ObjectiveTerm, WeightedObjective
from ..orbital.transfers import HohmannTransfer

TERM_METRICS: dict[str, tuple[str, str]] = {
    "delta_v": ("delta_v_total_m_s", "m/s"),
    "time": ("transfer_time_s", "s"),
    "energy": ("impulse_kinetic_energy_j_kg", "J/kg"),
    "risk": ("closest_approach_index", "dimensionless"),
}


def transfer_objective(reference: HohmannTransfer, *, delta_v: float = 1.0, time: float = 0.0,
                       energy: float = 0.0, risk: float = 0.0, rationale: str) -> WeightedObjective:
    """WeightedObjective with the four transfer terms, referenced to a Hohmann transfer.

    ``reference`` is orbital.transfers.hohmann_transfer(r1, r2, mu) for the
    same orbits as the problem. References: dv_ref = its total delta-v (m/s),
    t_ref = its transfer time (s), E_ref = (dv1^2 + dv2^2)/2 of its two burns
    (J/kg) and R_ref = 1 (its closest approach is the lower orbit).
    ``rationale`` explains the weights; it is required and recorded in every
    term. Invalid weights or an empty rationale raise ValueError.
    """
    if not isinstance(reference, HohmannTransfer):
        raise ValueError("reference must be a HohmannTransfer")
    if not isinstance(rationale, str) or not rationale.strip():
        raise ValueError("rationale must be a nonempty string explaining the weights")
    references = {
        "delta_v": (reference.total_delta_v, "Hohmann total delta-v of the same orbits"),
        "time": (reference.transfer_time, "Hohmann transfer time of the same orbits"),
        "energy": ((reference.first_delta_v**2 + reference.second_delta_v**2) / 2,
                   "impulse kinetic energy of the two Hohmann burns"),
        "risk": (1.0, "closest-approach index of the Hohmann transfer"),
    }
    weights = {"delta_v": delta_v, "time": time, "energy": energy, "risk": risk}
    terms = []
    for name, (metric, unit) in TERM_METRICS.items():
        scale, source = references[name]
        terms.append(ObjectiveTerm(name, metric, weights[name], scale, unit,
                                   f"weight: {rationale.strip()}; reference: {source}"))
    return WeightedObjective(terms)
