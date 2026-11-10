"""Additive Kepler helpers (Phase 3, Q6): true -> eccentric/mean anomaly and time of flight.

References are independent identities or closed forms, not values produced by
the code under test. Tolerances: the helpers use ~10 float64 operations and
the angle maps have derivative at most sqrt((1+e)/(1-e)) (<= 14.1 for
e <= 0.99), so their error is below ~10*14*eps ~ 3e-14 rad; the angle
tolerance 1e-12 rad leaves a margin of about 30.
"""
import unittest

import numpy as np

from orbital_optimizer import (
    EARTH_MU, OrbitalElements, elements_to_state, hohmann_transfer, orbital_period,
    propagate_kepler, solve_kepler, state_to_elements,
)
from orbital_optimizer.orbital.kepler import (
    elliptic_time_of_flight, true_to_eccentric_anomaly, true_to_mean_anomaly,
)

ANGLE_TOL = 1e-12
ECCENTRICITIES = (0.0, 0.1, 0.5, 0.9, 0.99)


class AnomalyConversionTests(unittest.TestCase):
    def test_circular_orbit_anomalies_coincide(self) -> None:
        for nu in (-3.0, -1.0, 0.0, 0.5, 3.1, 7.0, -20.0):
            with self.subTest(nu=nu):
                self.assertAlmostEqual(true_to_eccentric_anomaly(nu, 0.0), nu, delta=ANGLE_TOL)
                self.assertAlmostEqual(true_to_mean_anomaly(nu, 0.0), nu, delta=ANGLE_TOL)

    def test_eccentric_anomaly_satisfies_geometric_identities(self) -> None:
        # cos E = (e + cos nu)/(1 + e cos nu), sin E = sqrt(1-e^2) sin nu/(1 + e cos nu).
        for e in ECCENTRICITIES:
            for nu in np.linspace(-np.pi, np.pi, 37):
                with self.subTest(e=e, nu=nu):
                    eccentric = true_to_eccentric_anomaly(nu, e)
                    denominator = 1 + e * np.cos(nu)
                    self.assertAlmostEqual(np.cos(eccentric), (e + np.cos(nu)) / denominator, delta=ANGLE_TOL)
                    self.assertAlmostEqual(np.sin(eccentric), np.sqrt(1 - e * e) * np.sin(nu) / denominator,
                                           delta=ANGLE_TOL)
                    self.assertLessEqual(abs(eccentric), np.pi)

    def test_round_trip_with_solve_kepler(self) -> None:
        rng = np.random.Generator(np.random.PCG64(42))
        for e in (0.0, 0.1, 0.5, 0.9):
            for mean in rng.uniform(-np.pi, np.pi, 50):
                with self.subTest(e=e, mean=mean):
                    eccentric = solve_kepler(mean, e)
                    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(eccentric / 2), np.sqrt(1 - e) * np.cos(eccentric / 2))
                    self.assertAlmostEqual(true_to_mean_anomaly(nu, e), mean, delta=ANGLE_TOL)

    def test_symmetry_and_revolutions(self) -> None:
        for e in ECCENTRICITIES:
            with self.subTest(e=e):
                self.assertEqual(true_to_mean_anomaly(0.0, e), 0.0)
                self.assertAlmostEqual(true_to_mean_anomaly(np.pi, e), np.pi, delta=ANGLE_TOL)
                for nu in (0.3, 2.0, 3.0):
                    self.assertAlmostEqual(true_to_mean_anomaly(-nu, e), -true_to_mean_anomaly(nu, e), delta=ANGLE_TOL)
                    for turns in (1, -2, 5):
                        self.assertAlmostEqual(true_to_mean_anomaly(nu + 2 * np.pi * turns, e),
                                               true_to_mean_anomaly(nu, e) + 2 * np.pi * turns, delta=1e-11)

    def test_invalid_inputs(self) -> None:
        for function in (true_to_eccentric_anomaly, true_to_mean_anomaly):
            for nu, e in ((np.nan, 0.1), (0.1, 1.0), (0.1, -0.1), ([0.1, 0.2], 0.1), (0.1, np.inf)):
                with self.subTest(function=function.__name__, nu=nu, e=e), self.assertRaises(ValueError):
                    function(nu, e)


class TimeOfFlightTests(unittest.TestCase):
    def test_hohmann_half_ellipse(self) -> None:
        r1, r2 = 6_678_000.0, 42_164_000.0
        a, e = (r1 + r2) / 2, (r2 - r1) / (r2 + r1)
        expected = np.pi * np.sqrt(a**3 / EARTH_MU)   # independent closed form
        self.assertAlmostEqual(elliptic_time_of_flight(0.0, np.pi, a, e), expected, places=6)
        self.assertAlmostEqual(elliptic_time_of_flight(0.0, np.pi, a, e), hohmann_transfer(r1, r2).transfer_time, places=6)

    def test_circular_orbit_is_uniform(self) -> None:
        a = 7_000_000.0
        period = orbital_period(a)
        self.assertAlmostEqual(elliptic_time_of_flight(0.25, 1.75, a, 0.0), 1.5 / (2 * np.pi) * period, places=6)

    def test_forward_motion_and_period_closure(self) -> None:
        a, e = 20_000_000.0, 0.6
        period = orbital_period(a)
        outbound = elliptic_time_of_flight(0.0, np.pi, a, e)
        inbound = elliptic_time_of_flight(np.pi, 0.0, a, e)
        self.assertAlmostEqual(outbound + inbound, period, places=6)
        self.assertAlmostEqual(outbound, inbound, places=6)
        self.assertEqual(elliptic_time_of_flight(1.0, 1.0, a, e), 0.0)
        # Going backwards in angle means almost a full revolution forward.
        self.assertAlmostEqual(elliptic_time_of_flight(1.0, 0.9, a, e) + elliptic_time_of_flight(0.9, 1.0, a, e),
                               period, places=6)
        # Equivalent angles (mod 2 pi) give the same first arrival.
        self.assertAlmostEqual(elliptic_time_of_flight(0.3 + 2 * np.pi, 2.5 - 4 * np.pi, a, e),
                               elliptic_time_of_flight(0.3, 2.5, a, e), places=6)

    def test_agrees_with_analytic_propagation(self) -> None:
        # Propagating for the computed time must arrive at the requested true anomaly.
        # Bound: Kepler timing error 1.5*t*4*eps*(3+e)/(1-e) (tests/test_kepler.py model) times
        # dnu/dM <= (1+e)^2/(1-e^2)^1.5 gives ~1.4e-11 rad for e=0.9 over one period; 1e-9 leaves ~70x.
        for e, start, end in ((0.6, 0.3, 2.5), (0.2, 5.0, 1.0), (0.9, -0.5, 0.4)):
            with self.subTest(e=e, start=start, end=end):
                orbit = OrbitalElements(15_000_000.0, e, 0.4, 1.0, 2.0, start)
                time = elliptic_time_of_flight(start, end, orbit.semi_major_axis, e)
                arrival = state_to_elements(propagate_kepler(elements_to_state(orbit), time).final_state)
                difference = (arrival.true_anomaly - end + np.pi) % (2 * np.pi) - np.pi
                self.assertLess(abs(difference), 1e-9)

    def test_invalid_inputs(self) -> None:
        for args in ((0.0, 1.0, 0.0, 0.1), (0.0, 1.0, -1.0, 0.1), (0.0, 1.0, 1e7, 1.0),
                     (np.nan, 1.0, 1e7, 0.1), (0.0, 1.0, [1e7], 0.1)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                elliptic_time_of_flight(*args)
        with self.assertRaises(ValueError):
            elliptic_time_of_flight(0.0, 1.0, 1e7, 0.1, mu=0.0)


if __name__ == "__main__":
    unittest.main()
