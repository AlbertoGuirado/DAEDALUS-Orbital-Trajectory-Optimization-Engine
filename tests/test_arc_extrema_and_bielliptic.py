"""Radius extrema on a Keplerian arc and the bi-elliptic transfer (RF-2.05).

References independent of the implementation:
  - radius_extrema_on_arc: the conic equation r(nu) = p/(1 + e cos nu) at hand
    picked anomalies, and a brute-force minimum over a dense sampling of the
    arc obtained with propagate_kepler.
  - bielliptic_transfer: the Hohmann closed form (degenerate case), the
    bi-parabolic limit (sqrt(2) - 1)(v_c1 + v_c2) and the classical radius
    ratio thresholds 11.94 and 15.58 (docs/03 section 8.3).
"""
import math
import unittest

import numpy as np

from orbital_optimizer import (
    EARTH_MU, BiellipticTransfer, OrbitalElements, bielliptic_transfer, elements_to_state,
    hohmann_transfer, orbital_period, propagate_kepler,
)
from orbital_optimizer.orbital.kepler import elliptic_time_of_flight, radius_extrema_on_arc

A, E = 1.0e7, 0.5
P = A * (1 - E * E)


def radius(nu_deg: float) -> float:
    return P / (1 + E * math.cos(math.radians(nu_deg)))


def flight(start_deg: float, end_deg: float) -> float:
    return elliptic_time_of_flight(math.radians(start_deg), math.radians(end_deg), A, E)


class ArcExtremaTests(unittest.TestCase):
    def assert_extrema(self, start_deg: float, duration: float, low: float, high: float) -> None:
        got_low, got_high = radius_extrema_on_arc(math.radians(start_deg), duration, A, E)
        # r/a is O(1): the conic equation and the Kepler solution agree to a few ulps of a.
        self.assertAlmostEqual(got_low, low, delta=1e-6)
        self.assertAlmostEqual(got_high, high, delta=1e-6)

    def test_arc_without_an_apsis_uses_its_end_points(self) -> None:
        self.assert_extrema(20.0, flight(20.0, 100.0), radius(20.0), radius(100.0))      # outbound
        self.assert_extrema(200.0, flight(200.0, 300.0), radius(300.0), radius(200.0))   # inbound

    def test_arc_through_periapsis(self) -> None:
        # Both end points are valid and equal by symmetry; the minimum is in between.
        self.assert_extrema(300.0, flight(300.0, 60.0), A * (1 - E), radius(60.0))

    def test_arc_through_apoapsis(self) -> None:
        self.assert_extrema(150.0, flight(150.0, 210.0), radius(150.0), A * (1 + E))

    def test_arc_through_both_apsides(self) -> None:
        self.assert_extrema(170.0, flight(170.0, 10.0), A * (1 - E), A * (1 + E))
        self.assert_extrema(350.0, flight(350.0, 190.0), A * (1 - E), A * (1 + E))

    def test_end_points_on_an_apsis_and_degenerate_arcs(self) -> None:
        self.assert_extrema(0.0, flight(0.0, 90.0), A * (1 - E), radius(90.0))
        self.assert_extrema(90.0, flight(90.0, 180.0), radius(90.0), A * (1 + E))
        self.assert_extrema(40.0, 0.0, radius(40.0), radius(40.0))
        period = orbital_period(A)
        for duration in (period, 1.5 * period, 7 * period):
            self.assert_extrema(40.0, duration, A * (1 - E), A * (1 + E))
        self.assertEqual(radius_extrema_on_arc(1.0, 100.0, A, 0.0), (A, A))      # circular

    def test_matches_dense_sampling_of_the_propagated_arc(self) -> None:
        rng = np.random.Generator(np.random.PCG64(4))
        for _ in range(40):
            a = 10 ** rng.uniform(6.9, 8.0)
            e = rng.uniform(0.0, 0.9)
            nu = rng.uniform(0.0, 2 * np.pi)
            duration = rng.uniform(0.01, 1.3) * orbital_period(a)
            state = elements_to_state(OrbitalElements(a, e, 0.4, 1.0, 2.0, nu))
            samples = np.linspace(0.0, duration, 2001)
            radii = np.linalg.norm(propagate_kepler(state, duration, sample_times=samples).states[:, :3], axis=1)
            low, high = radius_extrema_on_arc(nu, duration, a, e)
            # Sampled extrema lie inside the exact ones; between samples the radius changes by at
            # most max|dr/dt| dt = e sqrt(mu/p) dt.
            step = e * math.sqrt(EARTH_MU / (a * (1 - e * e))) * duration / 2000 + 1e-6 * a
            with self.subTest(a=a, e=e, nu=nu, duration=duration):
                self.assertLessEqual(low, radii.min() + 1e-9 * a)
                self.assertGreaterEqual(high, radii.max() - 1e-9 * a)
                self.assertLessEqual(radii.min() - low, step)
                self.assertLessEqual(high - radii.max(), step)

    def test_invalid_inputs(self) -> None:
        for arguments in ((math.nan, 1.0, A, E), (0.0, -1.0, A, E), (0.0, math.inf, A, E), (0.0, 1.0, 0.0, E),
                          (0.0, 1.0, -A, E), (0.0, 1.0, A, 1.0), (0.0, 1.0, A, -0.1), (0.0, [1.0, 2.0], A, E)):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                radius_extrema_on_arc(*arguments)
        with self.assertRaises(ValueError):
            radius_extrema_on_arc(0.0, 1.0, A, E, mu=0.0)


def hohmann_cost(ratio: float) -> float:
    """Hohmann total delta-v in units of v_c1, for R = r2/r1 > 1 (docs/03 section 8.2)."""
    return math.sqrt(2 * ratio / (1 + ratio)) - 1 + (1 - math.sqrt(2 / (1 + ratio))) / math.sqrt(ratio)


def biparabolic_cost(ratio: float) -> float:
    """Limit of the bi-elliptic cost for rb -> infinity, in units of v_c1."""
    return (math.sqrt(2) - 1) * (1 + 1 / math.sqrt(ratio))


class BiellipticTests(unittest.TestCase):
    R1 = 6_678_000.0

    def cost(self, ratio: float, factor: float) -> float:
        """Bi-elliptic total delta-v in units of v_c1 with rb = factor * r2."""
        transfer = bielliptic_transfer(self.R1, ratio * self.R1, factor * ratio * self.R1)
        return transfer.total_delta_v / math.sqrt(EARTH_MU / self.R1)

    def test_degenerates_into_hohmann(self) -> None:
        for r1, r2 in ((self.R1, 42_164_000.0), (42_164_000.0, self.R1)):
            reference = hohmann_transfer(r1, r2)
            transfer = bielliptic_transfer(r1, r2, max(r1, r2))
            self.assertIsInstance(transfer, BiellipticTransfer)
            self.assertAlmostEqual(transfer.total_delta_v, reference.total_delta_v, places=9)
            # With rb on the larger circle, one "transfer ellipse" is that circle itself, flown for
            # half a revolution, and the other is the Hohmann ellipse.
            self.assertAlmostEqual(transfer.transfer_time,
                                   reference.transfer_time + orbital_period(max(r1, r2)) / 2, places=5)
            # One burn is null: the first when descending from rb = r1, the last when ascending to rb = r2.
            null = transfer.third_delta_v if r2 > r1 else transfer.first_delta_v
            self.assertEqual(null, 0.0)

    def test_closed_form_values_and_time(self) -> None:
        r1, r2, rb, mu = 7.0e6, 1.05e8, 2.1e8, EARTH_MU
        a1, a2 = (r1 + rb) / 2, (r2 + rb) / 2
        speed = lambda r, a: math.sqrt(mu * (2 / r - 1 / a))  # noqa: E731
        transfer = bielliptic_transfer(r1, r2, rb)
        self.assertAlmostEqual(transfer.first_delta_v, speed(r1, a1) - math.sqrt(mu / r1), places=9)
        self.assertAlmostEqual(transfer.second_delta_v, speed(rb, a2) - speed(rb, a1), places=9)
        self.assertAlmostEqual(transfer.third_delta_v, speed(r2, a2) - math.sqrt(mu / r2), places=9)
        self.assertAlmostEqual(transfer.total_delta_v,
                               transfer.first_delta_v + transfer.second_delta_v + transfer.third_delta_v, places=9)
        expected_time = math.pi * (math.sqrt(a1**3 / mu) + math.sqrt(a2**3 / mu))
        self.assertAlmostEqual(transfer.transfer_time / expected_time, 1.0, places=13)
        self.assertEqual((transfer.first_semi_major_axis, transfer.second_semi_major_axis), (a1, a2))

    def test_biparabolic_limit(self) -> None:
        for ratio in (2.0, 12.0, 50.0):
            # Leading correction is of relative order r2/rb = 1e-9.
            self.assertAlmostEqual(self.cost(ratio, 1e9), biparabolic_cost(ratio), delta=1e-7)

    def test_radius_ratio_thresholds(self) -> None:
        # 11.94: the bi-parabolic limit equals Hohmann. Located here by bisection on the closed forms.
        low, high = 5.0, 20.0
        for _ in range(60):
            middle = (low + high) / 2
            if biparabolic_cost(middle) > hohmann_cost(middle):
                low = middle
            else:
                high = middle
        self.assertAlmostEqual(low, 11.94, delta=0.005)
        factors = (1.001, 1.5, 10.0, 1e3, 1e6)
        # Below 11.94 Hohmann is cheaper for every intermediate radius.
        self.assertTrue(all(self.cost(11.5, factor) > hohmann_cost(11.5) for factor in factors))
        # Between the thresholds it depends on rb: worse just above r2, better far away.
        self.assertGreater(self.cost(13.0, 1.5), hohmann_cost(13.0))
        self.assertLess(self.cost(13.0, 1e3), hohmann_cost(13.0))
        # Above 15.58 the bi-elliptic transfer is cheaper for every rb > r2.
        self.assertTrue(all(self.cost(16.0, factor) < hohmann_cost(16.0) for factor in factors))
        # The second threshold is where rb slightly above r2 starts to pay off.
        self.assertGreater(self.cost(15.4, 1.001), hohmann_cost(15.4))
        self.assertLess(self.cost(15.8, 1.001), hohmann_cost(15.8))

    def test_normalized_costs_match_the_implementation(self) -> None:
        for ratio in (1.5, 6.3, 11.94, 40.0):
            reference = hohmann_transfer(self.R1, ratio * self.R1)
            self.assertAlmostEqual(reference.total_delta_v / math.sqrt(EARTH_MU / self.R1), hohmann_cost(ratio), places=12)

    def test_invalid_inputs(self) -> None:
        for arguments in ((self.R1, 2 * self.R1, 1.5 * self.R1), (2 * self.R1, self.R1, 1.5 * self.R1),
                          (-1.0, self.R1, self.R1), (self.R1, 0.0, self.R1), (self.R1, self.R1, math.nan)):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                bielliptic_transfer(*arguments)
        with self.assertRaises(ValueError):
            bielliptic_transfer(self.R1, 2 * self.R1, 3 * self.R1, mu=-1.0)
        with self.assertRaises(RuntimeError):
            bielliptic_transfer(1e-10, 2e-10, 4e-10, mu=1e308)


if __name__ == "__main__":
    unittest.main()
