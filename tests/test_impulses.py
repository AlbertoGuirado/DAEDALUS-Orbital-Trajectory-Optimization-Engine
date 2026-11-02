"""RV-2.01: Hohmann burns applied to the numerical propagator reach the target circular orbit."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from orbital_optimizer import (
    EARTH_MU, hohmann_transfer, propagate, propagate_with_impulses, state_to_elements, total_delta_v,
)

LEO_RADIUS = 6_678_000.0
GEO_RADIUS = 42_164_000.0


class ImpulsivePropagationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.transfer = hohmann_transfer(LEO_RADIUS, GEO_RADIUS)
        v_leo = np.sqrt(EARTH_MU/LEO_RADIUS)
        self.state = np.array([LEO_RADIUS, 0, 0, 0, v_leo, 0])
        self.first_burn = np.array([0, self.transfer.first_delta_v, 0])

    def test_hohmann_burns_on_propagator_reach_geo_circular_orbit(self) -> None:
        # The second burn is taken along the velocity at apogee of the ballistic transfer,
        # which lies on the transfer ellipse; its magnitude is the Hohmann value.
        T = self.transfer.transfer_time
        apogee = propagate_with_impulses(self.state, [(0.0, self.first_burn)], T)
        self.assertAlmostEqual(np.linalg.norm(apogee[:3]), GEO_RADIUS, delta=1.0)
        second_burn = apogee[3:]/np.linalg.norm(apogee[3:])*self.transfer.second_delta_v

        final = propagate_with_impulses(
            self.state, [(0.0, self.first_burn), (T, second_burn)], T,
        )
        elements = state_to_elements(final, EARTH_MU)
        self.assertAlmostEqual(elements.semi_major_axis, GEO_RADIUS, delta=1.0)
        self.assertLess(elements.eccentricity, 1e-6)
        v_geo = np.sqrt(EARTH_MU/GEO_RADIUS)
        assert_allclose(final[:3], [-GEO_RADIUS, 0, 0], atol=1.0)
        assert_allclose(final[3:], [0, -v_geo, 0], atol=1e-3)

    def test_burn_costs_match_hohmann_total(self) -> None:
        second = np.array([0, -self.transfer.second_delta_v, 0])
        self.assertAlmostEqual(
            total_delta_v([self.first_burn, second]), self.transfer.total_delta_v, places=9,
        )

    def test_no_burns_equals_free_flight(self) -> None:
        duration = self.transfer.transfer_time
        free = propagate(self.state, duration).final_state
        assert_allclose(propagate_with_impulses(self.state, [], duration), free, rtol=1e-12, atol=1e-9)

    def test_burn_at_start_and_end_are_applied(self) -> None:
        duration = 100.0
        kicked = propagate_with_impulses(self.state, [(0.0, self.first_burn)], duration)
        expected = propagate(self.state + np.r_[0, 0, 0, self.first_burn], duration).final_state
        assert_allclose(kicked, expected, rtol=1e-12, atol=1e-9)
        final_kick = propagate_with_impulses(self.state, [(duration, self.first_burn)], duration)
        assert_allclose(final_kick, propagate(self.state, duration).final_state + np.r_[0, 0, 0, self.first_burn], atol=1e-6)

    def test_invalid_burn_schedules_are_rejected(self) -> None:
        duration = 100.0
        cases = [
            [(50.0, self.first_burn), (10.0, self.first_burn)],
            [(50.0, self.first_burn), (50.0, self.first_burn)],
            [(-1.0, self.first_burn)],
            [(101.0, self.first_burn)],
            [(10.0, [0, 1.0])],
        ]
        for burns in cases:
            with self.subTest(burns=burns), self.assertRaises(ValueError):
                propagate_with_impulses(self.state, burns, duration)

    def test_nonpositive_duration_is_rejected(self) -> None:
        for duration in (0.0, -1.0):
            with self.subTest(duration=duration), self.assertRaises(ValueError):
                propagate_with_impulses(self.state, [], duration)


if __name__ == "__main__":
    unittest.main()
