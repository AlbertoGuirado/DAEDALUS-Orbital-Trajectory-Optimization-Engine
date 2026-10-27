"""Hohmann and Tsiolkovsky references from docs/03 section 8.2 and 9.1, not the production code path."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from orbital_optimizer import EARTH_MU, hohmann_transfer, propellant_mass, STANDARD_GRAVITY

LEO_RADIUS = 6_678_000.0
GEO_RADIUS = 42_164_000.0


def hohmann_reference(r1: float, r2: float, mu: float = EARTH_MU) -> tuple[float, float, float]:
    """Closed-form Hohmann burns and time, written with the docs/03 formulas."""
    dv1 = np.sqrt(mu/r1)*(np.sqrt(2*r2/(r1+r2))-1)
    dv2 = np.sqrt(mu/r2)*(1-np.sqrt(2*r1/(r1+r2)))
    time = np.pi*np.sqrt(((r1+r2)/2)**3/mu)
    return abs(dv1), abs(dv2), time


class HohmannTests(unittest.TestCase):
    def test_burns_match_closed_form_for_leo_to_geo(self) -> None:
        transfer = hohmann_transfer(LEO_RADIUS, GEO_RADIUS)
        dv1, dv2, _ = hohmann_reference(LEO_RADIUS, GEO_RADIUS)
        # RF-2.02: each impulse is reported separately, plus their sum.
        self.assertAlmostEqual(transfer.first_delta_v, dv1, places=9)
        self.assertAlmostEqual(transfer.second_delta_v, dv2, places=9)
        self.assertAlmostEqual(transfer.total_delta_v, dv1+dv2, places=9)

    def test_transfer_time_is_half_period_of_independent_ellipse(self) -> None:
        # RV-2.02: the time equals pi*sqrt(a^3/mu) for a=(r1+r2)/2, computed here directly.
        transfer = hohmann_transfer(LEO_RADIUS, GEO_RADIUS)
        a = (LEO_RADIUS+GEO_RADIUS)/2
        self.assertAlmostEqual(transfer.transfer_time, np.pi*np.sqrt(a**3/EARTH_MU), places=6)
        self.assertAlmostEqual(transfer.transfer_semi_major_axis, a, places=6)

    def test_equal_radii_need_no_delta_v(self) -> None:
        # RV-2.03: as r2 approaches r1 the total cost tends to zero.
        self.assertEqual(hohmann_transfer(LEO_RADIUS, LEO_RADIUS).total_delta_v, 0.0)
        totals = [hohmann_transfer(LEO_RADIUS, LEO_RADIUS*(1+eps)).total_delta_v for eps in (1e-1, 1e-2, 1e-3)]
        self.assertTrue(totals[0] > totals[1] > totals[2] > 0)

    def test_direction_does_not_change_magnitudes(self) -> None:
        up = hohmann_transfer(LEO_RADIUS, GEO_RADIUS)
        down = hohmann_transfer(GEO_RADIUS, LEO_RADIUS)
        self.assertAlmostEqual(up.total_delta_v, down.total_delta_v, places=9)
        self.assertAlmostEqual(up.transfer_time, down.transfer_time, places=9)

    def test_invalid_radii_are_rejected(self) -> None:
        for bad in (0.0, -1.0, np.nan, np.inf, [1.0, 2.0]):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                hohmann_transfer(bad, GEO_RADIUS)
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                hohmann_transfer(LEO_RADIUS, bad)


class TsiolkovskyTests(unittest.TestCase):
    def test_propellant_matches_mass_ratio_formula(self) -> None:
        # RF-2.04: m_p = m0*(1-exp(-dv/(Isp*g0))), with the conventional g0.
        m0, dv, isp = 1000.0, 3000.0, 320.0
        expected = m0*(1-np.exp(-dv/(isp*9.80665)))
        self.assertAlmostEqual(propellant_mass(m0, dv, isp), expected, places=9)
        self.assertAlmostEqual(STANDARD_GRAVITY, 9.80665, places=12)

    def test_propellant_consumed_for_leo_to_geo(self) -> None:
        dv_total = hohmann_transfer(LEO_RADIUS, GEO_RADIUS).total_delta_v
        mass = propellant_mass(2000.0, dv_total, 320.0)
        final = 2000.0 - mass
        # Tsiolkovsky in its exponential form must return the same final mass.
        assert_allclose(final, 2000.0*np.exp(-dv_total/(320.0*STANDARD_GRAVITY)), rtol=1e-12)
        self.assertTrue(0 < mass < 2000.0)

    def test_zero_delta_v_consumes_no_propellant(self) -> None:
        self.assertEqual(propellant_mass(1000.0, 0.0, 300.0), 0.0)

    def test_invalid_inputs_are_rejected(self) -> None:
        cases = [(0.0, 1.0, 300.0), (-1.0, 1.0, 300.0), (1000.0, -1.0, 300.0), (1000.0, 1.0, 0.0), (1000.0, [1.0, 2.0], 300.0)]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                propellant_mass(*args)


if __name__ == "__main__":
    unittest.main()
