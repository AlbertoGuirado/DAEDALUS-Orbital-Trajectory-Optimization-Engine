"""Analytic physics checks. Units and predeclared tolerances: docs/phase1.md."""
import unittest
import numpy as np
from numpy.testing import assert_allclose

from orbital_optimizer import (
    EARTH_MU, delta_v, eccentricity_vector, node_vector, orbital_period,
    specific_angular_momentum, specific_energy, total_delta_v, two_body_dynamics,
)
from orbital_optimizer.dynamics.two_body import validate_mu, validate_state


class DynamicsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.r = 7e6
        self.v = np.sqrt(EARTH_MU/self.r)
        self.state = np.array([self.r, 0, 0, 0, self.v, 0])

    def test_circular_analytic_quantities(self) -> None:
        assert_allclose(two_body_dynamics(0, self.state),
                        [0, self.v, 0, -EARTH_MU/self.r**2, 0, 0], rtol=1e-14)
        self.assertAlmostEqual(specific_energy(self.state)/(-EARTH_MU/(2*self.r)), 1., places=14)
        assert_allclose(specific_angular_momentum(self.state), [0, 0, self.r*self.v], rtol=1e-14)
        assert_allclose(eccentricity_vector(self.state), 0, atol=1e-14)
        assert_allclose(node_vector(self.state), 0, atol=0)
        self.assertAlmostEqual(orbital_period(self.r)/(2*np.pi*self.r/self.v), 1., places=14)

    def test_custom_mu(self) -> None:
        assert_allclose(two_body_dynamics(0, [1,0,0,0,1,0], 1), [0,1,0,-1,0,0])
        self.assertAlmostEqual(orbital_period(1, 1), 2*np.pi)

    def test_impulse_costs_do_not_cancel(self) -> None:
        self.assertEqual(delta_v([0,0,0], [3,4,0]), 5.)
        self.assertEqual(total_delta_v([[3,4,0], [-3,-4,0]]), 10.)
        self.assertEqual(total_delta_v([]), 0.)
        self.assertEqual(total_delta_v(np.empty((0,3))), 0.)

    def test_invalid_states(self) -> None:
        for value in ([0]*6, [1]*5, [[1]*6], [np.nan]*6, [np.inf]*6,
                      np.ones(6, dtype=complex), ["1"]*6):
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_state(value)

    def test_invalid_scalars_and_costs(self) -> None:
        for value in (0, -1, np.nan, np.inf, 1j, [1], True, "1"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_mu(value)
        for value in (0, -1, np.nan, [1]):
            with self.subTest(a=value), self.assertRaises(ValueError):
                orbital_period(value)
        with self.assertRaises(ValueError):
            delta_v([1,2], [1,2])
        with self.assertRaises(ValueError):
            total_delta_v([1,2,3])

    def test_state_copy_and_numerical_overflow(self) -> None:
        copy = validate_state(self.state)
        copy[0] = 1
        self.assertEqual(self.state[0], self.r)
        with self.assertRaises(RuntimeError):
            specific_energy([1,0,0,1e308,0,0])
