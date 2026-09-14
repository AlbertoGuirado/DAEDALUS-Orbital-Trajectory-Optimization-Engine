"""Element geometry, degeneracies and independent vis-viva checks."""
import unittest
import numpy as np
from numpy.testing import assert_allclose

from orbital_optimizer import EARTH_MU, OrbitalElements, elements_to_state, state_to_elements


class ElementsTests(unittest.TestCase):
    def test_known_periapsis_and_apsides(self) -> None:
        a, e = 12e6, 0.3
        for nu, radius in [(0, a*(1-e)), (np.pi, a*(1+e))]:
            state = elements_to_state(OrbitalElements(a,e,0,true_anomaly=nu))
            self.assertAlmostEqual(np.linalg.norm(state[:3])/radius, 1., places=14)
            self.assertAlmostEqual(np.dot(state[3:], state[3:])/(EARTH_MU*(2/radius-1/a)), 1., places=14)
        state = elements_to_state(OrbitalElements(a,e,np.pi/2,raan=np.pi/2))
        assert_allclose(state[:3], [0,a*(1-e),0], atol=1e-8)
        assert_allclose(state[3:], [0,0,np.sqrt(EARTH_MU*(1+e)/(a*(1-e)))], atol=1e-10)

    def test_roundtrip_general_and_degenerate_geometry(self) -> None:
        for e in (0., 1e-9, 0.2, 0.8):
            for i in (0., 1e-9, 0.6, np.pi/2, np.pi-1e-9, np.pi):
                with self.subTest(e=e, i=i):
                    original = OrbitalElements(20e6,e,i,1.1,0.8,4.3)
                    state = elements_to_state(original)
                    recovered = state_to_elements(state)
                    rebuilt = elements_to_state(recovered)
                    assert_allclose(rebuilt[:3], state[:3], rtol=1e-12, atol=1e-7)
                    assert_allclose(rebuilt[3:], state[3:], rtol=1e-12, atol=1e-10)
                    if e == 0:
                        self.assertEqual(recovered.argument_of_periapsis, 0)
                    if i in (0, np.pi):
                        self.assertEqual(recovered.raan, 0)
        original = OrbitalElements(12e6,0.3,0.7,1.1,0.8,4.3)
        recovered = state_to_elements(elements_to_state(original))
        assert_allclose(list(vars(recovered).values()), list(vars(original).values()), rtol=1e-12, atol=1e-12)

    def test_invalid_elements_and_unsupported_states(self) -> None:
        for args in ((-1,0,0), (1,-.1,0), (1,1,0), (1,.1,-1), (1,.1,4), (np.nan,0,0)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                OrbitalElements(*args)
        for state in ([7e6,0,0,0,0,0], [7e6,0,0,100,0,0],
                      [7e6,0,0,0,np.sqrt(2*EARTH_MU/7e6),0], [7e6,0,0,0,20000,0]):
            with self.subTest(state=state), self.assertRaises(ValueError):
                state_to_elements(state)
