"""Kepler equation and analytic propagation; float64 roundoff error model.

Analytic results differ from exact ones only by roundoff, with two terms:

1. Evaluation: 1e-12 relative to a (position) or periapsis speed v_p
   (velocity) allows ~30 operations of eps=2.2e-16 amplified by
   dE/dM<=1/(1-e)<=20 over 2*pi of anomaly (~1.3e-13), with ~8x margin.
2. Initial-value conditioning: a=-mu/(2*epsilon) comes from the rounded
   state. delta_a/a <= 4*eps*(v^2/2+mu/r)/|epsilon| <= 4*eps*(3+e)/(1-e),
   which shifts the period and gives a timing error dt=1.5*t*delta_a/a,
   i.e. v_p*dt in position and (mu/r_p^2)*dt in velocity.

Term 2 was missing from the first declared bound; the e=0.95 case exceeded
it (4.5e-8 vs 2.3e-8 m/s at t=T), matching the term-2 prediction.
"""
import unittest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from orbital_optimizer import (
    EARTH_MU, OrbitalElements, eccentricity_vector, elements_to_state, orbital_period,
    propagate_kepler, solve_kepler, specific_angular_momentum, specific_energy,
)
from tests.test_propagator import kepler_reference

EPS = np.finfo(float).eps


def bisect_kepler(mean: float, e: float) -> float:
    """Independent E for M in [-pi,pi]: 80 bisections of the monotone [-pi,pi] bracket."""
    low, high = -np.pi, np.pi
    for _ in range(80):
        mid = (low+high)/2
        if mid-e*np.sin(mid) < mean:
            low = mid
        else:
            high = mid
    return (low+high)/2


def rotation(raan: float, inc: float, arg: float) -> np.ndarray:
    """Rz(raan) Rx(inc) Rz(arg), composed here from elementary rotations."""
    def rz(angle: float) -> np.ndarray:
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[c,-s,0],[s,c,0],[0,0,1]])
    c, s = np.cos(inc), np.sin(inc)
    return rz(raan) @ np.array([[1,0,0],[0,c,-s],[0,s,c]]) @ rz(arg)


class KeplerEquationTests(unittest.TestCase):
    ECCENTRICITIES = (0., 1e-9, .1, .5, .9, .99, 1-1e-6)

    def test_residual_and_bounded_correction(self) -> None:
        means = np.concatenate((np.linspace(-3*np.pi,3*np.pi,61), [0,np.pi,-np.pi,1e-12,-1e-12,40.]))
        for e in self.ECCENTRICITIES:
            for mean in means:
                with self.subTest(e=e, M=mean):
                    E = solve_kepler(mean, e)
                    self.assertLessEqual(abs(E-e*np.sin(E)-mean), 16*EPS*max(1,abs(mean)))
                    # |E-M| = e*|sin E| <= e keeps the revolution count of M.
                    self.assertLessEqual(abs(E-mean), e+16*EPS*max(1,abs(mean)))

    def test_exact_special_values(self) -> None:
        for mean in (-7., 0., 1.3, 50.):
            self.assertEqual(solve_kepler(mean, 0), mean)
        for e in self.ECCENTRICITIES:
            self.assertEqual(solve_kepler(0, e), 0)
            self.assertLessEqual(abs(solve_kepler(np.pi, e)-np.pi), 4*EPS*np.pi)

    def test_agreement_with_independent_bisection(self) -> None:
        for e in self.ECCENTRICITIES[:-1]:
            for mean in np.linspace(-np.pi,np.pi,41):
                with self.subTest(e=e, M=mean):
                    self.assertLess(abs(solve_kepler(mean,e)-bisect_kepler(mean,e)), 1e-12)

    def test_invalid_inputs(self) -> None:
        for e in (-.1, 1, 1.5, np.nan, np.inf, [.1], "0.1", 1j):
            with self.subTest(e=e), self.assertRaises(ValueError):
                solve_kepler(1., e)
        for mean in (np.nan, np.inf, [1,2], "1", 1j):
            with self.subTest(M=mean), self.assertRaises(ValueError):
                solve_kepler(mean, .1)


class KeplerPropagationTests(unittest.TestCase):
    def assert_close_to(self, result, expected: np.ndarray, a: float, e: float) -> None:
        """Apply the module error model at the latest sampled time."""
        periapsis_speed = np.sqrt(EARTH_MU/a*(1+e)/(1-e))
        periapsis_gravity = EARTH_MU/(a*(1-e))**2
        timing = 1.5*result.times[-1]*4*EPS*(3+e)/(1-e)
        self.assertLess(np.max(np.linalg.norm(result.states[:,:3]-expected[:,:3],axis=1)),
                        1e-12*a+periapsis_speed*timing)
        self.assertLess(np.max(np.linalg.norm(result.states[:,3:]-expected[:,3:],axis=1)),
                        1e-12*periapsis_speed+periapsis_gravity*timing)

    def test_planar_ellipse_against_independent_kepler_solution(self) -> None:
        a, e = 18e6, .6
        period = orbital_period(a)
        times = np.linspace(0,period,121)
        initial = [a*(1-e),0,0,0,np.sqrt(EARTH_MU*(1+e)/(a*(1-e))),0]
        result = propagate_kepler(initial,period,sample_times=times)
        self.assert_close_to(result, kepler_reference(times,a,e), a, e)

    def test_three_dimensional_and_degenerate_geometry(self) -> None:
        cases = [(9e6,.2,.7,.4,.9), (30e6,.95,1.2,2.,4.), (12e6,.4,np.pi,0.,1.),
                 (12e6,.4,0.,0.,2.5), (8e6,1e-9,.5,1.,0.)]
        for a, e, inc, raan, arg in cases:
            with self.subTest(a=a, e=e, i=inc):
                period = orbital_period(a)
                times = np.linspace(0,period,201)
                matrix = rotation(raan,inc,arg)
                planar = kepler_reference(times,a,e)
                expected = np.column_stack((planar[:,:3]@matrix.T, planar[:,3:]@matrix.T))
                result = propagate_kepler(elements_to_state(OrbitalElements(a,e,inc,raan,arg)),period,sample_times=times)
                self.assert_close_to(result, expected, a, e)

    def test_circular_orbits_against_rotation_solution(self) -> None:
        for inc in (0., .5, np.pi):
            with self.subTest(i=inc):
                initial = elements_to_state(OrbitalElements(7e6,0.,inc,.3,0.,1.))
                period = orbital_period(7e6)
                times = np.linspace(0,period,101)
                n = 2*np.pi/period
                phase = (n*times)[:,None]
                expected = np.column_stack((initial[:3]*np.cos(phase)+initial[3:]*np.sin(phase)/n,
                                            -n*initial[:3]*np.sin(phase)+initial[3:]*np.cos(phase)))
                self.assert_close_to(propagate_kepler(initial,period,sample_times=times), expected, 7e6, 0.)

    def test_invariants_and_multi_revolution_closure(self) -> None:
        initial = elements_to_state(OrbitalElements(9e6,.2,.7,.4,.9,1.2))
        period = orbital_period(9e6)
        result = propagate_kepler(initial,10*period,sample_times=np.linspace(0,10*period,1001))
        energy = np.array([specific_energy(s) for s in result.states])
        h = np.array([specific_angular_momentum(s) for s in result.states])
        ecc = np.array([eccentricity_vector(s) for s in result.states])
        self.assertLess(np.max(np.abs(energy-energy[0]))/abs(energy[0]),1e-12)
        self.assertLess(np.max(np.linalg.norm(h-h[0],axis=1))/np.linalg.norm(h[0]),1e-12)
        self.assertLess(np.max(np.linalg.norm(ecc-ecc[0],axis=1)),1e-12)
        for turns in (1,10,100):
            with self.subTest(turns=turns):
                closed = propagate_kepler(initial,turns*period)
                self.assert_close_to(closed, np.array([initial,initial]), 9e6, .2)

    def test_result_contract(self) -> None:
        initial = elements_to_state(OrbitalElements(9e6,.2,.7,.4,.9,1.2))
        result = propagate_kepler(initial,100.)
        self.assertEqual((result.method,result.nfev,result.mu),("Kepler",0,EARTH_MU))
        assert_array_equal(result.times,[0,100.])
        assert_array_equal(result.states[0],initial)
        with self.assertRaises(ValueError):
            result.states[0,0] = 0
        sampled = propagate_kepler(initial,10.,sample_times=[2.,5.])
        assert_allclose(sampled.times,[0,2,5,10],rtol=0,atol=0)
        unit = propagate_kepler([1,0,0,0,1,0],np.pi,mu=1)
        assert_allclose(unit.final_state,[-1,0,0,0,-1,0],atol=1e-14)

    def test_invalid_inputs(self) -> None:
        state = [7e6,0,0,0,np.sqrt(EARTH_MU/7e6),0]
        for duration in (0,-1,np.nan,np.inf,1j,True):
            with self.subTest(duration=duration), self.assertRaises(ValueError):
                propagate_kepler(state,duration)
        for options in ({"mu":-1},{"sample_times":[]},{"sample_times":[2,1]},{"sample_times":[0,11]}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                propagate_kepler(state,10,**options)
        for bad in ([7e6,0,0,0,0,0],[7e6,0,0,100,0,0],[7e6,0,0,0,20000,0],[0]*6,[1]*5):
            with self.subTest(state=bad), self.assertRaises(ValueError):
                propagate_kepler(bad,10)
