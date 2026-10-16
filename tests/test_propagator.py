"""Independent analytic references; acceptance bounds declared in docs/phase1.md."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose

from orbital_optimizer import (
    EARTH_MU, OrbitalElements, elements_to_state, propagate, orbital_period,
    specific_energy, specific_angular_momentum, eccentricity_vector, propagate_kepler,
)


def kepler_reference(times: np.ndarray, a: float, e: float) -> np.ndarray:
    """SI planar ellipse starting at periapsis; solve E-e*sin(E)=n*t by bisection.

    Independent of production conversions and ODE solver. 60 bisections reduce
    the [0,2*pi] bracket below float64 resolution. Test reference only.
    """
    n = np.sqrt(EARTH_MU/a**3)
    rows = []
    for t in times:
        mean = (n*t) % (2*np.pi)
        low, high = 0., 2*np.pi
        for _ in range(60):
            mid = (low+high)/2
            if mid-e*np.sin(mid) < mean:
                low = mid
            else:
                high = mid
        E = (low+high)/2
        q = np.sqrt(1-e**2)
        rate = n/(1-e*np.cos(E))
        rows.append([a*(np.cos(E)-e), a*q*np.sin(E), 0,
                     -a*np.sin(E)*rate, a*q*np.cos(E)*rate, 0])
    return np.array(rows)


class PropagatorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.r = 7e6
        self.v = np.sqrt(EARTH_MU/self.r)
        self.state = np.array([self.r,0,0,0,self.v,0])
        self.period = 2*np.pi*self.r/self.v

    def test_circular_analytic_motion_and_period_closure(self) -> None:
        times = np.linspace(0,self.period,101)
        for method in ("RK4", "RK45"):
            with self.subTest(method=method):
                result = propagate(self.state,self.period,method=method,sample_times=times,
                                   max_step=self.period/1000 if method=="RK4" else 60)
                phase = self.v/self.r*times
                expected_r = np.column_stack([self.r*np.cos(phase),self.r*np.sin(phase),np.zeros_like(phase)])
                expected_v = np.column_stack([-self.v*np.sin(phase),self.v*np.cos(phase),np.zeros_like(phase)])
                self.assertLess(np.max(np.linalg.norm(result.states[:,:3]-expected_r,axis=1)), 1.)
                self.assertLess(np.max(np.linalg.norm(result.states[:,3:]-expected_v,axis=1)), .001)
                self.assertLess(np.max(np.abs(np.linalg.norm(result.states[:,:3],axis=1)-self.r)), 1.)
                self.assertLess(np.max(np.abs(np.linalg.norm(result.states[:,3:],axis=1)-self.v)), .001)
                self.assertLess(np.linalg.norm(result.final_state[:3]-self.state[:3]), 1.)

    def test_eccentric_motion_against_independent_kepler_solution(self) -> None:
        a,e = 18e6,.6
        period = 2*np.pi*np.sqrt(a**3/EARTH_MU)
        times = np.linspace(0,period,121)
        initial = [a*(1-e),0,0,0,np.sqrt(EARTH_MU*(1+e)/(a*(1-e))),0]
        result = propagate(initial,period,sample_times=times,rtol=1e-12,atol=[1e-5]*3+[1e-8]*3)
        expected = kepler_reference(times,a,e)
        self.assertLess(np.max(np.linalg.norm(result.states[:,:3]-expected[:,:3],axis=1)),1.)
        self.assertLess(np.max(np.linalg.norm(result.states[:,3:]-expected[:,3:],axis=1)),.001)
        self.assertLess(abs(np.linalg.norm(result.states[60,:3])-a*(1+e)),1.)

    def test_ten_period_invariants_in_inclined_ellipse(self) -> None:
        initial = elements_to_state(OrbitalElements(9e6,.2,.7,.4,.9,1.2))
        result = propagate(initial,10*orbital_period(9e6))
        energy = np.array([specific_energy(s) for s in result.states])
        h = np.array([specific_angular_momentum(s) for s in result.states])
        e = np.array([eccentricity_vector(s) for s in result.states])
        self.assertLess(np.max(np.abs(energy-energy[0]))/abs(energy[0]),1e-7)
        self.assertLess(np.max(np.linalg.norm(h-h[0],axis=1))/np.linalg.norm(h[0]),1e-7)
        self.assertLess(np.max(np.linalg.norm(e-e[0],axis=1)),1e-7)
        # Vis-viva with the input a: v^2-mu(2/r-1/a) = 2*delta_epsilon, and
        # 2|epsilon|/v^2 <= (1+e)/(1-e) = 1.5, so the energy bound gives 1.5e-7.
        radius = np.linalg.norm(result.states[:,:3],axis=1)
        speed2 = np.sum(result.states[:,3:]**2,axis=1)
        self.assertLess(np.max(np.abs(speed2/(EARTH_MU*(2/radius-1/9e6))-1)),1.5e-7)

    def test_propagated_apsides_in_inclined_ellipse(self) -> None:
        a, e, nu = 9e6, .2, 1.2
        initial = elements_to_state(OrbitalElements(a,e,.7,.4,.9,nu))
        period = orbital_period(a)
        # Closed-form passage times: tan(E/2)=sqrt((1-e)/(1+e))tan(nu/2), M=E-e sin E.
        E0 = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(nu/2))
        periapsis_time = (2*np.pi-(E0-e*np.sin(E0)))/(2*np.pi)*period
        apoapsis_time = periapsis_time-period/2
        result = propagate(initial,period,sample_times=[apoapsis_time,periapsis_time])
        for index, expected in ((1,a*(1+e)),(2,a*(1-e))):
            r, v = result.states[index,:3], result.states[index,3:]
            self.assertLess(abs(np.linalg.norm(r)-expected),1.)
            # Flight-path sine; 1 m and 1e-3 m/s errors perturb it by ~1e-7.
            self.assertLess(abs(np.dot(r,v))/(np.linalg.norm(r)*np.linalg.norm(v)),1e-6)

    def test_numerical_against_analytic_kepler_propagation(self) -> None:
        initial = elements_to_state(OrbitalElements(9e6,.2,.7,.4,.9,1.2))
        period = orbital_period(9e6)
        times = np.linspace(0,period,201)
        expected = propagate_kepler(initial,period,sample_times=times).states
        for method, step in (("RK45",60.),("RK4",period/1000)):
            with self.subTest(method=method):
                result = propagate(initial,period,method=method,sample_times=times,max_step=step)
                self.assertLess(np.max(np.linalg.norm(result.states[:,:3]-expected[:,:3],axis=1)),1.)
                self.assertLess(np.max(np.linalg.norm(result.states[:,3:]-expected[:,3:],axis=1)),.001)

    def test_rk4_fourth_order_convergence(self) -> None:
        errors = []
        for count in (100,200,400):
            result = propagate(self.state,self.period,method="RK4",max_step=self.period/count)
            errors.append(np.linalg.norm(result.final_state[:3]-self.state[:3]))
        orders = np.log2(np.array(errors[:-1])/errors[1:])
        self.assertTrue(np.all((orders>3.7)&(orders<4.5)), (errors,orders))

    def test_agreement_between_integrators(self) -> None:
        # Samples fall on every 5th nominal RK4 step; float rounding still adds
        # some shortened sub-steps (nfev ~4.4e3, not 4e3), which only reduce error.
        ellipse = elements_to_state(OrbitalElements(9e6,.2,.7,.4,.9,1.2))
        for name, initial, period in (("circular",self.state,self.period),
                                      ("inclined ellipse",ellipse,orbital_period(9e6))):
            with self.subTest(orbit=name):
                times = np.linspace(0,period,201)
                rk4 = propagate(initial,period,method="RK4",max_step=period/1000,sample_times=times)
                rk45 = propagate(initial,period,rtol=1e-12,atol=1e-9,sample_times=times)
                self.assertLess(np.max(np.linalg.norm(rk4.states[:,:3]-rk45.states[:,:3],axis=1)),1.)
                self.assertLess(np.max(np.linalg.norm(rk4.states[:,3:]-rk45.states[:,3:],axis=1)),.001)

    def test_sampling_copies_and_small_interval(self) -> None:
        for method in ("RK4","RK45"):
            result = propagate(self.state,1.,method=method,sample_times=[.2,.7])
            assert_allclose(result.times,[0,.2,.7,1.],atol=0,rtol=0)
            with self.assertRaises(ValueError):
                result.states[0,0] = 0
            with self.assertRaises(ValueError):
                result.times[0] = 1
            final = result.final_state
            final[:] = 0
            self.assertNotEqual(result.final_state[0],0)
            self.assertEqual(result.times[-1],1.)
        tiny = propagate(self.state,1e-6)
        self.assertEqual(tiny.times[-1],1e-6)
        assert_allclose(self.state,[self.r,0,0,0,self.v,0])

    def test_custom_mu_propagation(self) -> None:
        result = propagate([1,0,0,0,1,0],2*np.pi,mu=1,max_step=.05,atol=1e-12)
        assert_allclose(result.final_state,[1,0,0,0,1,0],atol=1e-8)

    def test_invalid_inputs(self) -> None:
        for duration in (0,-1,np.nan,np.inf,1j,True):
            with self.subTest(duration=duration), self.assertRaises(ValueError):
                propagate(self.state,duration)
        for options in ({"max_step":0},{"mu":-1},{"method":"Euler"},{"rtol":0},
                        {"rtol":1e-20},{"atol":0},{"atol":[1,1]}, {"atol":[1,1,1,1,1,np.nan]},
                        {"sample_times":[]},{"sample_times":[1,1]}, {"sample_times":[2,1]},
                        {"sample_times":[-1,1]}, {"sample_times":[0,11]}):
            with self.subTest(options=options),self.assertRaises(ValueError):
                propagate(self.state,10,**options)

    def test_solver_failure_is_explicit(self) -> None:
        path = "orbital_optimizer.simulation.propagator.solve_ivp"
        with patch(path,return_value=SimpleNamespace(success=False,message="test failure")):
            with self.assertRaisesRegex(RuntimeError,"test failure"):
                propagate(self.state,10)
