"""Trajectory-keeping impulsive propagation (DT-11), figure helpers, colours and provenance.

References: the validated ``propagate_with_impulses`` (same final state), the
analytic Kepler propagation (RK45 against exact two-body motion, with the
Phase 2 arrival tolerances) and the exact arc extrema of orbital.kepler.
"""
import importlib.util
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from orbital_optimizer import EARTH_MU, OrbitalElements, elements_to_state, hohmann_transfer, propagate_with_impulses
from orbital_optimizer.constants import EARTH_EQUATORIAL_RADIUS, STANDARD_GRAVITY
from orbital_optimizer.orbital.kepler import radius_extrema_on_arc
from orbital_optimizer.orbital.elements import state_to_elements
from orbital_optimizer.provenance import code_version, environment_record
from orbital_optimizer.simulation.impulsive import ImpulsiveTrajectory, propagate_impulsive_arcs
from orbital_optimizer.visualization.style import (
    CATEGORICAL, INK, algorithm_color, algorithm_family, algorithm_label,
)
from orbital_optimizer.visualization.transfer_plots import orbital_plane_basis, plane_coordinates

LEO, GEO = 6_678_000.0, 42_164_000.0


def hohmann_case() -> tuple[np.ndarray, list, float]:
    transfer = hohmann_transfer(LEO, GEO)
    state = np.array([LEO, 0.0, 0.0, 0.0, math.sqrt(EARTH_MU / LEO), 0.0])
    impulses = [(0.0, [0.0, transfer.first_delta_v, 0.0]), (transfer.transfer_time, [0.0, -transfer.second_delta_v, 0.0])]
    return state, impulses, transfer.transfer_time


class ImpulsiveArcsTests(unittest.TestCase):
    def test_final_state_equals_the_validated_wrapper(self) -> None:
        state, impulses, duration = hohmann_case()
        for method in ("RK45", "RK4"):
            with self.subTest(method=method):
                trajectory = propagate_impulsive_arcs(state, impulses, duration, method=method)
                self.assertIsInstance(trajectory, ImpulsiveTrajectory)
                assert_array_equal(trajectory.final_state, propagate_with_impulses(state, impulses, duration, method=method))
        # A burn inside the flight splits it into two arcs; the final state is still the wrapper's.
        mid = [(600.0, [5.0, 10.0, -2.0]), (4000.0, [0.0, -3.0, 1.0])]
        trajectory = propagate_impulsive_arcs(state, mid, 9000.0)
        assert_array_equal(trajectory.final_state, propagate_with_impulses(state, mid, 9000.0))
        self.assertEqual(len(trajectory.arcs), 3)
        assert_array_equal(trajectory.arc_start_times, [0.0, 600.0, 4000.0])

    def test_analytic_arcs_agree_with_numerical_integration(self) -> None:
        state, impulses, duration = hohmann_case()
        analytic = propagate_impulsive_arcs(state, impulses, duration, method="Kepler").final_state
        numerical = propagate_impulsive_arcs(state, impulses, duration, method="RK45").final_state
        # Phase 2 arrival tolerances (tests/test_impulses.py): 1 m and 1e-3 m/s.
        assert_allclose(numerical[:3], analytic[:3], rtol=0, atol=1.0)
        assert_allclose(numerical[3:], analytic[3:], rtol=0, atol=1e-3)
        elements = state_to_elements(analytic)
        self.assertAlmostEqual(elements.semi_major_axis, GEO, delta=1e-3)
        self.assertLess(elements.eccentricity, 1e-9)

    def test_burns_at_the_ends_create_no_empty_arc(self) -> None:
        state, impulses, duration = hohmann_case()
        trajectory = propagate_impulsive_arcs(state, impulses, duration, method="Kepler", samples_per_arc=11)
        self.assertEqual(len(trajectory.arcs), 1)
        assert_array_equal(trajectory.arc_start_times, [0.0])
        self.assertEqual(trajectory.arcs[0].states.shape, (11, 6))
        # The arc starts after the first burn and ends before the last one.
        assert_array_equal(trajectory.arcs[0].states[0, 3:], state[3:] + np.array(impulses[0][1]))
        assert_array_equal(trajectory.final_state[3:], trajectory.arcs[0].states[-1, 3:] + np.array(impulses[1][1]))
        assert_array_equal(trajectory.final_state[:3], trajectory.arcs[0].states[-1, :3])
        self.assertEqual((trajectory.duration, trajectory.mu), (duration, EARTH_MU))
        coast = propagate_impulsive_arcs(state, [], 1000.0, method="Kepler")
        self.assertEqual((len(coast.arcs), coast.burn_times.shape, coast.burn_vectors.shape), (1, (0,), (0, 3)))

    def test_junctions_are_continuous_in_position_and_jump_by_the_burn_in_velocity(self) -> None:
        state, _, _ = hohmann_case()
        burns = [(500.0, [3.0, 20.0, -1.0]), (2500.0, [-4.0, 7.0, 2.0])]
        trajectory = propagate_impulsive_arcs(state, burns, 6000.0, method="Kepler", samples_per_arc=50)
        for index, (_, vector) in enumerate(burns):
            before, after = trajectory.arcs[index].states[-1], trajectory.arcs[index + 1].states[0]
            assert_array_equal(after[:3], before[:3])
            assert_array_equal(after[3:], before[3:] + np.array(vector))
        times = trajectory.times()
        self.assertEqual(times.shape, (150,))
        self.assertEqual((times[0], times[49], times[50], times[-1]), (0.0, 500.0, 500.0, 6000.0))
        self.assertEqual(trajectory.states().shape, (150, 6))
        assert_array_equal(trajectory.burn_times, [500.0, 2500.0])
        with self.assertRaises(ValueError):
            trajectory.burn_vectors[0, 0] = 1.0

    def test_sampled_extrema_lie_inside_the_exact_arc_extrema(self) -> None:
        elements = OrbitalElements(2.0e7, 0.6, 0.3, 1.0, 2.0, 4.0)       # inbound, before periapsis
        state = elements_to_state(elements)
        trajectory = propagate_impulsive_arcs(state, [], 20_000.0, method="Kepler", samples_per_arc=7)
        low, high = trajectory.sampled_radius_extrema()
        exact_low, exact_high = radius_extrema_on_arc(elements.true_anomaly, 20_000.0, 2.0e7, 0.6)
        self.assertAlmostEqual(exact_low, 2.0e7 * 0.4, delta=1e-3)       # the arc contains the periapsis
        # Seven samples miss the periapsis by hundreds of kilometres: sampling is not a bound (RF-4.05).
        self.assertGreater(low - exact_low, 1e5)
        self.assertLessEqual(high, exact_high + 1e-3)

    def test_invalid_inputs(self) -> None:
        state, impulses, duration = hohmann_case()
        for kwargs in ({"method": "Euler"}, {"samples_per_arc": 1}, {"samples_per_arc": 2.0}, {"samples_per_arc": True}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                propagate_impulsive_arcs(state, impulses, duration, **kwargs)
        for bad_impulses in ([(10.0, [0, 1, 0]), (10.0, [0, 1, 0])], [(-1.0, [0, 1, 0])], [(duration + 1, [0, 1, 0])],
                             [(1.0, [0, 1])]):
            with self.subTest(impulses=bad_impulses), self.assertRaises(ValueError):
                propagate_impulsive_arcs(state, bad_impulses, duration)
        with self.assertRaises(ValueError):
            propagate_impulsive_arcs(state, impulses, 0.0)
        with self.assertRaises(ValueError):
            propagate_impulsive_arcs(state, impulses, duration, mu=-1.0)
        # A burn that makes the next arc hyperbolic is outside the validated domain.
        with self.assertRaises(ValueError):
            propagate_impulsive_arcs(state, [(0.0, [0.0, 5000.0, 0.0])], 1000.0, method="Kepler")


class PlaneProjectionTests(unittest.TestCase):
    def test_basis_is_orthonormal_and_spans_the_orbital_plane(self) -> None:
        state = elements_to_state(OrbitalElements(1.2e7, 0.3, 1.1, 0.7, 2.0, 0.9))
        first, second = orbital_plane_basis(state)
        momentum = np.cross(state[:3], state[3:])
        self.assertAlmostEqual(float(first @ first), 1.0, places=14)
        self.assertAlmostEqual(float(second @ second), 1.0, places=14)
        self.assertAlmostEqual(float(first @ second), 0.0, places=14)
        self.assertAlmostEqual(float(first @ momentum) / np.linalg.norm(momentum), 0.0, places=14)
        self.assertAlmostEqual(float(second @ momentum) / np.linalg.norm(momentum), 0.0, places=14)
        assert_allclose(first, state[:3] / np.linalg.norm(state[:3]), atol=1e-15)
        self.assertGreater(float(second @ state[3:]), 0.0)                 # prograde
        # In-plane vectors keep their length; the first sample lies on the first axis.
        positions = np.array([state[:3], 2.5 * state[:3] + 100.0 * state[3:]])
        projected = plane_coordinates(positions, (first, second))
        assert_allclose(np.linalg.norm(projected, axis=1), np.linalg.norm(positions, axis=1), rtol=1e-13)
        assert_allclose(projected[0], [np.linalg.norm(state[:3]), 0.0], atol=1e-6)
        with self.assertRaises(ValueError):
            orbital_plane_basis(np.array([1.0e7, 0.0, 0.0, 100.0, 0.0, 0.0]))     # radial

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib is optional")
    def test_transfer_figure_is_written(self) -> None:
        from orbital_optimizer.visualization.transfer_plots import save_transfer_figure

        state, impulses, duration = hohmann_case()
        trajectory = propagate_impulsive_arcs(state, impulses, duration, method="Kepler", samples_per_arc=181)
        with tempfile.TemporaryDirectory() as folder:
            path = save_transfer_figure(trajectory, LEO, GEO, Path(folder) / "figures" / "transfer.png",
                                        body_radius=EARTH_EQUATORIAL_RADIUS, title="LEO to GEO")
            self.assertGreater(path.stat().st_size, 10_000)


class StyleAndProvenanceTests(unittest.TestCase):
    def test_algorithm_colours_follow_the_family(self) -> None:
        names = {"DE/rand/1/bin@1": "DE", "PSO/gbest/inertia@1": "PSO", "RandomSearch/uniform@1": "RandomSearch",
                 "GA/tournament2/SBX/polynomial/elitist@1": "GA", "NelderMead/Lagarias1998@1": "NelderMead",
                 "Hybrid(DE/rand/1/bin@1 -> NelderMead/Lagarias1998@1)@1": "Hybrid"}
        for name, family in names.items():
            self.assertEqual(algorithm_family(name), family)
        colours = [algorithm_color(name) for name in names]
        self.assertEqual(len(set(colours)), len(colours))                  # one colour per family
        self.assertEqual(algorithm_color("DE/rand/1/bin@1"), CATEGORICAL["light"][0])
        self.assertEqual(algorithm_color("DE/best/2/exp@7"), CATEGORICAL["light"][0])      # same family, same colour
        self.assertEqual(algorithm_color("RandomSearch/uniform@1"), INK["light"]["muted"])  # de-emphasised baseline
        self.assertEqual(algorithm_color("DE/rand/1/bin@1", "dark"), CATEGORICAL["dark"][0])
        self.assertEqual(algorithm_label("PSO/gbest/inertia@1"), "Particle swarm")
        self.assertEqual(algorithm_label("Custom/x@1"), "Custom")
        self.assertEqual(len(set(CATEGORICAL["light"])), 8)
        with self.assertRaises(ValueError):
            algorithm_color("DE/rand/1/bin@1", "sepia")
        with self.assertRaises(ValueError):
            algorithm_family("")

    def test_constants_have_one_definition(self) -> None:
        from orbital_optimizer import constants
        from orbital_optimizer.dynamics import two_body
        from orbital_optimizer.orbital import mechanics

        self.assertIs(two_body.EARTH_MU, constants.EARTH_MU)
        self.assertIs(mechanics.STANDARD_GRAVITY, constants.STANDARD_GRAVITY)
        self.assertEqual((EARTH_MU, STANDARD_GRAVITY, EARTH_EQUATORIAL_RADIUS), (3.986004418e14, 9.80665, 6_378_137.0))

    def test_provenance_records(self) -> None:
        environment = environment_record()
        self.assertEqual(environment["numpy"], np.__version__)
        self.assertEqual(set(environment), {"python", "python_implementation", "numpy", "scipy", "platform", "machine"})
        version = code_version()
        self.assertEqual(set(version), {"package_version", "git_commit", "git_branch", "git_dirty"})
        if version["git_commit"] is not None:
            self.assertRegex(version["git_commit"], r"^[0-9a-f]{40}$")
            self.assertIsInstance(version["git_dirty"], bool)
        with tempfile.TemporaryDirectory() as folder:                      # not a repository
            outside = code_version(folder)
            self.assertEqual((outside["git_commit"], outside["git_branch"], outside["git_dirty"]), (None, None, None))


if __name__ == "__main__":
    unittest.main()
