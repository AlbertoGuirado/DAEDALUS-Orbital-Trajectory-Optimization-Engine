"""Comparison optimizers (Phase 3 increment d): contracts shared with DE, operators and draw order.

Scope: RandomSearch, ParticleSwarm, GeneticAlgorithm, NelderMead and
HybridOptimizer. Differential Evolution has its own module.

Performance thresholds were declared before the first run and derive from the
null baseline, not from observed results:
  - Sphere: the best of N independent uniform samples in a box of volume V
    has f ~ eps_RS = (V / (N V_D))^(2/D), V_D the volume of the unit D-ball
    (the radius at which one sample is expected inside the ball). A directed
    search must reach 1e-2 * eps_RS on every guard seed.
  - Every function: the median final gap of a global method over the guard
    seeds must be lower than the median of RandomSearch with the same seeds
    and budget.
They detect a broken optimizer; they are not statistical performance evidence
(that is the multi-seed benchmark).
"""
import json
import math
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from orbital_optimizer.optimization import (
    Bounds, DifferentialEvolution, Evaluation, EvaluationStatus, FunctionProblem, GeneticAlgorithm,
    GeneticAlgorithmSettings, HybridOptimizer, HybridSettings, NelderMead, NelderMeadSettings,
    OptimizationResult, Optimizer, ParticleSwarm, ParticleSwarmSettings, RandomSearch,
    RandomSearchSettings, TerminationCriteria, TerminationReason, rastrigin_benchmark,
    rosenbrock_benchmark, sphere, sphere_benchmark,
)
from orbital_optimizer.optimization.genetic_algorithm import _polynomial_step, _sbx_spread
from orbital_optimizer.optimization.optimizer import signed_unit_difference, wrap_unit

GUARD_SEEDS = (0, 1, 2, 3, 4)


def unit_bounds(dimension: int, periodic=None) -> Bounds:
    return Bounds(np.zeros(dimension), np.ones(dimension),
                  periodic=None if periodic is None else np.array(periodic, dtype=bool))


class RecordingProblem:
    """Delegate to FunctionProblem and keep a copy of every evaluated x."""

    def __init__(self, function, bounds: Bounds) -> None:
        self.inner = FunctionProblem(function, bounds)
        self.bounds = bounds
        self.name = "recording"
        self.points: list[np.ndarray] = []

    def evaluate(self, x: np.ndarray) -> Evaluation:
        self.points.append(np.array(x, copy=True))
        return self.inner.evaluate(x)


def global_optimizers() -> list:
    return [RandomSearch(), ParticleSwarm(), GeneticAlgorithm(), HybridOptimizer()]


def all_optimizers() -> list:
    return global_optimizers() + [NelderMead()]


def comparable(result: OptimizationResult) -> dict:
    data = result.to_dict()
    data.pop("wall_time_s")
    return data


def partly_invalid(x: np.ndarray) -> float:
    if x[0] < 0.2:
        return -math.inf
    if x[0] > 0.9:
        return math.nan
    return float(np.sum((x - 0.5) ** 2))


class SharedContractTests(unittest.TestCase):
    """RF-3.10, RF-3.16, RF-3.17, RF-3.18 for every additional optimizer."""

    def setUp(self) -> None:
        self.benchmark = rastrigin_benchmark(3, lower=-2.56, upper=5.12)

    def test_protocol_and_names(self) -> None:
        names = set()
        for optimizer in all_optimizers() + [DifferentialEvolution()]:
            with self.subTest(optimizer=optimizer.name):
                self.assertIsInstance(optimizer, Optimizer)
                self.assertRegex(optimizer.name, r"@\d+$")
                names.add(optimizer.name)
        self.assertEqual(len(names), 6)

    def test_budget_is_respected_exactly(self) -> None:
        for optimizer in all_optimizers():
            for budget in (30, 31, 97, 400):
                with self.subTest(optimizer=optimizer.name, budget=budget):
                    result = optimizer.minimize(self.benchmark.problem(), TerminationCriteria(budget), seed=4)
                    self.assertEqual(result.nfev, budget)
                    self.assertIs(result.termination_reason, TerminationReason.MAX_EVALUATIONS)
                    self.assertEqual(len(result.history.best_so_far), budget)
                    self.assertEqual(len(result.history.status), budget)
                    self.assertEqual(result.history.generation_nfev[-1], budget)
                    self.assertEqual(sum(result.status_counts.values()), budget)
                    self.assertTrue(np.all(np.diff(result.history.generation_nfev) > 0))
                    self.assertEqual(result.generations, len(result.history.generation_nfev) - 1)

    def test_same_seed_is_bitwise_identical_and_seeds_differ(self) -> None:
        for optimizer in all_optimizers():
            with self.subTest(optimizer=optimizer.name):
                first = optimizer.minimize(self.benchmark.problem(), TerminationCriteria(600), seed=123)
                second = optimizer.minimize(self.benchmark.problem(), TerminationCriteria(600), seed=123)
                self.assertEqual(comparable(first), comparable(second))
                self.assertEqual(first.x_best.tobytes(), second.x_best.tobytes())
                self.assertEqual(first.history.best_so_far.tobytes(), second.history.best_so_far.tobytes())
                other = optimizer.minimize(self.benchmark.problem(), TerminationCriteria(600), seed=124)
                self.assertNotEqual(first.history.best_so_far.tobytes(), other.history.best_so_far.tobytes())

    def test_result_is_consistent_with_the_evaluated_points(self) -> None:
        bounds = Bounds([-2.0, -1.0], [3.0, 4.0], names=("a", "b"), units=("m", "s"))
        for optimizer in all_optimizers():
            with self.subTest(optimizer=optimizer.name):
                problem = RecordingProblem(sphere, bounds)
                result = optimizer.minimize(problem, TerminationCriteria(300), seed=8)
                points = np.array(problem.points)
                values = np.array([sphere(p) for p in points])
                first_best = int(np.argmin(values))
                self.assertEqual(result.f_best, values[first_best])
                assert_array_equal(result.x_best, points[first_best])
                assert_array_equal(result.x_best, bounds.from_unit(result.u_best))
                assert_array_equal(result.history.best_so_far, np.minimum.accumulate(values))
                self.assertTrue(all(bounds.contains(p) for p in points))
                self.assertEqual(result.algorithm, optimizer.name)
                self.assertEqual(result.seed, 8)
                self.assertEqual(result.problem_info["bounds"]["units"], ["m", "s"])
                json.dumps(result.to_dict(), allow_nan=False)
                with self.assertRaises(ValueError):
                    result.x_best[0] = 0.0

    def test_invalid_evaluations_never_become_best(self) -> None:
        for optimizer in all_optimizers():
            with self.subTest(optimizer=optimizer.name):
                problem = RecordingProblem(partly_invalid, unit_bounds(2))
                result = optimizer.minimize(problem, TerminationCriteria(600), seed=5)
                points = np.array(problem.points)
                expected_invalid = int(np.count_nonzero((points[:, 0] < 0.2) | (points[:, 0] > 0.9)))
                self.assertEqual(result.status_counts["non_finite"], expected_invalid)
                self.assertEqual(result.n_invalid, expected_invalid)
                self.assertTrue(result.found_valid)
                self.assertTrue(0.2 <= result.x_best[0] <= 0.9)
                finite = [v for v in (partly_invalid(p) for p in points) if math.isfinite(v)]
                self.assertEqual(result.f_best, min(finite))

    def test_all_invalid_problem(self) -> None:
        for optimizer in all_optimizers():
            with self.subTest(optimizer=optimizer.name):
                result = optimizer.minimize(FunctionProblem(lambda x: math.nan, unit_bounds(2)),
                                            TerminationCriteria(75), seed=1)
                self.assertFalse(result.found_valid)
                self.assertIsNone(result.x_best)
                self.assertIsNone(result.best_evaluation)
                self.assertEqual(result.f_best, math.inf)
                self.assertEqual(result.nfev, 75)
                self.assertEqual(result.status_counts["non_finite"], 75)
                json.dumps(result.to_dict(), allow_nan=False)

    def test_other_termination_criteria(self) -> None:
        problem = FunctionProblem(sphere, Bounds([-2.0, -1.0], [3.0, 4.0]))
        for optimizer in (ParticleSwarm(), GeneticAlgorithm(), NelderMead(), HybridOptimizer()):
            with self.subTest(optimizer=optimizer.name, criterion="target"):
                result = optimizer.minimize(problem, TerminationCriteria(50_000, target_objective=1e-3), seed=0)
                self.assertIs(result.termination_reason, TerminationReason.TARGET_REACHED)
                self.assertLessEqual(result.f_best, 1e-3)
                self.assertLess(result.nfev, 50_000)
        for optimizer in (RandomSearch(), ParticleSwarm(), GeneticAlgorithm(), NelderMead()):
            with self.subTest(optimizer=optimizer.name, criterion="max_generations"):
                result = optimizer.minimize(problem, TerminationCriteria(50_000, max_generations=3), seed=0)
                self.assertIs(result.termination_reason, TerminationReason.MAX_GENERATIONS)
                self.assertEqual(result.generations, 3)
            with self.subTest(optimizer=optimizer.name, criterion="stagnation"):
                constant = FunctionProblem(lambda x: 1.0, unit_bounds(2))
                result = optimizer.minimize(constant, TerminationCriteria(5_000, stagnation_generations=3), seed=0)
                self.assertIs(result.termination_reason, TerminationReason.STAGNATION)
                self.assertEqual(result.generations, 3)
        result = NelderMead().minimize(problem, TerminationCriteria(50_000, x_tol=1e-6), seed=0)
        self.assertIs(result.termination_reason, TerminationReason.POPULATION_CONVERGED)
        self.assertLessEqual(result.history.generation_spread[-1], 1e-6)

    def test_invalid_arguments_are_rejected(self) -> None:
        problem = FunctionProblem(sphere, unit_bounds(2))
        for optimizer in all_optimizers():
            for seed in (-1, True, 1.5, None):
                with self.subTest(optimizer=optimizer.name, seed=seed), self.assertRaises(ValueError):
                    optimizer.minimize(problem, TerminationCriteria(100), seed=seed)  # type: ignore[arg-type]
            with self.subTest(optimizer=optimizer.name), self.assertRaises(TypeError):
                optimizer.minimize(object(), TerminationCriteria(100), seed=0)  # type: ignore[arg-type]
            with self.subTest(optimizer=optimizer.name), self.assertRaises(ValueError):
                optimizer.minimize(problem, 100, seed=0)  # type: ignore[arg-type]
        for optimizer, minimum in ((ParticleSwarm(), 20), (GeneticAlgorithm(), 20), (NelderMead(), 3),
                                   (HybridOptimizer(), 20)):
            with self.subTest(optimizer=optimizer.name), self.assertRaises(ValueError):
                optimizer.minimize(problem, TerminationCriteria(minimum - 1), seed=0)
        self.assertEqual(RandomSearch().minimize(problem, TerminationCriteria(1), seed=0).nfev, 1)


class KnownOptimaTests(unittest.TestCase):
    """RV-3.02 for the additional global optimizers, against the null baseline."""

    @staticmethod
    def gaps(optimizer, benchmark, budget: int) -> list[float]:
        values = []
        for seed in GUARD_SEEDS:
            result = optimizer.minimize(benchmark.problem(), TerminationCriteria(budget), seed=seed)
            assert result.found_valid and result.f_best >= benchmark.minimum
            assert benchmark.bounds.contains(result.x_best)
            values.append(result.f_best - benchmark.minimum)
        return values

    def test_sphere_reaches_a_fraction_of_the_random_sampling_scale(self) -> None:
        dimension, lower, upper, budget = 3, -2.56, 5.12, 3_000
        unit_ball = math.pi ** (dimension / 2) / math.gamma(dimension / 2 + 1)
        random_scale = ((upper - lower) ** dimension / (budget * unit_ball)) ** (2 / dimension)
        self.assertAlmostEqual(random_scale, 0.109, places=3)
        benchmark = sphere_benchmark(dimension, lower, upper)
        for optimizer in (ParticleSwarm(), GeneticAlgorithm(), HybridOptimizer()):
            with self.subTest(optimizer=optimizer.name):
                gaps = self.gaps(optimizer, benchmark, budget)
                self.assertLessEqual(max(gaps), 1e-2 * random_scale, gaps)

    def test_global_methods_beat_random_search_medians(self) -> None:
        cases = ((sphere_benchmark(3, -2.56, 5.12), 3_000), (rosenbrock_benchmark(2, -2.048, 4.096), 4_000),
                 (rastrigin_benchmark(3, -2.56, 5.12), 6_000))
        for benchmark, budget in cases:
            baseline = float(np.median(self.gaps(RandomSearch(), benchmark, budget)))
            self.assertGreater(baseline, 0.0)
            for optimizer in (DifferentialEvolution(), ParticleSwarm(), GeneticAlgorithm(), HybridOptimizer()):
                with self.subTest(function=benchmark.name, optimizer=optimizer.name):
                    self.assertLess(float(np.median(self.gaps(optimizer, benchmark, budget))), baseline)


class RandomSearchTests(unittest.TestCase):
    def test_samples_follow_the_documented_draws_and_prefix_property(self) -> None:
        bounds = Bounds([-2.0, 0.0, 1.0], [3.0, 1.0, 5.0])
        short, long = RecordingProblem(sphere, bounds), RecordingProblem(sphere, bounds)
        RandomSearch(RandomSearchSettings(batch_size=7)).minimize(short, TerminationCriteria(20), seed=9)
        result = RandomSearch(RandomSearchSettings(batch_size=7)).minimize(long, TerminationCriteria(53), seed=9)
        rng = np.random.Generator(np.random.PCG64(9))
        expected = np.concatenate([rng.random((7, 3)) for _ in range(8)])[:53]
        assert_array_equal(np.array(long.points), bounds.from_unit(expected))
        assert_array_equal(np.array(short.points), np.array(long.points)[:20])
        assert_array_equal(result.history.generation_nfev, [7, 14, 21, 28, 35, 42, 49, 53])
        self.assertEqual(result.settings["batch_size"], 7)

    def test_settings_validation(self) -> None:
        for kwargs in ({"batch_size": 0}, {"batch_size": True}, {"batch_size": 2.0}, {"record_population": 1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                RandomSearchSettings(**kwargs)
        with self.assertRaises(ValueError):
            RandomSearch("settings")  # type: ignore[arg-type]
        self.assertEqual(RandomSearchSettings().resolved_batch_size(3), 30)


def replicate_pso(seed: int, size: int, dimension: int, inertia: float, cognitive: float, social: float,
                  limit: float, periodic: list[bool], function, generations: int) -> tuple[list[np.ndarray], int]:
    """Independent scalar re-implementation of the documented PSO draws and update."""
    rng = np.random.Generator(np.random.PCG64(seed))
    x = rng.random((size, dimension))
    v = rng.random((size, dimension)) - 0.5
    best, best_values = x.copy(), [function(row) for row in x]
    batches, walls = [x.copy()], 0
    for _ in range(generations):
        leader = best_values.index(min(best_values))
        r1 = rng.random((size, dimension))
        r2 = rng.random((size, dimension))
        moved = x.copy()
        for i in range(size):
            for j in range(dimension):
                to_personal = best[i, j] - x[i, j]
                to_leader = best[leader, j] - x[i, j]
                if periodic[j]:
                    to_personal = (to_personal + 0.5) % 1.0 - 0.5
                    to_leader = (to_leader + 0.5) % 1.0 - 0.5
                speed = inertia * v[i, j] + cognitive * r1[i, j] * to_personal + social * r2[i, j] * to_leader
                speed = max(-limit, min(limit, speed))
                position = x[i, j] + speed
                if periodic[j]:
                    position = position % 1.0
                    position = 0.0 if position >= 1.0 else position
                elif position < 0.0:
                    position, speed, walls = 0.0, 0.0, walls + 1
                elif position > 1.0:
                    position, speed, walls = 1.0, 0.0, walls + 1
                v[i, j], moved[i, j] = speed, position
        x = moved
        batches.append(x.copy())
        for i in range(size):
            value = function(x[i])
            if value <= best_values[i]:
                best_values[i], best[i] = value, x[i]
    return batches, walls


class ParticleSwarmTests(unittest.TestCase):
    def test_evaluated_points_match_independent_replication(self) -> None:
        size, dimension, periodic = 6, 3, [False, True, False]
        settings = ParticleSwarmSettings(swarm_size=size, inertia=0.9, cognitive=2.0, social=2.0, max_velocity=0.8)
        # The optimum sits on a wall and near the seam, so walls and wrapping are exercised.
        def function(u: np.ndarray) -> float:
            return float(u[0] + (1.0 - math.cos(2 * math.pi * (u[1] - 0.97))) + (1.0 - u[2]))

        problem = RecordingProblem(function, unit_bounds(dimension, periodic))
        result = ParticleSwarm(settings).minimize(problem, TerminationCriteria(size * 6), seed=11)
        batches, walls = replicate_pso(11, size, dimension, 0.9, 2.0, 2.0, 0.8, periodic, function, 5)
        self.assertGreater(walls, 0)
        assert_array_equal(np.array(problem.points), np.concatenate(batches))
        self.assertEqual(result.generations, 5)

    def test_periodic_helpers(self) -> None:
        periodic = np.array([True, False])
        assert_allclose(signed_unit_difference(np.array([0.05, 0.05]), np.array([0.95, 0.95]), periodic),
                        [0.1, -0.9], atol=1e-15)
        assert_allclose(signed_unit_difference(np.array([0.95, 0.2]), np.array([0.05, 0.7]), periodic),
                        [-0.1, -0.5], atol=1e-15)
        wrapped = wrap_unit(np.array([[1.25, 1.25], [-0.25, -0.25], [-1e-20, 0.5]]), periodic)
        assert_array_equal(wrapped, [[0.25, 1.25], [0.75, -0.25], [0.0, 0.5]])

    def test_personal_bests_are_recorded_and_bounded(self) -> None:
        result = ParticleSwarm(ParticleSwarmSettings(record_population=True)).minimize(
            FunctionProblem(sphere, Bounds([-2.0, -1.0], [3.0, 4.0])), TerminationCriteria(400), seed=3)
        populations = result.history.populations
        self.assertEqual(populations.shape, (result.generations + 1, 20, 2))
        self.assertTrue(np.all((populations >= 0.0) & (populations <= 1.0)))
        self.assertTrue(np.all(np.diff(result.history.generation_best) <= 0))

    def test_settings_validation(self) -> None:
        for kwargs in ({"swarm_size": 1}, {"swarm_size": True}, {"swarm_size": 5.0}, {"inertia": -0.1},
                       {"inertia": 1.1}, {"inertia": True}, {"cognitive": -1.0}, {"cognitive": 0.0, "social": 0.0},
                       {"social": math.inf}, {"max_velocity": 0.0}, {"max_velocity": 1.5}, {"record_population": 1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                ParticleSwarmSettings(**kwargs)
        self.assertEqual(ParticleSwarmSettings().resolved_swarm_size(3), 30)


def replicate_ga(seed: int, size: int, dimension: int, p_cross: float, eta_cross: float, p_mutation: float,
                 eta_mutation: float, periodic: list[bool], function, generations: int) -> list[np.ndarray]:
    """Independent scalar re-implementation of the documented GA draws and operators."""
    rng = np.random.Generator(np.random.PCG64(seed))
    population = rng.random((size, dimension))
    values = [function(row) for row in population]
    batches = [population.copy()]
    for _ in range(generations):
        tournaments = rng.integers(0, size, size=(size, 2))
        recombine = rng.random(size // 2)
        take_part = rng.random((size // 2, dimension))
        spread_draws = rng.random((size // 2, dimension))
        mutate = rng.random((size, dimension))
        step_draws = rng.random((size, dimension))
        winners = [int(a) if values[a] <= values[b] else int(b) for a, b in tournaments]
        children = np.empty((size, dimension))
        for k in range(size // 2):
            a, b = population[winners[2 * k]], population[winners[2 * k + 1]]
            for j in range(dimension):
                image = a[j] + ((b[j] - a[j] + 0.5) % 1.0 - 0.5) if periodic[j] else b[j]
                if recombine[k] < p_cross and take_part[k, j] < 0.5:
                    u = spread_draws[k, j]
                    beta = (2 * u) ** (1 / (eta_cross + 1)) if u <= 0.5 else (1 / (2 * (1 - u))) ** (1 / (eta_cross + 1))
                    first = 0.5 * ((1 + beta) * a[j] + (1 - beta) * image)
                    second = 0.5 * ((1 - beta) * a[j] + (1 + beta) * image)
                else:
                    first, second = a[j], b[j]
                children[2 * k, j], children[2 * k + 1, j] = first, second
        for i in range(size):
            for j in range(dimension):
                if mutate[i, j] < p_mutation:
                    u = step_draws[i, j]
                    children[i, j] += ((2 * u) ** (1 / (eta_mutation + 1)) - 1 if u < 0.5
                                       else 1 - (2 * (1 - u)) ** (1 / (eta_mutation + 1)))
                if periodic[j]:
                    children[i, j] = children[i, j] % 1.0
                    children[i, j] = 0.0 if children[i, j] >= 1.0 else children[i, j]
                else:
                    children[i, j] = min(1.0, max(0.0, children[i, j]))
        batches.append(children.copy())
        child_values = [function(row) for row in children]
        elite = min(range(size), key=lambda index: (values[index], index))
        elite_value, elite_point = values[elite], population[elite].copy()
        population, values = children, child_values
        if elite_value < min(values):
            worst = max(range(size), key=lambda index: (values[index], index))
            population[worst], values[worst] = elite_point, elite_value
    return batches


class GeneticAlgorithmTests(unittest.TestCase):
    def test_evaluated_points_match_independent_replication(self) -> None:
        size, dimension, periodic = 8, 3, [False, True, False]
        settings = GeneticAlgorithmSettings(population_size=size, crossover_probability=0.8, crossover_index=2.0,
                                            mutation_probability=0.4, mutation_index=3.0)

        def function(u: np.ndarray) -> float:
            return float(u[0] + (1.0 - math.cos(2 * math.pi * (u[1] - 0.97))) + (1.0 - u[2]))

        problem = RecordingProblem(function, unit_bounds(dimension, periodic))
        result = GeneticAlgorithm(settings).minimize(problem, TerminationCriteria(size * 5), seed=17)
        batches = replicate_ga(17, size, dimension, 0.8, 2.0, 0.4, 3.0, periodic, function, 4)
        # Array and scalar pow may differ in the last bit; an operator or draw-order error is O(1).
        assert_allclose(np.array(problem.points), np.concatenate(batches), rtol=0.0, atol=1e-12)
        self.assertEqual(result.generations, 4)

    def test_sbx_and_polynomial_operators(self) -> None:
        uniforms = np.array([0.0, 0.1, 0.5, 0.9, 0.999999])
        for index in (0.0, 2.0, 15.0):
            beta = _sbx_spread(uniforms, index)
            self.assertEqual(beta[0], 0.0)
            self.assertEqual(beta[2], 1.0)                      # u = 1/2: children equal the parents
            self.assertTrue(np.all(np.diff(beta) > 0))
            # beta(u) and beta(1-u) are reciprocal: contraction and expansion are symmetric.
            self.assertAlmostEqual(beta[1] * beta[3], 1.0, places=12)
            delta = _polynomial_step(uniforms, index)
            self.assertEqual(delta[0], -1.0)
            self.assertEqual(delta[2], 0.0)
            self.assertTrue(np.all((delta >= -1.0) & (delta < 1.0)))
            self.assertAlmostEqual(delta[1], -delta[3], places=12)
        a, b, beta = 0.3, 0.8, _sbx_spread(np.array([0.2]), 15.0)[0]
        first, second = 0.5 * ((1 + beta) * a + (1 - beta) * b), 0.5 * ((1 - beta) * a + (1 + beta) * b)
        self.assertAlmostEqual(first + second, a + b, places=15)        # mean preserving
        self.assertAlmostEqual(abs(second - first), beta * abs(b - a), places=15)

    def test_elitism_keeps_the_best_member_in_the_population(self) -> None:
        bounds = Bounds([-2.0, -1.0], [3.0, 4.0])
        result = GeneticAlgorithm(GeneticAlgorithmSettings(record_population=True, mutation_probability=1.0)).minimize(
            FunctionProblem(sphere, bounds), TerminationCriteria(600), seed=6)
        for generation, population in enumerate(result.history.populations):
            best_in_population = min(sphere(bounds.from_unit(row)) for row in population)
            self.assertEqual(best_in_population, result.history.generation_best[generation])

    def test_settings_validation(self) -> None:
        for kwargs in ({"population_size": 3}, {"population_size": 5}, {"population_size": True},
                       {"crossover_probability": 1.1}, {"crossover_probability": -0.1}, {"crossover_index": -1.0},
                       {"mutation_probability": 1.5}, {"mutation_index": math.inf}, {"mutation_index": True},
                       {"record_population": 0}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                GeneticAlgorithmSettings(**kwargs)
        settings = GeneticAlgorithmSettings()
        self.assertEqual(settings.resolved_population_size(3), 30)
        self.assertAlmostEqual(settings.resolved_mutation_probability(4), 0.25)


class NelderMeadTests(unittest.TestCase):
    def test_converges_to_the_minimum_of_a_convex_quadratic(self) -> None:
        # f = (x-1)^2 + 10 (y+0.5)^2 on [-2,3]x[-1,4]: unique minimum at (1, -0.5), value 0.
        bounds = Bounds([-2.0, -1.0], [3.0, 4.0])
        problem = FunctionProblem(lambda x: float((x[0] - 1.0) ** 2 + 10.0 * (x[1] + 0.5) ** 2), bounds)
        start = tuple(bounds.to_unit([2.5, 3.0]))
        result = NelderMead(NelderMeadSettings(initial_point=start)).minimize(
            problem, TerminationCriteria(5_000, x_tol=1e-9), seed=0)
        self.assertIs(result.termination_reason, TerminationReason.POPULATION_CONVERGED)
        # Simplex diameter 1e-9 of a range of 5: the vertices are within 5e-9 of each other.
        assert_allclose(result.x_best, [1.0, -0.5], atol=1e-7)
        self.assertLess(result.f_best, 1e-13)
        self.assertEqual(result.settings["start"], "given initial point")
        operations = result.settings["operations"]
        self.assertEqual(sum(operations.values()), result.generations)
        self.assertGreater(operations["reflection"] + operations["expansion"], 0)
        self.assertGreater(operations["inside_contraction"] + operations["outside_contraction"], 0)

    def test_minimum_on_the_boundary_is_reached_by_mirroring(self) -> None:
        # Corner minimum of a linear function and face minimum of a convex quadratic whose
        # unconstrained minimum (0.3, -0.2) lies outside the box: constrained minimum (0.3, 0), f = 0.04.
        cases = (
            (lambda u: float(u[0] + 2.0 * u[1]), (0.4, 0.3), [0.0, 0.0], 0.0),
            (lambda u: float((u[0] - 0.3) ** 2 + (u[1] + 0.2) ** 2), (0.8, 0.7), [0.3, 0.0], 0.04),
        )
        for function, start, minimizer, minimum in cases:
            with self.subTest(minimizer=minimizer):
                problem = RecordingProblem(function, unit_bounds(2))
                result = NelderMead(NelderMeadSettings(initial_point=start)).minimize(
                    problem, TerminationCriteria(5_000, x_tol=1e-10), seed=0)
                self.assertIs(result.termination_reason, TerminationReason.POPULATION_CONVERGED)
                # The simplex ends with diameter 1e-10; the face minimum is quadratic along the
                # face, so its abscissa is resolved to about sqrt(1e-10 * diameter) ~ 1e-8.
                assert_allclose(result.x_best, minimizer, atol=1e-7)
                self.assertAlmostEqual(result.f_best, minimum, places=9)
                points = np.array(problem.points)
                self.assertTrue(np.all((points >= 0.0) & (points <= 1.0)))

    def test_mirror(self) -> None:
        from orbital_optimizer.optimization.nelder_mead import _mirror
        assert_allclose(_mirror(np.array([-0.2, 0.0, 0.4, 1.0, 1.3, -1.5, 2.7])),
                        [0.2, 0.0, 0.4, 1.0, 0.7, 0.5, 0.0], atol=1e-15)

    def test_periodic_variable_crosses_the_seam(self) -> None:
        # Minimum at u = 0.98; from u0 = 0.05 the short way goes through the seam at 0 = 1.
        def function(u: np.ndarray) -> float:
            return float(1.0 - math.cos(2 * math.pi * (u[0] - 0.98)) + (u[1] - 0.5) ** 2)

        problem = RecordingProblem(function, unit_bounds(2, [True, False]))
        result = NelderMead(NelderMeadSettings(initial_point=(0.05, 0.5))).minimize(
            problem, TerminationCriteria(2_000, x_tol=1e-10), seed=0)
        assert_allclose(result.x_best, [0.98, 0.5], atol=1e-6)
        points = np.array(problem.points)
        self.assertTrue(np.all((points[:, 0] >= 0.0) & (points[:, 0] < 1.0)))
        # The search never visits the long way round (0.1 < u < 0.9).
        self.assertFalse(np.any((points[:, 0] > 0.1) & (points[:, 0] < 0.9)))

    def test_initial_simplex_and_random_start(self) -> None:
        problem = RecordingProblem(sphere, unit_bounds(3, [False, True, False]))
        NelderMead(NelderMeadSettings(initial_point=(0.5, 0.99, 0.98), initial_step=0.05)).minimize(
            problem, TerminationCriteria(4), seed=0)
        expected = [[0.5, 0.99, 0.98], [0.55, 0.99, 0.98], [0.5, 0.04, 0.98], [0.5, 0.99, 0.93]]
        assert_allclose(np.array(problem.points), expected, atol=1e-15)
        random_start = RecordingProblem(sphere, unit_bounds(3))
        result = NelderMead().minimize(random_start, TerminationCriteria(4), seed=21)
        assert_array_equal(random_start.points[0], np.random.Generator(np.random.PCG64(21)).random(3))
        self.assertEqual(result.settings["start"], "seeded uniform random point")

    def test_settings_validation(self) -> None:
        for kwargs in ({"initial_point": (0.5, 1.5)}, {"initial_point": (math.nan,)}, {"initial_point": "ab"},
                       {"initial_point": ()}, {"initial_step": 0.0}, {"initial_step": 0.6}, {"initial_step": True},
                       {"record_population": 1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                NelderMeadSettings(**kwargs)
        with self.assertRaises(ValueError):       # start of the wrong dimension
            NelderMead(NelderMeadSettings(initial_point=(0.5,))).minimize(
                FunctionProblem(sphere, unit_bounds(2)), TerminationCriteria(100), seed=0)


class HybridTests(unittest.TestCase):
    def setUp(self) -> None:
        self.benchmark = rosenbrock_benchmark(2, -2.048, 4.096)

    def test_budget_split_and_prefix_of_the_global_phase(self) -> None:
        budget, share = 1_000, 800
        hybrid = HybridOptimizer().minimize(self.benchmark.problem(), TerminationCriteria(budget), seed=7)
        alone = DifferentialEvolution().minimize(self.benchmark.problem(), TerminationCriteria(share), seed=7)
        self.assertEqual(hybrid.nfev, budget)
        self.assertEqual(hybrid.settings["global_evaluations"], share)
        self.assertEqual(hybrid.settings["local_evaluations"], budget - share)
        self.assertEqual(hybrid.settings["local_phase"], "executed")
        assert_array_equal(hybrid.history.best_so_far[:share], alone.history.best_so_far)
        assert_array_equal(hybrid.history.status[:share], alone.history.status)
        self.assertLessEqual(hybrid.f_best, alone.f_best)
        self.assertTrue(np.all(np.diff(hybrid.history.best_so_far) <= 0))
        self.assertTrue(np.all(np.diff(hybrid.history.generation_best) <= 0))
        self.assertTrue(np.all(np.diff(hybrid.history.generation_nfev) > 0))
        self.assertEqual(hybrid.history.generation_nfev[-1], budget)
        self.assertEqual(hybrid.f_best, hybrid.history.best_so_far[-1])
        self.assertEqual(self.benchmark.function(hybrid.x_best), hybrid.f_best)
        self.assertEqual(sum(hybrid.status_counts.values()), budget)
        local = hybrid.settings["local_settings"]
        self.assertEqual(local["start"], "given initial point")
        assert_allclose(local["initial_point"], alone.u_best)
        spread = float(alone.history.generation_spread[-1])
        self.assertEqual(local["initial_step"], min(0.05, max(1e-6, spread)))

    def test_local_phase_is_skipped_when_it_cannot_help(self) -> None:
        problem = FunctionProblem(sphere, Bounds([-2.0, -1.0], [3.0, 4.0]))
        reached = HybridOptimizer().minimize(problem, TerminationCriteria(50_000, target_objective=1e-2), seed=0)
        self.assertIs(reached.termination_reason, TerminationReason.TARGET_REACHED)
        self.assertIn("target reached", reached.settings["local_phase"])
        self.assertLess(reached.nfev, 50_000)
        tight = HybridOptimizer(settings=HybridSettings(global_fraction=0.99)).minimize(
            problem, TerminationCriteria(200), seed=0)
        self.assertEqual(tight.nfev, 198)                       # floor(0.99*200); 2 < D + 1 remain
        self.assertIn("fewer than D + 1", tight.settings["local_phase"])
        small = HybridOptimizer().minimize(problem, TerminationCriteria(22), seed=0)
        self.assertEqual(small.nfev, 22)                        # share 17 < NP: global phase takes all
        self.assertIn("skipped", small.settings["local_phase"])

    def test_other_global_optimizer_and_validation(self) -> None:
        hybrid = HybridOptimizer(ParticleSwarm(), HybridSettings(global_fraction=0.5))
        self.assertEqual(hybrid.name, "Hybrid(PSO/gbest/inertia@1 -> NelderMead/Lagarias1998@1)@1")
        result = hybrid.minimize(self.benchmark.problem(), TerminationCriteria(400), seed=2)
        self.assertEqual((result.nfev, result.settings["global_evaluations"]), (400, 200))
        self.assertEqual(result.algorithm, hybrid.name)
        for fraction in (0.0, 1.0, -0.5, True, "0.5"):
            with self.subTest(fraction=fraction), self.assertRaises(ValueError):
                HybridSettings(global_fraction=fraction)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            HybridOptimizer(global_optimizer=object())  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            HybridOptimizer(settings="settings")  # type: ignore[arg-type]


class ConstantStatusProblem:
    """Half the box is INVALID_DOMAIN; used to check key-based comparisons with ties."""

    name = "half"
    bounds = unit_bounds(2)

    def evaluate(self, x):
        if x[0] < 0.5:
            return Evaluation(x, math.nan, EvaluationStatus.INVALID_DOMAIN, "outside domain")
        return Evaluation(x, float(x[1]), EvaluationStatus.VALID)


class DomainStatusTests(unittest.TestCase):
    def test_invalid_domain_is_counted_and_avoided(self) -> None:
        for optimizer in global_optimizers():
            with self.subTest(optimizer=optimizer.name):
                result = optimizer.minimize(ConstantStatusProblem(), TerminationCriteria(300), seed=2)
                self.assertGreater(result.status_counts["invalid_domain"], 0)
                self.assertGreaterEqual(result.x_best[0], 0.5)

    def test_local_search_treats_invalid_points_as_worst(self) -> None:
        """Started next to the invalid half, the simplex probes it and never keeps such a point.

        Only the contract is asserted (best point valid, objective improved on
        the start value 0.6). Convergence to the minimum x1 = 0 is *not*
        asserted: on 2026-11-10 this run stalled at f = 0.033 because the
        simplex collapsed against the invalid half-plane, the known weakness
        of Nelder-Mead next to a discontinuous barrier. That limitation is
        documented in optimization/nelder_mead.py.
        """
        result = NelderMead(NelderMeadSettings(initial_point=(0.52, 0.6))).minimize(
            ConstantStatusProblem(), TerminationCriteria(400), seed=0)
        self.assertGreater(result.status_counts["invalid_domain"], 0)
        self.assertGreaterEqual(result.x_best[0], 0.5)
        self.assertLess(result.f_best, 0.6)
        self.assertEqual(result.f_best, result.x_best[1])

    def test_local_search_started_in_an_invalid_region_reports_no_valid_point(self) -> None:
        # A local method has no mechanism to leave a region where every point ties as invalid.
        result = NelderMead(NelderMeadSettings(initial_point=(0.1, 0.5))).minimize(
            ConstantStatusProblem(), TerminationCriteria(200), seed=0)
        self.assertFalse(result.found_valid)
        self.assertEqual(result.status_counts["invalid_domain"], 200)


if __name__ == "__main__":
    unittest.main()
