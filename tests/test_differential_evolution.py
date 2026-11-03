"""DE/rand/1/bin: operators, documented draw order, invalid handling, budget, termination, reproducibility."""
import hashlib
import json
import math
import subprocess
import sys
import unittest

import numpy as np
from numpy.testing import assert_array_equal

from orbital_optimizer.optimization import (
    ALGORITHM, TERMINATION_PRECEDENCE, Bounds, DifferentialEvolution, DifferentialEvolutionSettings,
    Evaluation, EvaluationStatus, FunctionProblem, OptimizationResult, TerminationCriteria,
    TerminationReason, rastrigin_benchmark, rosenbrock_benchmark, sphere,
)
from orbital_optimizer.optimization.differential_evolution import (
    _binomial_crossover, _repair, _select_donors, _trial_wins,
)


class RecordingProblem:
    """Delegate to FunctionProblem and keep a copy of every evaluated x."""

    def __init__(self, function, bounds: Bounds, **kwargs) -> None:
        self.inner = FunctionProblem(function, bounds, **kwargs)
        self.bounds = bounds
        self.name = "recording"
        self.points: list[np.ndarray] = []

    def evaluate(self, x: np.ndarray) -> Evaluation:
        self.points.append(np.array(x, copy=True))
        return self.inner.evaluate(x)


def unit_bounds(dimension: int, periodic=None) -> Bounds:
    """[0,1]^D, for which the physical and normalized vectors coincide exactly."""
    return Bounds(np.zeros(dimension), np.ones(dimension),
                  periodic=None if periodic is None else np.array(periodic, dtype=bool))


def run(problem, budget: int, seed: int = 0, **settings) -> OptimizationResult:
    return DifferentialEvolution(DifferentialEvolutionSettings(**settings)).minimize(
        problem, TerminationCriteria(max_evaluations=budget), seed=seed)


def constant(x: np.ndarray) -> float:
    return 1.0


def replicate(seed: int, size: int, dimension: int, mutation: float, crossover: float,
              periodic: list[bool], function, generations: int) -> tuple[list[np.ndarray], int]:
    """Independent scalar re-implementation of the documented draws and operators."""
    rng = np.random.Generator(np.random.PCG64(seed))
    population = rng.random((size, dimension))
    values = [function(p) for p in population]
    batches, repairs = [population.copy()], 0
    for _ in range(generations):
        donors = []
        for i in range(size):
            choice = rng.choice(size - 1, size=3, replace=False)
            donors.append([c + 1 if c >= i else c for c in choice])
        jrand = rng.integers(0, dimension, size=size)
        uniforms = rng.random((size, dimension))
        trials = population.copy()
        for i in range(size):
            r1, r2, r3 = donors[i]
            for j in range(dimension):
                if uniforms[i, j] <= crossover or j == jrand[i]:
                    value = population[r1, j] + mutation * (population[r2, j] - population[r3, j])
                    if periodic[j]:
                        repairs += not 0.0 <= value < 1.0
                        value = value % 1.0
                        value = 0.0 if value >= 1.0 else value
                    elif value < 0.0:
                        repairs += 1
                        value = population[i, j] / 2.0
                    elif value > 1.0:
                        repairs += 1
                        value = (population[i, j] + 1.0) / 2.0
                    trials[i, j] = value
        batches.append(trials.copy())
        for i in range(size):
            trial_value = function(trials[i])
            if trial_value <= values[i]:
                population[i], values[i] = trials[i], trial_value
    return batches, repairs


class OperatorTests(unittest.TestCase):
    def test_donors_are_distinct_and_differ_from_target(self) -> None:
        rng = np.random.Generator(np.random.PCG64(3))
        for size in (4, 5, 20):
            for _ in range(50):
                donors = _select_donors(rng, size)
                self.assertEqual(donors.shape, (size, 3))
                for i, row in enumerate(donors):
                    self.assertEqual(len(set(row.tolist()) | {i}), 4)
                    self.assertTrue(np.all((0 <= row) & (row < size)))

    def test_binomial_crossover_extremes(self) -> None:
        parents, mutants = np.zeros((6, 4)), np.ones((6, 4))
        uniforms = np.full((6, 4), 0.5)
        jrand = np.array([0, 1, 2, 3, 0, 1])
        assert_array_equal(_binomial_crossover(parents, mutants, uniforms, jrand, 1.0), mutants)
        trials = _binomial_crossover(parents, mutants, uniforms, jrand, 0.0)
        assert_array_equal(trials.sum(axis=1), np.ones(6))
        assert_array_equal(np.argmax(trials, axis=1), jrand)
        # R_ij <= CR selects the mutant component (inclusive comparison).
        assert_array_equal(_binomial_crossover(parents, mutants, uniforms, jrand, 0.5), mutants)

    def test_boundary_repair(self) -> None:
        parents = np.array([[0.2, 0.6, 0.3, 0.3]])
        trials = np.array([[-0.5, 1.7, -0.25, 1.25]])
        repaired = _repair(trials, parents, np.array([False, False, True, True]))
        assert_array_equal(repaired, [[0.1, 0.8, 0.75, 0.25]])
        tiny = _repair(np.array([[-1e-20]]), np.array([[0.5]]), np.array([True]))
        self.assertEqual(tiny[0, 0], 0.0)
        inside = np.array([[0.0, 1.0, 0.4]])
        assert_array_equal(_repair(inside, inside, np.array([False, False, True])), inside)

    def test_selection_key(self) -> None:
        self.assertTrue(_trial_wins(True, 1.0, True, 1.0))     # ties favour the trial
        self.assertTrue(_trial_wins(True, 0.5, True, 1.0))
        self.assertFalse(_trial_wins(True, 2.0, True, 1.0))
        self.assertTrue(_trial_wins(True, 1e300, False, math.inf))
        self.assertFalse(_trial_wins(False, math.inf, True, 1e300))
        self.assertTrue(_trial_wins(False, math.inf, False, math.inf))  # neutral drift


class DrawOrderTests(unittest.TestCase):
    def test_evaluated_points_match_independent_replication(self) -> None:
        # Two non-periodic and one periodic variable; F=0.9 forces boundary repairs.
        size, dimension, mutation, crossover = 6, 3, 0.9, 0.7
        periodic = [False, True, False]
        problem = RecordingProblem(sphere, unit_bounds(dimension, periodic))
        result = DifferentialEvolution(DifferentialEvolutionSettings(size, mutation, crossover)).minimize(
            problem, TerminationCriteria(max_evaluations=size * 4), seed=11)
        batches, repairs = replicate(11, size, dimension, mutation, crossover, periodic, sphere, 3)
        self.assertGreater(repairs, 0)
        assert_array_equal(np.array(problem.points), np.concatenate(batches))
        self.assertEqual(result.generations, 3)


class InvalidEvaluationTests(unittest.TestCase):
    @staticmethod
    def partly_invalid(x: np.ndarray) -> float:
        if x[0] < 0.2:
            return -math.inf
        if x[0] > 0.9:
            return math.nan
        if x[1] > 0.9:
            return math.inf
        return float(np.sum((x - 0.5) ** 2))

    def test_non_finite_values_never_become_best_and_are_counted(self) -> None:
        problem = RecordingProblem(self.partly_invalid, unit_bounds(2))
        result = run(problem, 600, seed=5)
        points = np.array(problem.points)
        expected_invalid = int(np.count_nonzero((points[:, 0] < 0.2) | (points[:, 0] > 0.9) | (points[:, 1] > 0.9)))
        self.assertGreater(expected_invalid, 0)
        self.assertEqual(result.status_counts["non_finite"], expected_invalid)
        self.assertEqual(result.n_invalid, expected_invalid)
        self.assertEqual(sum(result.status_counts.values()), result.nfev)
        self.assertTrue(result.found_valid)
        self.assertTrue(math.isfinite(result.f_best) and result.f_best >= 0.0)
        self.assertTrue(0.2 <= result.x_best[0] <= 0.9 and result.x_best[1] <= 0.9)
        valid_values = [self.partly_invalid(p) for p in points]
        finite = [v for v in valid_values if math.isfinite(v)]
        self.assertEqual(result.f_best, min(finite))

    def test_all_invalid_problem(self) -> None:
        result = run(FunctionProblem(lambda x: math.nan, unit_bounds(2)), 75, seed=1)
        self.assertFalse(result.found_valid)
        self.assertIsNone(result.x_best)
        self.assertIsNone(result.u_best)
        self.assertIsNone(result.best_evaluation)
        self.assertEqual(result.f_best, math.inf)
        self.assertIs(result.termination_reason, TerminationReason.MAX_EVALUATIONS)
        self.assertEqual(result.nfev, 75)
        self.assertEqual(result.status_counts["non_finite"], 75)
        self.assertTrue(np.all(np.isinf(result.history.best_so_far)))
        json.dumps(result.to_dict(), allow_nan=False)

    def test_domain_and_numerical_statuses_are_counted_separately(self) -> None:
        class DomainProblem:
            name = "domain"
            bounds = unit_bounds(2)

            def evaluate(self, x):
                if x[0] < 0.5:
                    return Evaluation(x, math.nan, EvaluationStatus.INVALID_DOMAIN, "outside domain")
                return Evaluation(x, float(x[1]), EvaluationStatus.VALID)

        result = run(DomainProblem(), 200, seed=2)
        self.assertGreater(result.status_counts["invalid_domain"], 0)
        self.assertGreaterEqual(result.x_best[0], 0.5)

        def overflow(x):
            raise OverflowError("too large")

        failing = FunctionProblem(overflow, unit_bounds(2), invalid_exceptions=(OverflowError,))
        self.assertEqual(run(failing, 40).status_counts["numerical_failure"], 40)

    def test_unexpected_exceptions_and_bad_returns_abort(self) -> None:
        def broken(x):
            raise RuntimeError("bug")

        with self.assertRaises(RuntimeError):
            run(FunctionProblem(broken, unit_bounds(2)), 40)

        class BadProblem:
            name = "bad"
            bounds = unit_bounds(2)

            def evaluate(self, x):
                return 1.0

        with self.assertRaises(TypeError):
            run(BadProblem(), 40)


class ProblemContractTests(unittest.TestCase):
    """The budget checks the problem's Evaluation and isolates its own copies of x."""

    def test_evaluation_for_another_point_is_rejected(self) -> None:
        class ShiftedProblem:
            name = "shifted"
            bounds = unit_bounds(2)

            def evaluate(self, x):
                return Evaluation(np.asarray(x) * 0.5, 1.0, EvaluationStatus.VALID)

        class TruncatedProblem:
            name = "truncated"
            bounds = unit_bounds(2)

            def evaluate(self, x):
                return Evaluation(np.asarray(x)[:1], 1.0, EvaluationStatus.VALID)

        for problem, message in ((ShiftedProblem(), "differs"), (TruncatedProblem(), "shape")):
            with self.subTest(problem=problem.name), self.assertRaisesRegex(ValueError, message):
                run(problem, 40)

    def test_problem_mutating_its_argument_cannot_alter_the_result(self) -> None:
        class MutatingProblem:
            """Re-enables writing on every argument it received and overwrites them later."""

            name = "mutating"
            bounds = Bounds([-2.0, -1.0], [3.0, 4.0])

            def __init__(self) -> None:
                self.kept: list[np.ndarray] = []

            def evaluate(self, x):
                for old in self.kept:
                    old.setflags(write=True)
                    old[:] = 99.0
                evaluation = Evaluation(x, sphere(x), EvaluationStatus.VALID)
                self.kept.append(x)
                return evaluation

        problem = MutatingProblem()
        result = run(problem, 200, seed=3)
        for old in problem.kept:
            old.setflags(write=True)
            old[:] = 99.0
        self.assertTrue(all(np.all(old == 99.0) for old in problem.kept))   # the mutation really happened
        self.assertTrue(np.all(result.x_best != 99.0))
        assert_array_equal(result.x_best, problem.bounds.from_unit(result.u_best))
        assert_array_equal(result.best_evaluation.x, result.x_best)
        self.assertEqual(sphere(result.x_best), result.f_best)
        self.assertFalse(any(np.shares_memory(old, result.x_best) for old in problem.kept))


class ScriptedProblem:
    """Status by call index: listed calls are INVALID_DOMAIN, every other call is VALID with value 1.0."""

    name = "scripted"

    def __init__(self, dimension: int, invalid_calls: set[int]) -> None:
        self.bounds = unit_bounds(dimension)
        self.invalid_calls = invalid_calls
        self.calls = 0

    def evaluate(self, x):
        call, self.calls = self.calls, self.calls + 1
        if call in self.invalid_calls:
            return Evaluation(x, math.nan, EvaluationStatus.INVALID_DOMAIN, "scripted")
        return Evaluation(x, 1.0, EvaluationStatus.VALID)


class BudgetAndTerminationTests(unittest.TestCase):
    def setUp(self) -> None:
        # Asymmetric box: the optimum at the origin is not the centre of the box.
        self.problem = FunctionProblem(sphere, Bounds([-2.0, -1.0], [3.0, 4.0]))

    def test_budget_is_respected_exactly(self) -> None:
        for budget in (20, 21, 57, 200):
            with self.subTest(budget=budget):
                result = run(self.problem, budget, seed=4)
                self.assertEqual(result.nfev, budget)
                self.assertIs(result.termination_reason, TerminationReason.MAX_EVALUATIONS)
                self.assertEqual(len(result.history.best_so_far), budget)
                self.assertEqual(len(result.history.status), budget)
                self.assertEqual(result.history.generation_nfev[-1], budget)
        partial = run(self.problem, 57)
        self.assertEqual(partial.generations, 2)            # 20 initial + 20 + 17
        assert_array_equal(partial.history.generation_nfev, [20, 40, 57])
        with self.assertRaises(ValueError):
            run(self.problem, 19)

    def test_max_generations(self) -> None:
        result = DifferentialEvolution().minimize(
            self.problem, TerminationCriteria(max_evaluations=10_000, max_generations=3), seed=0)
        self.assertIs(result.termination_reason, TerminationReason.MAX_GENERATIONS)
        self.assertEqual((result.generations, result.nfev), (3, 80))

    def test_population_convergence(self) -> None:
        result = DifferentialEvolution().minimize(
            self.problem, TerminationCriteria(max_evaluations=50_000, x_tol=1e-3), seed=0)
        self.assertIs(result.termination_reason, TerminationReason.POPULATION_CONVERGED)
        spreads = result.history.generation_spread
        self.assertLessEqual(spreads[-1], 1e-3)
        self.assertTrue(np.all(spreads[:-1] > 1e-3))
        self.assertLess(result.nfev, 50_000)

    def test_objective_convergence_requires_an_all_valid_population(self) -> None:
        def half_invalid(x):
            return math.nan if x[0] > 0.5 else 1.0

        result = DifferentialEvolution().minimize(
            FunctionProblem(half_invalid, unit_bounds(2)),
            TerminationCriteria(max_evaluations=5_000, f_tol=0.0), seed=0)
        self.assertLess(result.history.generation_valid[0], 20)
        self.assertIs(result.termination_reason, TerminationReason.OBJECTIVE_CONVERGED)
        self.assertGreater(result.generations, 0)
        self.assertEqual(result.history.generation_valid[-1], 20)
        self.assertTrue(np.all(result.history.generation_valid[:-1] < 20))

    def test_stagnation(self) -> None:
        result = DifferentialEvolution().minimize(
            FunctionProblem(constant, unit_bounds(2)),
            TerminationCriteria(max_evaluations=5_000, stagnation_generations=3), seed=0)
        self.assertIs(result.termination_reason, TerminationReason.STAGNATION)
        self.assertEqual(result.generations, 3)

    def test_target_reached(self) -> None:
        result = DifferentialEvolution().minimize(
            self.problem, TerminationCriteria(max_evaluations=50_000, target_objective=1e-6), seed=0)
        self.assertIs(result.termination_reason, TerminationReason.TARGET_REACHED)
        self.assertLessEqual(result.f_best, 1e-6)
        self.assertGreater(result.history.generation_best[-2], 1e-6)

    def test_precedence_when_several_criteria_hold(self) -> None:
        self.assertEqual(TERMINATION_PRECEDENCE, (
            TerminationReason.MAX_EVALUATIONS, TerminationReason.TARGET_REACHED,
            TerminationReason.POPULATION_CONVERGED, TerminationReason.OBJECTIVE_CONVERGED,
            TerminationReason.STAGNATION, TerminationReason.MAX_GENERATIONS,
        ))
        problem = FunctionProblem(constant, unit_bounds(2))
        cases = [
            (TerminationCriteria(20, target_objective=2.0), TerminationReason.MAX_EVALUATIONS),
            (TerminationCriteria(1000, target_objective=2.0, x_tol=1.0, f_tol=0.0), TerminationReason.TARGET_REACHED),
            (TerminationCriteria(1000, x_tol=1.0, f_tol=0.0), TerminationReason.POPULATION_CONVERGED),
            (TerminationCriteria(1000, stagnation_generations=2, max_generations=2), TerminationReason.STAGNATION),
        ]
        for criteria, expected in cases:
            with self.subTest(expected=expected):
                result = DifferentialEvolution().minimize(problem, criteria, seed=0)
                self.assertIs(result.termination_reason, expected)

    def test_objective_convergence_precedes_stagnation_and_max_generations_after_generation_zero(self) -> None:
        # Call 0 is invalid, so f_tol cannot hold at generation 0 (19/20 valid). In generation 1
        # every trial is valid with f = 1.0: the population becomes all valid with zero spread,
        # and the best (1.0) did not decrease, so stagnation (1 generation) and
        # max_generations=1 hold at the same check as f_tol.
        def scripted_run(**criteria) -> OptimizationResult:
            return DifferentialEvolution().minimize(
                ScriptedProblem(2, {0}), TerminationCriteria(1000, **criteria), seed=0)

        for competing, alone in (({"stagnation_generations": 1}, TerminationReason.STAGNATION),
                                 ({"max_generations": 1}, TerminationReason.MAX_GENERATIONS)):
            with self.subTest(competing=competing):
                control = scripted_run(**competing)
                self.assertIs(control.termination_reason, alone)
                self.assertEqual(control.generations, 1)
                both = scripted_run(f_tol=0.0, **competing)
                self.assertIs(both.termination_reason, TerminationReason.OBJECTIVE_CONVERGED)
                self.assertEqual(both.generations, 1)
                assert_array_equal(both.history.generation_valid, [19, 20])
                assert_array_equal(both.history.generation_best, [1.0, 1.0])

    def test_invalid_configuration_is_rejected(self) -> None:
        for kwargs in ({"mutation": 0.0}, {"mutation": 2.1}, {"mutation": math.nan}, {"mutation": True},
                       {"crossover": -0.1}, {"crossover": 1.1}, {"population_size": 3},
                       {"population_size": True}, {"population_size": 4.0}, {"record_population": 1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                DifferentialEvolutionSettings(**kwargs)
        for kwargs in ({"max_evaluations": 0}, {"max_evaluations": True}, {"max_evaluations": 1.5},
                       {"max_evaluations": 10, "max_generations": 0}, {"max_evaluations": 10, "target_objective": math.nan},
                       {"max_evaluations": 10, "x_tol": -1.0}, {"max_evaluations": 10, "f_tol": math.inf},
                       {"max_evaluations": 10, "stagnation_generations": 0}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                TerminationCriteria(**kwargs)
        for seed in (-1, True, 1.5, None, "1"):
            with self.subTest(seed=seed), self.assertRaises(ValueError):
                DifferentialEvolution().minimize(self.problem, TerminationCriteria(100), seed=seed)  # type: ignore[arg-type]
        with self.assertRaises(TypeError):
            DifferentialEvolution().minimize(object(), TerminationCriteria(100), seed=0)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            DifferentialEvolution().minimize(self.problem, 100, seed=0)  # type: ignore[arg-type]
        self.assertEqual(DifferentialEvolutionSettings().resolved_population_size(3), 30)


class ReproducibilityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.benchmark = rastrigin_benchmark(3, lower=-2.56, upper=5.12)

    @staticmethod
    def comparable(result: OptimizationResult) -> dict:
        data = result.to_dict()
        data.pop("wall_time_s")
        return data

    def test_same_seed_is_bitwise_identical(self) -> None:
        first = run(self.benchmark.problem(), 900, seed=123, record_population=True)
        second = run(self.benchmark.problem(), 900, seed=123, record_population=True)
        self.assertEqual(self.comparable(first), self.comparable(second))
        self.assertEqual(first.x_best.tobytes(), second.x_best.tobytes())
        for name in ("best_so_far", "status", "generation_best", "populations"):
            self.assertEqual(getattr(first.history, name).tobytes(), getattr(second.history, name).tobytes())

    def test_different_seeds_differ(self) -> None:
        first = RecordingProblem(self.benchmark.function, self.benchmark.bounds)
        second = RecordingProblem(self.benchmark.function, self.benchmark.bounds)
        run(first, 300, seed=1)
        run(second, 300, seed=2)
        self.assertFalse(np.array_equal(np.array(first.points), np.array(second.points)))

    def test_prefix_property_across_budgets(self) -> None:
        short = RecordingProblem(self.benchmark.function, self.benchmark.bounds)
        long = RecordingProblem(self.benchmark.function, self.benchmark.bounds)
        a, b = run(short, 200, seed=9), run(long, 537, seed=9)
        assert_array_equal(np.array(short.points), np.array(long.points)[:200])
        assert_array_equal(a.history.best_so_far, b.history.best_so_far[:200])
        assert_array_equal(a.history.status, b.history.status[:200])


CROSS_PROCESS_SCRIPT = r'''
import hashlib, json
import numpy
from orbital_optimizer.optimization import DifferentialEvolution, TerminationCriteria, rosenbrock_benchmark
benchmark = rosenbrock_benchmark(2, -2.048, 4.096)
result = DifferentialEvolution().minimize(benchmark.problem(), TerminationCriteria(max_evaluations=600), seed=2026)
record = result.to_dict()
record.pop("wall_time_s")
full = hashlib.sha256(json.dumps(record, sort_keys=True, allow_nan=False).encode()).hexdigest()
record.pop("environment")
portable = hashlib.sha256(json.dumps(record, sort_keys=True, allow_nan=False).encode()).hexdigest()
print(json.dumps({"full": full, "portable": portable, "f_best": result.f_best.hex(), "numpy": numpy.__version__}))
'''

# Drift detector, not a correctness oracle: values produced by this implementation on
# 2026-11-03 (Windows 11, Python 3.12.1, NumPy 2.3.3) with CROSS_PROCESS_SCRIPT. The
# "portable" hash excludes wall time and the environment record. Rosenbrock-2D uses only
# elementwise IEEE arithmetic (no libm transcendental call, no multi-element reduction),
# so the pin is expected to hold on other platforms with the same NumPy minor version.
PINNED_NUMPY = (2, 3)
PINNED_F_BEST_HEX = "0x1.58c21da13dab6p-5"   # 0.04208474909858266
PINNED_PORTABLE_SHA256 = "9718f9e875ca0b5bed6106761073b4fe8fc94b4a4c0064a98d0592cdf40e0a59"


class CrossProcessReproducibilityTests(unittest.TestCase):
    runs: list[dict] = []

    @classmethod
    def setUpClass(cls) -> None:
        cls.runs = []
        for _ in range(2):
            completed = subprocess.run([sys.executable, "-B", "-c", CROSS_PROCESS_SCRIPT],
                                       capture_output=True, text=True, timeout=300, check=False)
            if completed.returncode != 0:
                raise AssertionError(f"subprocess failed:\n{completed.stderr}")
            cls.runs.append(json.loads(completed.stdout.strip().splitlines()[-1]))

    def test_two_processes_produce_identical_records(self) -> None:
        first, second = self.runs
        self.assertEqual(first, second)
        self.assertEqual(len(first["full"]), 64)

    def test_matches_in_process_run(self) -> None:
        record = run(rosenbrock_benchmark(2, -2.048, 4.096).problem(), 600, seed=2026).to_dict()
        record.pop("wall_time_s")
        record.pop("environment")
        digest = hashlib.sha256(json.dumps(record, sort_keys=True, allow_nan=False).encode()).hexdigest()
        self.assertEqual(digest, self.runs[0]["portable"])

    def test_pinned_reference(self) -> None:
        major_minor = tuple(int(part) for part in np.__version__.split(".")[:2])
        if major_minor != PINNED_NUMPY:
            self.skipTest(f"pin recorded with NumPy {PINNED_NUMPY[0]}.{PINNED_NUMPY[1]}.x; running {np.__version__}")
        self.assertEqual(self.runs[0]["f_best"], PINNED_F_BEST_HEX)
        self.assertEqual(self.runs[0]["portable"], PINNED_PORTABLE_SHA256)


class ResultTests(unittest.TestCase):
    def test_result_metadata_and_immutability(self) -> None:
        bounds = Bounds([-2.0, -1.0], [3.0, 4.0], names=("a", "b"), units=("m", "s"))
        problem = RecordingProblem(sphere, bounds)
        result = run(problem, 300, seed=8, record_population=True)
        self.assertEqual(result.algorithm, ALGORITHM)
        self.assertEqual(result.seed, 8)
        self.assertEqual(result.settings["population_size"], 20)
        self.assertEqual(result.problem_info["bounds"]["units"], ["m", "s"])
        self.assertEqual(result.environment["numpy"], np.__version__)
        self.assertEqual(result.environment["bit_generator"], "PCG64")
        assert_array_equal(result.x_best, bounds.from_unit(result.u_best))
        points = np.array(problem.points)
        values = np.array([sphere(p) for p in points])
        first_best = int(np.argmin(values))
        self.assertEqual(result.f_best, values[first_best])
        assert_array_equal(result.x_best, points[first_best])
        self.assertTrue(np.all(np.diff(result.history.best_so_far) <= 0))
        self.assertTrue(np.all(np.diff(result.history.generation_best) <= 0))
        populations = result.history.populations
        self.assertEqual(populations.shape, (result.generations + 1, 20, 2))
        assert_array_equal(bounds.from_unit(populations[0]), points[:20])
        for array in (result.x_best, result.u_best, result.history.best_so_far, populations):
            with self.assertRaises(ValueError):
                array[0] = 0
        with self.assertRaises(TypeError):
            result.settings["mutation"] = 1.0  # type: ignore[index]
        encoded = json.dumps(result.to_dict(), allow_nan=False)
        self.assertIn('"termination_reason": "max_evaluations"', encoded)
        self.assertNotIn("history", result.to_dict(include_history=False))
        assert_array_equal(bounds.lower, [-2.0, -1.0])


if __name__ == "__main__":
    unittest.main()
