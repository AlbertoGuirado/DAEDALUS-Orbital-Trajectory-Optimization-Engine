"""Benchmark infrastructure (Phase 3 increment c): registry, runner, report and experiment.

The runner must add nothing to a run: every record is compared with a direct
call of the optimizer. Aggregates are recounted from the records. Budgets in
these tests are small; they exercise the machinery, not optimizer performance.
"""
import csv
import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
from numpy.testing import assert_array_equal

from experiments.benchmark_optimizers import check_reproducibility, headline, run, scenario_table
from orbital_optimizer import bielliptic_transfer, hohmann_transfer
from orbital_optimizer.benchmarking import (
    SCHEMA, Scenario, build_optimizer, build_scenario, evaluate_acceptance, record_run, records_to_rows,
    run_benchmark, run_once, summarize_benchmark, validate_config,
)
from orbital_optimizer.benchmarking.report import format_number
from orbital_optimizer.optimization import (
    DifferentialEvolution, HybridOptimizer, NelderMead, ParticleSwarm, RandomSearch, TerminationCriteria,
    sphere_benchmark,
)
from orbital_optimizer.problems import CoplanarMultiImpulseProblem

REPOSITORY = Path(__file__).resolve().parents[1]
LEO, GEO = 6_678_000.0, 42_164_000.0

SPHERE = {"id": "sphere-2d", "type": "function", "function": "sphere", "dimension": 2, "lower": -2.56,
          "upper": 5.12, "max_evaluations": 400, "success_threshold": 1e-3}
TRANSFER = {"id": "leo-geo", "type": "coplanar_transfer", "initial_radius_m": LEO, "final_radius_m": GEO,
            "max_evaluations": 300, "success_threshold_m_s": 100.0}


def config(algorithms=None, scenarios=None, count: int = 3, **extra) -> dict:
    return {
        "schema": SCHEMA, "name": "test", "seeds": {"first": 5, "count": count},
        "algorithms": algorithms or [{"id": "de", "type": "differential_evolution"},
                                     {"id": "rs", "type": "random_search"}],
        "scenarios": scenarios or [SPHERE, TRANSFER],
        **extra,
    }


class RegistryTests(unittest.TestCase):
    def test_optimizers_are_built_with_their_settings(self) -> None:
        cases = {
            "differential_evolution": ({"mutation": 0.7}, "DE/rand/1/bin@1"),
            "particle_swarm": ({"swarm_size": 12}, "PSO/gbest/inertia@1"),
            "genetic_algorithm": ({"population_size": 10}, "GA/tournament2/SBX/polynomial/elitist@1"),
            "random_search": ({"batch_size": 5}, "RandomSearch/uniform@1"),
            "nelder_mead": ({"initial_point": [0.2, 0.3], "initial_step": 0.1}, "NelderMead/Lagarias1998@1"),
        }
        for kind, (settings, name) in cases.items():
            with self.subTest(kind=kind):
                optimizer = build_optimizer({"type": kind, "settings": settings})
                self.assertEqual(optimizer.name, name)
                for key, value in settings.items():
                    expected = tuple(value) if isinstance(value, list) else value
                    self.assertEqual(getattr(optimizer.settings, key), expected)
        self.assertIsInstance(build_optimizer({"type": "differential_evolution"}), DifferentialEvolution)
        hybrid = build_optimizer({"type": "hybrid", "settings": {"global": {"type": "particle_swarm"},
                                                                 "global_fraction": 0.6}})
        self.assertIsInstance(hybrid, HybridOptimizer)
        self.assertIsInstance(hybrid.global_optimizer, ParticleSwarm)
        self.assertEqual(hybrid.settings.global_fraction, 0.6)
        self.assertIsInstance(build_optimizer({"type": "hybrid"}).global_optimizer, DifferentialEvolution)

    def test_invalid_optimizer_configurations(self) -> None:
        for bad in ({"type": "simulated_annealing"}, {"settings": {}}, {"type": "random_search", "extra": 1},
                    {"type": "differential_evolution", "settings": {"mutaton": 0.5}},
                    {"type": "differential_evolution", "settings": {"mutation": 9.0}},
                    {"type": "hybrid", "settings": {"local": {}}}, {"type": "random_search", "settings": [1]},
                    "differential_evolution"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                build_optimizer(bad)  # type: ignore[arg-type]

    def test_function_scenario(self) -> None:
        scenario = build_scenario(SPHERE)
        self.assertIsInstance(scenario, Scenario)
        self.assertEqual((scenario.id, scenario.max_evaluations, scenario.objective_unit), ("sphere-2d", 400, "dimensionless"))
        self.assertEqual((scenario.reference, scenario.success_threshold, scenario.lower_threshold), (0.0, 1e-3, 0.0))
        self.assertIsNone(scenario.verifier)
        first, second = scenario.problem_factory(), scenario.problem_factory()
        self.assertIsNot(first, second)
        self.assertEqual(first.bounds.lower.tolist(), [-2.56, -2.56])

    def test_transfer_scenario_references(self) -> None:
        scenario = build_scenario(TRANSFER)
        self.assertEqual(scenario.reference, hohmann_transfer(LEO, GEO).total_delta_v)
        self.assertEqual((scenario.objective_unit, scenario.lower_threshold), ("m/s", 1e-6))
        self.assertIsInstance(scenario.problem_factory(), CoplanarMultiImpulseProblem)
        self.assertIsNotNone(scenario.verifier)
        self.assertIn("delta_v_total_m_s", scenario.metrics)
        three = build_scenario({"id": "x", "type": "coplanar_transfer", "initial_radius_m": LEO,
                                "final_radius_m": 20 * LEO, "impulses": 3, "max_evaluations": 100,
                                "reference": {"type": "bielliptic", "intermediate_radius_m": 60 * 20 * LEO}})
        self.assertEqual(three.reference, bielliptic_transfer(LEO, 20 * LEO, 1200 * LEO).total_delta_v)
        self.assertEqual(three.problem_factory().impulses, 3)
        self.assertIsNone(three.success_threshold)

    def test_invalid_scenario_configurations(self) -> None:
        for change in ({"type": "unknown"}, {"id": ""}, {"function": "ackley"}, {"dimension": 0},
                       {"max_evaluations": 0}, {"max_evaluations": 10.0}, {"lower": "a"}, {"budget": 5},
                       {"success_threshold": math.inf}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                build_scenario(SPHERE | change)
        for change in ({"final_radius_m": LEO}, {"impulses": 1}, {"reference": {"type": "lambert"}},
                       {"reference": {"type": "bielliptic"}}, {"verification": {"tolerance": 1.0}},
                       {"initial_radius_m": -1.0}, {"mu_m3_s2": 0.0}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                build_scenario(TRANSFER | change)
        with self.assertRaises(ValueError):
            build_scenario({"id": "missing-type"})


class RunnerTests(unittest.TestCase):
    def test_configuration_validation(self) -> None:
        scenarios, algorithms, seeds = validate_config(config())
        self.assertEqual(([s["id"] for s in scenarios], [a["id"] for a in algorithms], seeds),
                         (["sphere-2d", "leo-geo"], ["de", "rs"], [5, 6, 7]))
        base = config()
        for change in ({"schema": "daedalus.benchmark/0"}, {"seeds": {"first": 0}}, {"seeds": {"first": -1, "count": 2}},
                       {"seeds": {"first": 0, "count": 0}}, {"seeds": {"first": 0, "count": 2.0}},
                       {"algorithms": []}, {"scenarios": []}, {"algorithms": "de"}, {"budget": 100},
                       {"algorithms": [{"type": "random_search"}]},
                       {"algorithms": [{"id": "a", "type": "random_search"}, {"id": "a", "type": "random_search"}]},
                       {"scenarios": [SPHERE, SPHERE]}, {"algorithms": [{"id": "a", "type": "nope"}]}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate_config(base | change)
        without_name = dict(base)
        del without_name["name"]
        with self.assertRaises(ValueError):
            validate_config(without_name)
        with self.assertRaises(ValueError):
            validate_config("config")  # type: ignore[arg-type]

    def test_records_equal_direct_optimizer_calls_with_equal_budgets_and_seeds(self) -> None:
        records = run_benchmark(config())
        self.assertEqual([(r.scenario, r.algorithm, r.seed) for r in records],
                         [(s, a, seed) for s in ("sphere-2d", "leo-geo") for a in ("de", "rs") for seed in (5, 6, 7)])
        optimizers = {"de": DifferentialEvolution(), "rs": RandomSearch()}
        problems = {"sphere-2d": lambda: sphere_benchmark(2, -2.56, 5.12).problem(),
                    "leo-geo": lambda: CoplanarMultiImpulseProblem(LEO, GEO)}
        budgets = {"sphere-2d": 400, "leo-geo": 300}
        for record in records:
            with self.subTest(scenario=record.scenario, algorithm=record.algorithm, seed=record.seed):
                self.assertEqual(record.nfev, budgets[record.scenario])           # RF-6A.03
                self.assertEqual(record.termination_reason, "max_evaluations")
                direct = optimizers[record.algorithm].minimize(
                    problems[record.scenario](), TerminationCriteria(budgets[record.scenario]), seed=record.seed)
                self.assertEqual(record.f_best, direct.f_best)
                assert_array_equal(record.best_so_far, direct.history.best_so_far)
                self.assertEqual(record.x_best, tuple(direct.x_best))
                self.assertEqual(record.algorithm_name, optimizers[record.algorithm].name)
                with self.assertRaises(ValueError):
                    record.best_so_far[0] = 0.0

    def test_worker_processes_give_the_same_records(self) -> None:
        serial = run_benchmark(config(count=2))
        calls: list[tuple[int, int]] = []
        parallel = run_benchmark(config(count=2), jobs=2, progress=lambda done, total, record: calls.append((done, total)))
        self.assertEqual(calls, [(index, len(serial)) for index in range(1, len(serial) + 1)])
        for first, second in zip(serial, parallel):
            left, right = first.to_row(), second.to_row()
            left.pop("wall_time_s"), right.pop("wall_time_s")
            self.assertEqual(left, right)
            self.assertEqual(first.best_so_far.tobytes(), second.best_so_far.tobytes())
        for jobs in (0, True, 1.5):
            with self.subTest(jobs=jobs), self.assertRaises(ValueError):
                run_benchmark(config(), jobs=jobs)  # type: ignore[arg-type]

    def test_success_gap_and_anomaly_follow_the_scenario_criteria(self) -> None:
        scenario = build_scenario(SPHERE)
        record = run_once(scenario, "de", DifferentialEvolution(), 0)
        self.assertEqual(record.gap, record.f_best - 0.0)
        self.assertEqual(record.success, record.f_best <= 1e-3)
        self.assertFalse(record.anomaly)
        if record.success:
            first = int(np.nonzero(record.best_so_far <= 1e-3)[0][0]) + 1
            self.assertEqual(record.evaluations_to_success, first)
        # A reference above what the optimizer reaches: every run is "better than the optimum".
        impossible = Scenario("sphere", "wrong reference", scenario.problem_factory, 400, "dimensionless",
                              reference=5.0, success_threshold=0.0, lower_threshold=0.0)
        flagged = run_once(impossible, "de", DifferentialEvolution(), 0)
        self.assertTrue(flagged.anomaly)
        self.assertLess(flagged.gap, 0.0)
        no_criteria = Scenario("sphere", "no reference", scenario.problem_factory, 400, "dimensionless")
        plain = run_once(no_criteria, "de", DifferentialEvolution(), 0)
        self.assertEqual((plain.gap, plain.success, plain.anomaly, plain.evaluations_to_success), (None, None, False, None))

    def test_transfer_runs_are_verified_and_report_metrics(self) -> None:
        scenario = build_scenario(TRANSFER)
        record = run_once(scenario, "de", DifferentialEvolution(), 1)
        self.assertTrue(record.verification["passed"], record.verification["message"])
        self.assertEqual(record.verification["tolerances"], {"radius_m": 1.0, "eccentricity": 1e-6})
        self.assertEqual((record.verification["captured_value_errors"], record.verification["captured_runtime_errors"]), (0, 0))
        self.assertAlmostEqual(record.metrics["delta_v_total_m_s"], record.f_best, places=9)
        self.assertGreater(record.metrics["transfer_time_s"], 0.0)
        self.assertEqual(record.to_row()["metric_delta_v_total_m_s"], record.metrics["delta_v_total_m_s"])
        self.assertTrue(record.to_row()["verification_passed"])
        # A verifier that fails turns a numerically good run into a failure.
        problem = scenario.problem_factory()
        result = DifferentialEvolution().minimize(problem, TerminationCriteria(300), seed=1)
        failing = Scenario(scenario.id, "", scenario.problem_factory, 300, "m/s", reference=scenario.reference,
                           success_threshold=1e9, lower_threshold=1e-6,
                           verifier=lambda problem, result: {"passed": False, "message": "forced"})
        self.assertFalse(record_run(failing, "de", DifferentialEvolution(), 1, problem, result).success)
        passing = Scenario(scenario.id, "", scenario.problem_factory, 300, "m/s", reference=scenario.reference,
                           success_threshold=1e9, lower_threshold=1e-6, verifier=scenario.verifier)
        self.assertTrue(record_run(passing, "de", DifferentialEvolution(), 1, problem, result).success)

    def test_run_without_a_valid_point(self) -> None:
        scenario = Scenario("local", "", build_scenario(TRANSFER).problem_factory, 60, "m/s",
                            reference=hohmann_transfer(LEO, GEO).total_delta_v, success_threshold=1.0,
                            lower_threshold=1e-6, verifier=build_scenario(TRANSFER).verifier)
        # dv1 = 0: the vehicle never leaves the initial orbit; a local search started there finds nothing valid.
        optimizer = build_optimizer({"type": "nelder_mead", "settings": {"initial_point": [0.0, 0.5], "initial_step": 0.01}})
        self.assertIsInstance(optimizer, NelderMead)
        record = run_once(scenario, "nm", optimizer, 0)
        self.assertFalse(record.found_valid)
        self.assertEqual((record.f_best, record.gap), (math.inf, math.inf))
        self.assertEqual((record.success, record.anomaly, record.verification, record.x_best), (False, False, None, None))
        self.assertEqual(record.metrics, {})


class ReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = config(algorithms=[{"id": "de", "type": "differential_evolution"},
                                        {"id": "rs", "type": "random_search"},
                                        {"id": "pso", "type": "particle_swarm"}], count=5)
        cls.records = run_benchmark(cls.config)
        cls.summary = summarize_benchmark(cls.config, cls.records)

    def runs(self, scenario: str, algorithm: str) -> list:
        return [r for r in self.records if r.scenario == scenario and r.algorithm == algorithm]

    def test_cells_recount_the_records(self) -> None:
        for scenario in ("sphere-2d", "leo-geo"):
            for algorithm in ("de", "rs", "pso"):
                with self.subTest(scenario=scenario, algorithm=algorithm):
                    runs, cell = self.runs(scenario, algorithm), self.summary["scenarios"][scenario][algorithm]
                    finals = sorted(run.f_best for run in runs)
                    self.assertEqual((cell["runs"], cell["seeds"]), (5, [5, 6, 7, 8, 9]))
                    self.assertEqual((cell["objective"]["best"], cell["objective"]["median"], cell["objective"]["worst"]),
                                     (finals[0], finals[2], finals[4]))
                    self.assertEqual(cell["gap"]["median"], sorted(run.gap for run in runs)[2])
                    successes = sum(run.success for run in runs)
                    self.assertEqual(cell["success"]["successes"], successes)
                    self.assertEqual(cell["success"]["rate"], successes / 5)
                    low, high = cell["success"]["wilson_95"]
                    self.assertTrue(low <= successes / 5 <= high)
                    self.assertEqual(cell["termination_reasons"], {"max_evaluations": 5})
                    self.assertEqual(cell["anomalies"], 0)
                    band = cell["convergence"]
                    self.assertEqual(band["evaluations"][-1], runs[0].nfev)
                    self.assertEqual(band["median"][-1], finals[2])            # history ends at the final value
                    # Non-increasing; the history is +inf until the first valid evaluation.
                    self.assertTrue(all(later <= earlier for earlier, later in zip(band["median"], band["median"][1:])))
                    self.assertTrue(all(a <= b <= c for a, b, c in zip(band["q25"], band["median"], band["q75"])))
                    self.assertEqual(cell["found_valid"], sum(run.found_valid for run in runs))
        transfer = self.summary["scenarios"]["leo-geo"]["de"]
        # Only runs that found a valid point have a plan to verify; every such plan must pass.
        valid_runs = sum(run.found_valid for run in self.runs("leo-geo", "de"))
        self.assertGreater(valid_runs, 0)
        self.assertEqual(transfer["verification"], {"checked": valid_runs, "passed": valid_runs})
        self.assertIn("delta_v_total_m_s", transfer["metrics"])
        self.assertIsNone(self.summary["scenarios"]["sphere-2d"]["de"]["verification"])

    def test_pairwise_comparisons_and_overall_test(self) -> None:
        for scenario in ("sphere-2d", "leo-geo"):
            comparisons = self.summary["comparisons"][scenario]
            self.assertEqual([(c["first"], c["second"]) for c in comparisons], [("de", "rs"), ("de", "pso"), ("rs", "pso")])
            for comparison in comparisons:
                self.assertGreaterEqual(comparison["p_value_holm"], comparison["p_value"])
                self.assertTrue(0.0 <= comparison["a12"] <= 1.0)
        # DE against random search on the sphere: all five DE runs are better (A12 = 1, exact p = 2/252).
        de_rs = self.summary["comparisons"]["sphere-2d"][0]
        self.assertEqual(de_rs["a12"], 1.0)
        self.assertAlmostEqual(de_rs["p_value"], 2 / math.comb(10, 5), places=12)
        overall = self.summary["overall"]
        self.assertEqual((overall["test"], overall["algorithms"], overall["n_problems"]), ("friedman", ["de", "rs", "pso"], 2))
        self.assertEqual(len(overall["average_ranks"]), 3)

    def test_two_algorithms_use_the_signed_rank_test_and_one_has_no_overall_test(self) -> None:
        two = config(count=2)
        self.assertEqual(summarize_benchmark(two, run_benchmark(two))["overall"]["test"], "wilcoxon-signed-rank")
        one = config(algorithms=[{"id": "rs", "type": "random_search"}], count=2)
        summary = summarize_benchmark(one, run_benchmark(one))
        self.assertIsNone(summary["overall"])
        self.assertEqual(summary["comparisons"]["sphere-2d"], [])

    def test_records_must_match_the_configuration(self) -> None:
        with self.assertRaises(ValueError):
            summarize_benchmark(self.config, self.records[:-1])
        with self.assertRaises(ValueError):
            summarize_benchmark(config(count=5), self.records)

    def test_acceptance_criteria(self) -> None:
        vacuous = evaluate_acceptance(None, self.summary)
        self.assertTrue(vacuous["passed"])
        self.assertEqual(vacuous["note"], "no acceptance criteria declared")
        met = evaluate_acceptance({"max_anomalies": 0, "require_verification": True,
                                   "minimum_success_rate": {"sphere-2d": {"de": 0.2}}}, self.summary)
        self.assertEqual(met["checks"], {"anomalies": True, "verification": True, "success_rate:sphere-2d/de": True})
        self.assertTrue(met["passed"])
        observed = self.summary["scenarios"]["sphere-2d"]["rs"]["success"]["rate"]
        self.assertLess(observed, 1.0)
        unmet = evaluate_acceptance({"minimum_success_rate": {"sphere-2d": {"rs": 1.0}}}, self.summary)
        self.assertFalse(unmet["passed"])
        self.assertEqual(unmet["details"]["success_rate:sphere-2d/rs"]["observed"], observed)
        with self.assertRaises(ValueError):
            evaluate_acceptance({"minimum_gap": 1.0}, self.summary)
        with self.assertRaises(ValueError):
            evaluate_acceptance({"minimum_success_rate": {"sphere-2d": {"ga": 0.5}}}, self.summary)

    def test_rows_and_number_formatting(self) -> None:
        rows = records_to_rows(self.records)
        self.assertEqual(len(rows), len(self.records))
        self.assertTrue(all(list(row) == list(rows[0]) for row in rows))      # one stable set of columns
        self.assertIn("metric_delta_v_total_m_s", rows[0])
        self.assertIsNone(rows[0]["metric_delta_v_total_m_s"])                # the sphere has no such metric
        self.assertEqual([format_number(v) for v in (None, True, "x", math.nan, math.inf, -math.inf, 1234.5678)],
                         ["-", "True", "x", "nan", "inf", "-inf", "1235"])


class ExperimentTests(unittest.TestCase):
    def test_artifacts_reproducibility_and_acceptance(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            report = run(config(count=2, acceptance={"max_anomalies": 0, "require_verification": True}),
                         output, quiet=True)
            self.assertTrue(report["passed"])
            self.assertEqual(report["reproducibility"], report["reproducibility"] | {"checked": 4, "identical": 4, "passed": True})
            stored = json.loads((output / "report.json").read_text(encoding="utf-8"),
                                parse_constant=lambda token: self.fail(f"non-strict JSON token {token}"))
            self.assertEqual(stored["configuration"]["seeds"], {"first": 5, "count": 2})
            self.assertEqual(stored["execution"]["runs"], 8)
            self.assertEqual(set(stored["environment"]), {"python", "python_implementation", "numpy", "scipy", "platform", "machine"})
            self.assertIn("git_commit", stored["code"])
            self.assertEqual(stored["scenarios"]["leo-geo"]["reference"], hohmann_transfer(LEO, GEO).total_delta_v)
            self.assertTrue(stored["scenarios"]["leo-geo"]["verified_by_propagation"])
            with (output / "runs.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 8)
            self.assertEqual(json.loads(rows[0]["x_best"]), list(json.loads(rows[0]["x_best"])))
            with (output / "convergence.csv").open(encoding="utf-8", newline="") as handle:
                band = list(csv.DictReader(handle))
            self.assertEqual({row["scenario"] for row in band}, {"sphere-2d", "leo-geo"})
            self.assertEqual({row["objective_unit"] for row in band if row["scenario"] == "leo-geo"}, {"m/s"})
            self.assertEqual(set(headline(report)), {"sphere-2d", "leo-geo"})

    def test_unmet_acceptance_fails_the_experiment(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            report = run(config(count=2, acceptance={"minimum_success_rate": {"sphere-2d": {"rs": 1.0}}}),
                         Path(folder) / "run", quiet=True)
            self.assertFalse(report["passed"])
            self.assertFalse(report["acceptance"]["checks"]["success_rate:sphere-2d/rs"])
            self.assertTrue(report["reproducibility"]["passed"])

    def test_reproducibility_check_detects_a_different_run(self) -> None:
        settings = config(count=1)
        records = run_benchmark(settings)
        self.assertTrue(check_reproducibility(settings, records)["passed"])
        other = run_benchmark(settings | {"seeds": {"first": 6, "count": 1}})
        tampered = [type(record)(**{**record.__dict__, "seed": 5}) for record in other]
        result = check_reproducibility(settings, tampered)
        self.assertFalse(result["passed"])
        self.assertEqual(len(result["different"]), 4)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib is optional")
    def test_figures_are_exported(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            report = run(config(count=2), output, plots=True, quiet=True)
            for name in ("convergence.png", "final_distribution.png", "success_rate.png", "cost.png",
                         "transfer_leo-geo.png"):
                with self.subTest(figure=name):
                    self.assertIn(name, report["artifacts"])
                    self.assertGreater((output / name).stat().st_size, 10_000)

    def test_versioned_benchmark_configurations_are_valid(self) -> None:
        paths = sorted((REPOSITORY / "configs" / "benchmarks").glob("*.json"))
        self.assertTrue(paths)
        for path in paths:
            with self.subTest(config=path.name):
                settings = json.loads(path.read_text(encoding="utf-8"))
                scenarios, algorithms, seeds = validate_config(settings)
                self.assertGreaterEqual(len(seeds), 30)                       # RF-6A.02
                table = scenario_table(settings)
                self.assertTrue(all(row["reference"] is not None for row in table.values()))
                evaluate_acceptance(settings.get("acceptance"), {"scenarios": {
                    scenario["id"]: {algorithm["id"]: {"anomalies": 0, "verification": None,
                                                       "success": {"rate": 1.0, "wilson_95": [1.0, 1.0]}}
                                     for algorithm in algorithms} for scenario in scenarios}})


if __name__ == "__main__":
    unittest.main()
