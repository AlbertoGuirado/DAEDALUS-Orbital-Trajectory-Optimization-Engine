"""Verify reproducible settings and explicit experiment acceptance outcomes."""
import json
from pathlib import Path
import unittest
from experiments.validate_two_body import run_experiment


class ExperimentTests(unittest.TestCase):
    def setUp(self) -> None:
        path = Path(__file__).resolve().parents[1]/"configs"/"earth_orbit.json"
        self.config = json.loads(path.read_text(encoding="utf-8"))

    def test_reference_experiment_and_strict_failure(self) -> None:
        report, trajectory, _ = run_experiment(self.config)
        self.assertTrue(report["validation_passed"])
        self.assertEqual(len(trajectory.times),self.config["samples"])
        self.assertEqual(report["configuration"],self.config)
        json.dumps(report,allow_nan=False)
        self.config["acceptance"]["max_position_error_m"] = 1e-15
        failed, _, _ = run_experiment(self.config)
        self.assertFalse(failed["validation_passed"])
        self.assertFalse(failed["checks"]["max_position_error_m"])

    def test_bad_configuration_rejected(self) -> None:
        self.config["periods"] = 1.5
        with self.assertRaises(ValueError):
            run_experiment(self.config)
