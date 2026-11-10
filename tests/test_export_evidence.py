"""Evidence export: the versioned record must be derived from a report, unchanged except for thinning."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from experiments.export_evidence import SCHEMA, evidence_record, thin


class ThinningTests(unittest.TestCase):
    def test_long_series_are_thinned_together_and_keep_their_ends(self) -> None:
        band = {"evaluations": list(range(1, 101)), "median": [float(100 - k) for k in range(100)],
                "note": "text", "short": [1, 2, 3]}
        result = thin(band, 5)
        # Positions round(linspace(0, 99, 5)) = 0, 25, 50, 74, 99 (49.5 rounds to the even 50).
        self.assertEqual(result["evaluations"], [1, 26, 51, 75, 100])
        self.assertEqual(result["median"], [100.0, 75.0, 50.0, 26.0, 1.0])     # same positions as the counts
        self.assertEqual((result["note"], result["short"]), ("text", [1, 2, 3]))
        self.assertEqual(len(band["evaluations"]), 100)                         # the input is not modified

    def test_structure_and_scalars_are_preserved(self) -> None:
        report = {"a": {"b": [{"c": list(range(30)), "d": None}], "flag": True}, "names": ["x", "y"], "value": 1.5}
        result = thin(report, 10)
        self.assertEqual(result["a"]["b"][0]["c"], [0, 3, 6, 10, 13, 16, 19, 23, 26, 29])
        self.assertEqual((result["a"]["b"][0]["d"], result["a"]["flag"], result["names"], result["value"]),
                         (None, True, ["x", "y"], 1.5))
        mixed = {"series": list(range(20)), "other": list(range(25))}
        thinned = thin(mixed, 10)
        self.assertEqual((len(thinned["series"]), len(thinned["other"])), (10, 10))
        self.assertEqual((thinned["series"][-1], thinned["other"][-1]), (19, 24))


class EvidenceRecordTests(unittest.TestCase):
    def test_record_identifies_its_source(self) -> None:
        report = {"experiment": "optimizer_benchmark", "created_at_utc": "2026-11-10T12:00:00+00:00",
                  "passed": True, "summary": {"band": {"evaluations": list(range(200)), "median": [1.0] * 200}}}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "report.json"
            path.write_text(json.dumps(report), encoding="utf-8")
            record = evidence_record(path, 20)
            self.assertEqual(record["schema"], SCHEMA)
            self.assertEqual(record["source"]["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual((record["source"]["file"], record["source"]["experiment"]),
                             ("report.json", "optimizer_benchmark"))
            self.assertEqual(record["source"]["created_at_utc"], "2026-11-10T12:00:00+00:00")
            self.assertEqual(record["report"]["passed"], True)
            self.assertEqual(len(record["report"]["summary"]["band"]["evaluations"]), 20)
            self.assertIn("recorded execution", record["kind"])
            json.dumps(record, allow_nan=False)
            for bad in (1, True, 2.5):
                with self.subTest(max_points=bad), self.assertRaises(ValueError):
                    evidence_record(path, bad)  # type: ignore[arg-type]
            path.write_text(json.dumps({"not": "a report"}), encoding="utf-8")
            with self.assertRaises(ValueError):
                evidence_record(path)

    def test_versioned_evidence_records_are_well_formed(self) -> None:
        folder = Path(__file__).resolve().parents[1] / "docs" / "evidence"
        for path in sorted(folder.glob("*.json")):
            with self.subTest(record=path.name):
                record = json.loads(path.read_text(encoding="utf-8"),
                                    parse_constant=lambda token: self.fail(f"non-strict JSON token {token}"))
                self.assertEqual(record["schema"], SCHEMA)
                self.assertRegex(record["source"]["sha256"], r"^[0-9a-f]{64}$")
                self.assertIn("created_at_utc", record["report"])
                self.assertIn("environment", record["report"])
                self.assertIn("configuration", record["report"])


if __name__ == "__main__":
    unittest.main()
