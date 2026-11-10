"""Export a compact, versionable evidence record from an experiment report.

    python -m experiments.export_evidence --report outputs/<run>/report.json \
        --output docs/evidence/<name>.json

``outputs/`` holds regenerable artifacts and is not versioned (RR-07). The
documentation and the Verification Console need a small, versioned record of
what a recorded execution measured. This tool derives that record from the
report itself -- nothing is typed by hand -- and stores the SHA-256 of the
source report so the record can be matched to the full artifact.

The record keeps the configuration, environment, code version, statistics,
reproducibility and acceptance sections unchanged. Long numerical series
(lists of more than ``--max-points`` numbers, such as convergence bands) are
thinned to at most that many entries, always keeping the first and the last;
sibling series of equal length are thinned at the same positions, so a band
stays aligned with its evaluation counts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

SCHEMA = "daedalus.evidence/1"


def _is_series(value: Any) -> bool:
    return isinstance(value, list) and len(value) > 0 and all(
        isinstance(item, (int, float, str)) and not isinstance(item, bool) for item in value)


def thin(value: Any, max_points: int) -> Any:
    """Copy of a JSON value with every long numerical series reduced to ``max_points`` entries."""
    if isinstance(value, dict):
        lengths = {len(item) for item in value.values() if _is_series(item) and len(item) > max_points}
        result = {}
        for key, item in value.items():
            if _is_series(item) and len(item) in lengths:
                keep = np.unique(np.round(np.linspace(0, len(item) - 1, max_points)).astype(int))
                result[key] = [item[index] for index in keep]
            else:
                result[key] = thin(item, max_points)
        return result
    if isinstance(value, list):
        return [thin(item, max_points) for item in value]
    return value


def evidence_record(report_path: Path, max_points: int = 60) -> dict[str, Any]:
    """Evidence record of the report at ``report_path`` (see the module docstring)."""
    if isinstance(max_points, bool) or not isinstance(max_points, int) or max_points < 2:
        raise ValueError("max_points must be an integer >= 2")
    raw = report_path.read_bytes()
    report = json.loads(raw.decode("utf-8"))
    if not isinstance(report, dict) or "experiment" not in report or "created_at_utc" not in report:
        raise ValueError("the file is not an experiment report (missing 'experiment' or 'created_at_utc')")
    return {
        "schema": SCHEMA,
        "kind": "recorded execution: values measured when the experiment ran, not recomputed on reading",
        "exported_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {"file": report_path.name, "sha256": hashlib.sha256(raw).hexdigest(),
                   "experiment": report["experiment"], "created_at_utc": report["created_at_utc"]},
        "series_thinned_to": max_points,
        "report": thin(report, max_points),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Export a compact evidence record from an experiment report.")
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-points", type=int, default=60)
    args = parser.parse_args()
    try:
        record = evidence_record(args.report, args.max_points)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Export failed: {exc}\n")
    print(json.dumps({"output": str(args.output.resolve()), "bytes": args.output.stat().st_size,
                      "source_sha256": record["source"]["sha256"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
