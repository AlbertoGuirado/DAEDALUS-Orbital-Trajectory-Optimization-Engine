"""Reproducible SI orbital validation experiment."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
from time import perf_counter
from typing import Any

import numpy as np
import scipy
from orbital_optimizer import (
    OrbitalElements, elements_to_state, orbital_period, propagate,
    specific_energy, specific_angular_momentum, eccentricity_vector, PropagationResult,
)
from orbital_optimizer.dynamics.two_body import real_array, validate_mu

METRICS = (
    "max_position_error_m", "max_velocity_error_m_s", "max_relative_energy_drift",
    "max_relative_angular_momentum_drift", "max_absolute_eccentricity_drift",
)


def run_experiment(config: dict[str, Any]) -> tuple[dict[str, Any], PropagationResult, np.ndarray]:
    """Propagate a configured SI circular orbit and compare r(t), v(t) analytically.

    Reference: r(t)=r0*cos(nt)+v0*sin(nt)/n; v(t)=-n*r0*sin(nt)+v0*cos(nt).
    n=sqrt(mu/radius^3), in rad/s. Return JSON-compatible metrics, sampled SI
    trajectory, and analytic (N,6) states. Thresholds are explicit in config.
    """
    keys = {"radius_m", "inclination_rad", "mu_m3_s2", "periods", "samples",
            "method", "max_step_s", "rtol", "atol", "acceptance"}
    if not isinstance(config, dict) or set(config) != keys:
        raise ValueError(f"configuration must contain exactly these keys: {sorted(keys)}")
    for name, minimum in (("periods", 1), ("samples", 2)):
        if type(config[name]) is not int or config[name] < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    bounds = config["acceptance"]
    if not isinstance(bounds, dict) or set(bounds) != set(METRICS):
        raise ValueError(f"acceptance must contain exactly {list(METRICS)}")
    for name, value in bounds.items():
        scalar = real_array(value, name)
        if scalar.shape != () or scalar <= 0:
            raise ValueError(f"{name} must be a positive finite scalar")
    gm = validate_mu(config["mu_m3_s2"])
    orbit = OrbitalElements(config["radius_m"], 0., config["inclination_rad"])
    initial = elements_to_state(orbit, gm)
    period = orbital_period(orbit.semi_major_axis, gm)
    duration = period*config["periods"]
    samples = np.linspace(0., duration, config["samples"])
    start = perf_counter()
    trajectory = propagate(
        initial, duration, gm, method=config["method"], sample_times=samples,
        max_step=config["max_step_s"], rtol=config["rtol"], atol=config["atol"],
    )
    elapsed = perf_counter()-start
    n = 2*np.pi/period
    phase = trajectory.times*n
    cosine, sine = np.cos(phase)[:,None], np.sin(phase)[:,None]
    reference = np.column_stack((initial[:3]*cosine+initial[3:]*sine/n,
                                 -n*initial[:3]*sine+initial[3:]*cosine))
    energy = np.array([specific_energy(s,gm) for s in trajectory.states])
    h = np.array([specific_angular_momentum(s) for s in trajectory.states])
    e = np.array([eccentricity_vector(s,gm) for s in trajectory.states])
    metrics = {
        "max_position_error_m": float(np.max(np.linalg.norm(trajectory.states[:,:3]-reference[:,:3],axis=1))),
        "max_velocity_error_m_s": float(np.max(np.linalg.norm(trajectory.states[:,3:]-reference[:,3:],axis=1))),
        "max_relative_energy_drift": float(np.max(np.abs(energy-energy[0]))/abs(energy[0])),
        "max_relative_angular_momentum_drift": float(np.max(np.linalg.norm(h-h[0],axis=1))/np.linalg.norm(h[0])),
        "max_absolute_eccentricity_drift": float(np.max(np.linalg.norm(e-e[0],axis=1))),
    }
    passed = {name: metrics[name] <= bounds[name] for name in METRICS}
    report = {
        "experiment": "phase1_circular_orbit", "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "configuration": config,
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "platform": platform.platform()},
        "model": "Newtonian two-body point mass; no perturbations or burns",
        "randomness": "none; deterministic experiment, no seed required",
        "period_s": period, "propagation_duration_s": duration,
        "initial_state_si": initial.tolist(), "final_state_si": trajectory.final_state.tolist(),
        "computation_time_s": elapsed, "function_evaluations": trajectory.nfev,
        "metrics": metrics, "checks": passed, "validation_passed": all(passed.values()),
        "limitations": ["Errors are maxima over saved samples, not continuous-time bounds.",
                        "Computation time covers one propagation run, not a statistical benchmark.",
                        "Integrator tolerances are local; acceptance compares global sampled errors."],
    }
    return report, trajectory, reference


def main() -> int:
    """CLI: write JSON/CSV (SI) and optional PNGs; exit 1 on failed acceptance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/earth_orbit.json"))
    parser.add_argument("--output", type=Path, default=Path("outputs/earth_orbit"))
    parser.add_argument("--plots", action="store_true", help="Export scientific PNG figures with Matplotlib")
    args = parser.parse_args()
    try:
        config = json.loads(args.config.read_text(encoding="utf-8"))
        report, trajectory, reference = run_experiment(config)
        args.output.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.output/"trajectory.csv", np.column_stack((trajectory.times,trajectory.states)),
                   delimiter=",", header="time_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s", comments="")
        if args.plots:
            from orbital_optimizer.visualization.plots import save_validation_plots
            save_validation_plots(trajectory, reference, args.output)
        (args.output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        parser.exit(2, f"Experiment failed: {exc}\n")
    print(json.dumps({"validation_passed": report["validation_passed"], "metrics": report["metrics"],
                      "function_evaluations": report["function_evaluations"],
                      "report": str((args.output/"report.json").resolve())},indent=2))
    return 0 if report["validation_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
