"""Computations and figures for the Streamlit interface; no Streamlit calls here.

Physics stays in SI; km, minutes and degrees appear only at this presentation boundary.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import altair as alt
import numpy as np
import pandas as pd
import plotly.graph_objects as go

ROOT = Path(__file__).resolve().parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from experiments.validate_two_body import METRICS, run_experiment  # noqa: E402
from orbital_optimizer import (  # noqa: E402
    EARTH_MU, eccentricity_vector, hohmann_transfer, orbital_period, propagate,
    propagate_with_impulses, propellant_mass, specific_angular_momentum,
    specific_energy, state_to_elements,
)

# WGS84 equatorial radius; used only to convert altitudes and draw the Earth.
EARTH_RADIUS_M = 6_378_137.0

# Categorical slots 1-3 of the dataviz reference palette, plus a neutral reference.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
NEUTRAL = "#8a8986"

METRIC_LABELS = {
    "max_position_error_m": ("Error máximo de posición", "m"),
    "max_velocity_error_m_s": ("Error máximo de velocidad", "m/s"),
    "max_relative_energy_drift": ("Deriva relativa de energía", "-"),
    "max_relative_angular_momentum_drift": ("Deriva relativa de momento angular", "-"),
    "max_absolute_eccentricity_drift": ("Deriva absoluta de excentricidad", "-"),
}

PHASE2_LIMITS = {"position_m": 1.0, "velocity_m_s": 1e-3, "eccentricity": 1e-6}

TEST_SUITES = {
    "Todas": None,
    "Fase 1": ["tests.test_dynamics", "tests.test_elements", "tests.test_experiment", "tests.test_propagator"],
    "Fase 2": ["tests.test_transfers", "tests.test_impulses"],
}


# --- Phase 1 -------------------------------------------------------------

def run_phase1(config: dict[str, Any]) -> dict[str, Any]:
    """Run the official Phase 1 experiment and tabulate its time series."""
    report, trajectory, reference = run_experiment(config)
    states, mu = trajectory.states, trajectory.mu
    energy = np.array([specific_energy(s, mu) for s in states])
    h = np.array([specific_angular_momentum(s) for s in states])
    e = np.array([eccentricity_vector(s, mu) for s in states])
    series = pd.DataFrame({
        "time_min": trajectory.times / 60,
        "position_error_m": np.linalg.norm(states[:, :3] - reference[:, :3], axis=1),
        "velocity_error_m_s": np.linalg.norm(states[:, 3:] - reference[:, 3:], axis=1),
        "energy": np.abs(energy - energy[0]) / abs(energy[0]),
        "angular_momentum": np.linalg.norm(h - h[0], axis=1) / np.linalg.norm(h[0]),
        "eccentricity": np.linalg.norm(e - e[0], axis=1),
    })
    trajectory_table = pd.DataFrame(
        np.column_stack((trajectory.times, states)),
        columns=["time_s", "x_m", "y_m", "z_m", "vx_m_s", "vy_m_s", "vz_m_s"],
    )
    checks = pd.DataFrame([
        {
            "Métrica": METRIC_LABELS[name][0],
            "Valor": report["metrics"][name],
            "Límite": config["acceptance"][name],
            "Unidad": METRIC_LABELS[name][1],
            "Cumple": report["checks"][name],
        }
        for name in METRICS
    ])
    return {
        "report": report, "series": series, "trajectory": trajectory_table, "checks": checks,
        "orbit_km": states[:, :3] / 1000, "reference_km": reference[:, :3] / 1000,
    }


# --- Phase 2 -------------------------------------------------------------

def _circular_state(radius: float, mu: float) -> np.ndarray:
    return np.array([radius, 0.0, 0.0, 0.0, np.sqrt(mu / radius), 0.0])


def run_phase2(
    initial_altitude_km: float, final_altitude_km: float, initial_mass_kg: float,
    specific_impulse_s: float, method: str, max_step_s: float, rtol: float,
) -> dict[str, Any]:
    """Hohmann baseline, Tsiolkovsky budget and RV-2.01 check by propagation.

    Orbits are equatorial and counterclockwise; the transfer starts on +x.
    """
    mu = EARTH_MU
    r1 = EARTH_RADIUS_M + initial_altitude_km * 1000
    r2 = EARTH_RADIUS_M + final_altitude_km * 1000
    transfer = hohmann_transfer(r1, r2, mu)
    options = {"method": method, "max_step": max_step_s, "rtol": rtol}
    sign = 1.0 if r2 >= r1 else -1.0
    T = transfer.transfer_time

    start = _circular_state(r1, mu)
    first_burn = sign * transfer.first_delta_v * start[3:] / np.linalg.norm(start[3:])
    after_first = start + np.r_[0, 0, 0, first_burn]
    arc_times = np.linspace(0, T, 400)
    arc = propagate(after_first, T, mu, sample_times=arc_times, **options)
    apse = arc.final_state
    second_burn = sign * transfer.second_delta_v * apse[3:] / np.linalg.norm(apse[3:])

    final = propagate_with_impulses(start, [(0.0, first_burn), (T, second_burn)], T, mu, **options)
    elements = state_to_elements(final, mu)
    target_position = np.array([-r2, 0.0, 0.0])
    target_velocity = np.array([0.0, -np.sqrt(mu / r2), 0.0])
    position_error = float(np.linalg.norm(final[:3] - target_position))
    velocity_error = float(np.linalg.norm(final[3:] - target_velocity))
    checks = pd.DataFrame([
        {"Comprobación": "Posición final frente a la órbita objetivo", "Valor": position_error,
         "Límite": PHASE2_LIMITS["position_m"], "Unidad": "m",
         "Cumple": position_error <= PHASE2_LIMITS["position_m"]},
        {"Comprobación": "Velocidad final frente a la circular objetivo", "Valor": velocity_error,
         "Límite": PHASE2_LIMITS["velocity_m_s"], "Unidad": "m/s",
         "Cumple": velocity_error <= PHASE2_LIMITS["velocity_m_s"]},
        {"Comprobación": "Excentricidad de la órbita final", "Valor": elements.eccentricity,
         "Límite": PHASE2_LIMITS["eccentricity"], "Unidad": "-",
         "Cumple": elements.eccentricity <= PHASE2_LIMITS["eccentricity"]},
    ])

    first_propellant = propellant_mass(initial_mass_kg, transfer.first_delta_v, specific_impulse_s)
    second_propellant = propellant_mass(
        initial_mass_kg - first_propellant, transfer.second_delta_v, specific_impulse_s,
    )
    total_propellant = propellant_mass(initial_mass_kg, transfer.total_delta_v, specific_impulse_s)

    circle = np.linspace(0, 2 * np.pi, 361)
    arc_radius = np.linalg.norm(arc.states[:, :3], axis=1)
    arc_speed = np.linalg.norm(arc.states[:, 3:], axis=1)
    profile = pd.DataFrame({
        "time_min": arc.times / 60,
        "altitude_km": (arc_radius - EARTH_RADIUS_M) / 1000,
        "speed_km_s": arc_speed / 1000,
    })

    # Hohmann is symmetric in direction, so the sweep starts from the lower orbit.
    lower, ratio = min(r1, r2), max(r1, r2) / min(r1, r2)
    ratios = np.geomspace(1.0, max(40.0, 1.5 * ratio), 300)
    sweep = pd.DataFrame({
        "ratio": ratios,
        "total_delta_v_km_s": [hohmann_transfer(lower, lower * k, mu).total_delta_v / 1000 for k in ratios],
    })

    result = {
        "inputs": {
            "initial_altitude_km": initial_altitude_km, "final_altitude_km": final_altitude_km,
            "initial_radius_m": r1, "final_radius_m": r2, "initial_mass_kg": initial_mass_kg,
            "specific_impulse_s": specific_impulse_s, "mu_m3_s2": mu, **options,
        },
        "hohmann": {
            "first_delta_v_m_s": transfer.first_delta_v,
            "second_delta_v_m_s": transfer.second_delta_v,
            "total_delta_v_m_s": transfer.total_delta_v,
            "transfer_time_s": transfer.transfer_time,
            "transfer_semi_major_axis_m": transfer.transfer_semi_major_axis,
            "initial_period_s": orbital_period(r1, mu), "final_period_s": orbital_period(r2, mu),
        },
        "propellant": {
            "first_burn_kg": first_propellant, "second_burn_kg": second_propellant,
            "total_kg": total_propellant, "final_mass_kg": initial_mass_kg - total_propellant,
        },
        "verification": {
            "final_state_si": final.tolist(),
            "final_semi_major_axis_m": elements.semi_major_axis,
            "final_eccentricity": elements.eccentricity,
            "position_error_m": position_error, "velocity_error_m_s": velocity_error,
            "passed": bool(checks["Cumple"].all()),
        },
    }
    return {
        "result": result, "checks": checks, "profile": profile, "sweep": sweep,
        "ratio": ratio,
        "initial_km": np.column_stack((r1 * np.cos(circle), r1 * np.sin(circle))) / 1000,
        "final_km": np.column_stack((r2 * np.cos(circle), r2 * np.sin(circle))) / 1000,
        "arc_km": arc.states[:, :2] / 1000,
        "burns_km": np.array([start[:2], apse[:2]]) / 1000,
    }


# --- Test suite ----------------------------------------------------------

_TEST_LINE = re.compile(r"^(test\w*) \(([\w.]+)\)\s*(?:\n.*?)?\.\.\. (ok|FAIL|ERROR|skipped.*|expected failure|unexpected success)$", re.M)


def run_test_suite(modules: list[str] | None) -> dict[str, Any]:
    """Run unittest in a subprocess and parse its verbose report."""
    command = [sys.executable, "-m", "unittest"]
    command += ["-v", *modules] if modules else ["discover", "-s", "tests", "-t", ".", "-v"]
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src") + os.pathsep + os.environ.get("PYTHONPATH", "")}
    completed = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
    log = completed.stderr + completed.stdout
    rows = []
    for name, dotted, outcome in _TEST_LINE.findall(log):
        parts = dotted.split(".")
        rows.append({
            "Módulo": parts[-3] if len(parts) >= 3 else dotted,
            "Clase": parts[-2] if len(parts) >= 2 else "",
            "Test": name,
            "Resultado": "Correcto" if outcome == "ok" else outcome,
        })
    summary = re.search(r"Ran (\d+) tests? in ([\d.]+)s", log)
    return {
        "table": pd.DataFrame(rows, columns=["Módulo", "Clase", "Test", "Resultado"]),
        "ran": int(summary.group(1)) if summary else len(rows),
        "seconds": float(summary.group(2)) if summary else 0.0,
        "passed": completed.returncode == 0,
        "log": log,
    }


# --- Figures -------------------------------------------------------------

def _earth_sphere() -> go.Surface:
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 20))
    r = EARTH_RADIUS_M / 1000
    return go.Surface(
        x=r * np.cos(u) * np.sin(v), y=r * np.sin(u) * np.sin(v), z=r * np.cos(v),
        colorscale=[[0, "#c9d6e3"], [1, "#c9d6e3"]], showscale=False, opacity=0.35,
        hoverinfo="skip", name="Tierra",
    )


def orbit_3d(orbit_km: np.ndarray, reference_km: np.ndarray, method: str) -> go.Figure:
    """Numerical orbit against the analytic circular reference, in km."""
    fig = go.Figure([
        _earth_sphere(),
        go.Scatter3d(x=reference_km[:, 0], y=reference_km[:, 1], z=reference_km[:, 2], mode="lines",
                     line={"color": NEUTRAL, "width": 2, "dash": "dash"}, name="Solución analítica"),
        go.Scatter3d(x=orbit_km[:, 0], y=orbit_km[:, 1], z=orbit_km[:, 2], mode="lines",
                     line={"color": SERIES[0], "width": 4}, name=f"Propagación {method}",
                     hovertemplate="x %{x:.1f} km<br>y %{y:.1f} km<br>z %{z:.1f} km<extra></extra>"),
        go.Scatter3d(x=[orbit_km[0, 0]], y=[orbit_km[0, 1]], z=[orbit_km[0, 2]], mode="markers",
                     marker={"size": 5, "color": SERIES[1]}, name="Estado inicial"),
    ])
    bound = float(np.max(np.abs(orbit_km))) * 1.1
    axis = {"range": [-bound, bound]}
    fig.update_layout(
        scene={"xaxis": {**axis, "title": "x (km)"}, "yaxis": {**axis, "title": "y (km)"},
               "zaxis": {**axis, "title": "z (km)"}, "aspectmode": "cube"},
        height=520, margin={"l": 0, "r": 0, "t": 10, "b": 0}, legend={"orientation": "h", "y": 1.02},
    )
    return fig


def transfer_plane(data: dict[str, Any]) -> go.Figure:
    """Initial orbit, transfer arc and target orbit in the orbital plane, in km."""
    r = EARTH_RADIUS_M / 1000
    fig = go.Figure()
    fig.add_shape(type="circle", x0=-r, y0=-r, x1=r, y1=r, fillcolor="#c9d6e3", line={"width": 0},
                  opacity=0.6, layer="below")
    for points, color, name, dash in (
        (data["initial_km"], SERIES[0], "Órbita inicial", "solid"),
        (data["final_km"], SERIES[2], "Órbita objetivo", "solid"),
        (data["arc_km"], SERIES[1], "Elipse de transferencia", "dash"),
    ):
        fig.add_trace(go.Scatter(
            x=points[:, 0], y=points[:, 1], mode="lines", name=name,
            line={"color": color, "width": 2, "dash": dash},
            hovertemplate="x %{x:.0f} km<br>y %{y:.0f} km<extra>" + name + "</extra>",
        ))
    burns = data["burns_km"]
    fig.add_trace(go.Scatter(
        x=burns[:, 0], y=burns[:, 1], mode="markers+text", name="Impulsos",
        text=["Δv₁", "Δv₂"], textposition="top center",
        marker={"size": 10, "color": "#0b0b0b", "line": {"width": 2, "color": "#fcfcfb"}},
        hovertemplate="%{text}<extra></extra>",
    ))
    bound = float(np.max(np.abs(data["final_km"]))) * 1.1
    bound = max(bound, float(np.max(np.abs(data["initial_km"]))) * 1.1)
    fig.update_layout(
        xaxis={"range": [-bound, bound], "title": "x (km)", "zeroline": False},
        yaxis={"range": [-bound, bound], "title": "y (km)", "scaleanchor": "x", "zeroline": False},
        height=560, margin={"l": 10, "r": 10, "t": 10, "b": 10}, legend={"orientation": "h", "y": 1.06},
    )
    return fig


def line_chart(frame: pd.DataFrame, x: str, y: str, x_title: str, y_title: str,
               color: str = SERIES[0], log: bool = False, log_x: bool = False) -> alt.Chart:
    """Single-series line with a crosshair tooltip."""
    data = frame[frame[y] > 0] if log else frame
    scale = alt.Scale(type="log") if log else alt.Scale(zero=False)
    x_scale = alt.Scale(type="log") if log_x else alt.Scale(zero=False)
    base = alt.Chart(data).encode(x=alt.X(x, title=x_title, scale=x_scale))
    line = base.mark_line(color=color, strokeWidth=2).encode(y=alt.Y(y, title=y_title, scale=scale))
    hover = alt.selection_point(fields=[x], nearest=True, on="pointerover", empty=False)
    points = base.mark_point(size=60, filled=True, color=color).encode(
        y=alt.Y(y, scale=scale),
        opacity=alt.condition(hover, alt.value(1), alt.value(0)),
        tooltip=[alt.Tooltip(x, title=x_title, format=".2f"), alt.Tooltip(y, title=y_title, format=".3e")],
    ).add_params(hover)
    return (line + points).properties(height=260)


def drift_chart(series: pd.DataFrame) -> alt.Chart:
    """Three dimensionless invariant drifts on one logarithmic axis."""
    names = {"energy": "Energía (relativa)", "angular_momentum": "Momento angular (relativa)",
             "eccentricity": "Excentricidad (absoluta)"}
    long = series.melt(id_vars="time_min", value_vars=list(names), var_name="key", value_name="drift")
    long = long[long["drift"] > 0].assign(Invariante=lambda d: d["key"].map(names))
    hover = alt.selection_point(fields=["time_min"], nearest=True, on="pointerover", empty=False)
    color = alt.Color("Invariante:N", scale=alt.Scale(domain=list(names.values()), range=SERIES),
                      legend=alt.Legend(orient="top", title=None))
    base = alt.Chart(long).encode(
        x=alt.X("time_min:Q", title="Tiempo (min)"),
        y=alt.Y("drift:Q", title="Deriva (adimensional)", scale=alt.Scale(type="log")),
        color=color,
    )
    points = base.mark_point(size=50, filled=True).encode(
        opacity=alt.condition(hover, alt.value(1), alt.value(0)),
        tooltip=[alt.Tooltip("Invariante:N"), alt.Tooltip("time_min:Q", title="Tiempo (min)", format=".2f"),
                 alt.Tooltip("drift:Q", title="Deriva", format=".3e")],
    ).add_params(hover)
    return (base.mark_line(strokeWidth=2) + points).properties(height=300)


def sweep_chart(sweep: pd.DataFrame, ratio: float, total_km_s: float) -> alt.Chart:
    """Hohmann total Δv against radius ratio, marking the current case."""
    line = line_chart(sweep, "ratio", "total_delta_v_km_s", "Relación de radios (mayor/menor, log)",
                      "Δv total (km/s)", log_x=True)
    current = pd.DataFrame({"ratio": [ratio], "total_delta_v_km_s": [total_km_s]})
    marker = alt.Chart(current).mark_point(size=120, filled=True, color=SERIES[1], stroke="#fcfcfb",
                                           strokeWidth=2).encode(
        x="ratio:Q", y="total_delta_v_km_s:Q",
        tooltip=[alt.Tooltip("ratio:Q", title="Relación de radios", format=".3f"),
                 alt.Tooltip("total_delta_v_km_s:Q", title="Δv total (km/s)", format=".4f")],
    )
    return line + marker


def to_json(data: dict[str, Any]) -> str:
    return json.dumps(data, indent=2, allow_nan=False, default=float) + "\n"
