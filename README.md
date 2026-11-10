# DAEDALUS — Orbital Trajectory Optimization Engine

[![Tests](https://github.com/AlbertoGuirado/DAEDALUS-Orbital-Trajectory-Optimization-Engine/actions/workflows/tests.yml/badge.svg)](https://github.com/AlbertoGuirado/DAEDALUS-Orbital-Trajectory-Optimization-Engine/actions/workflows/tests.yml)

Autonomous Orbital Transfer Optimizer: optimización matemática de trayectorias orbitales.

**Estado (2026-11-10): Fases 1 y 2 integradas en `develop`. Fase 3 implementada
y validada en `feature/phase-3-optimization`, pendiente de integración.**

- **Fase 1 — Orbital Dynamics Engine.** Motor de dos cuerpos en Python, con
  conversiones orbitales, RK4/RK45, ecuación de Kepler (`solve_kepler`) y
  propagación kepleriana analítica (`propagate_kepler`), validados frente a
  soluciones analíticas.
- **Fase 2 — Baseline analítico.** Transferencia de Hohmann entre órbitas
  circulares coplanarias (`hohmann_transfer`), propelente por Tsiolkovsky
  (`propellant_mass`) y propagación con impulsos instantáneos
  (`propagate_with_impulses`), verificada propagando la trayectoria.
  Transferencia bielíptica (`bielliptic_transfer`). Cambio de plano y Lambert
  no están implementados.
- **Fase 3 — Motor de optimización.** Interfaz común de optimizador con
  semilla, presupuesto e historia; Differential Evolution, enjambre de
  partículas, algoritmo genético, búsqueda aleatoria, Nelder–Mead e híbrido;
  objetivo ponderado adimensional; problemas de transferencia coplanaria de
  dos y de N impulsos verificados por propagación; benchmark multisemilla con
  estadística no paramétrica. Ver [API](docs/phase3-api.md) y
  [validación ejecutada](docs/phase3-validation.md).

Las restricciones de misión, las perturbaciones y el optimizador de misiones
(Fases 4–6) no están implementados en esta rama.

## Ejecutar en PowerShell

Desde `orbital-optimizer`, con Python 3.11 o posterior:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[plots]"
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/earth_orbit.json --output outputs/earth_orbit --plots
```

La dependencia opcional `plots` permite exportar figuras con Matplotlib.
Para usar solo el núcleo, instalar `-e .` y omitir `--plots`. En Linux/macOS
usar `.venv/bin/python`. Los tests usan `unittest`, incluido en Python.

El ejemplo guarda un informe JSON con errores e invariantes, trayectoria CSV y,
con `--plots`, figuras PNG de trayectoria y deriva numérica. El proceso devuelve
un código de salida no nulo si no cumple los criterios configurados.

```python
import numpy as np
from orbital_optimizer import OrbitalElements, elements_to_state, orbital_period, propagate

orbit = OrbitalElements(7_000_000.0, 0.0, np.deg2rad(28.5))
initial = elements_to_state(orbit)
trajectory = propagate(initial, orbital_period(orbit.semi_major_axis))
print(trajectory.final_state)  # m, m/s

from orbital_optimizer import hohmann_transfer, propellant_mass
transfer = hohmann_transfer(6_678_000.0, 42_164_000.0)  # radios en m: LEO -> GEO
print(transfer.total_delta_v, transfer.transfer_time)    # m/s, s
print(propellant_mass(2000.0, transfer.total_delta_v, 320.0))  # kg (m0, Δv, Isp)
```

Las posiciones son radios desde el centro, no altitudes; todas las operaciones
internas usan SI y radianes. Dominio validado: órbitas elípticas no radiales.
No se modelan superficie, colisiones, empuje finito, masa variable, drag ni J2;
los impulsos de la Fase 2 son instantáneos.
Los errores locales del integrador no garantizan por sí solos la precisión
global: el experimento mide invariantes y error contra la solución circular.

Consultar [diseño, ecuaciones, unidades y criterios de validación](docs/phase1.md).
La [validación ejecutada](docs/phase1-validation.md) registra los resultados medidos.
La [API de la Fase 2](docs/phase2-api.md) documenta Hohmann, Tsiolkovsky e impulsos.
`configs/experiments/` contiene configuraciones exploratorias, fuera de la referencia validada.
La implementación está en `src/orbital_optimizer/{dynamics,orbital,simulation,visualization}`.
`src/physics` conserva los imports de la estructura inicial sin duplicar física.
Los archivos vacíos originales de misiones/optimización/agent son reservas de
fases posteriores y no se publican como algoritmos implementados.

El experimento de la Fase 1 es determinista: no precisa seeds. El benchmark
de optimizadores usa semillas explícitas y se lanza con:

```powershell
.\.venv\Scripts\python.exe -m experiments.benchmark_optimizers --config configs/benchmarks/phase3_differential_evolution.json --output outputs/<ejecucion-nueva> --plots --jobs 8
```

Ningún optimizador garantiza un óptimo global; los resultados valen para los
escenarios, presupuestos y semillas configurados.

La interfaz Streamlit local se desarrolla aparte, en la rama `feature/streamlit-ui`.

`requirements-reproducible.txt` registra las versiones del núcleo usadas en la
validación; el informe añade las versiones efectivas y la configuración completa.

## Desarrollo

`main` contiene la versión validada; `develop` integra los siguientes cambios.
Las ramas `feature/*` se crean desde `develop` y se integran mediante pull requests.
La Fase 1 está identificada por `v0.1.0`. Consulta [CONTRIBUTING.md](CONTRIBUTING.md)
para los comandos, las comprobaciones y la organización de las próximas fases.
