# Fase 1: API y contratos de datos

Interfaces contrastadas con `src/orbital_optimizer` el 2026-10-16, incluida
la ampliación local en `feature/phase-1-kepler` aún sin commit.
Instalación y entorno: [README del proyecto](../README.md).
Ecuaciones y convenciones singulares: [phase1.md](phase1.md).

## Estado y elementos

Todas las funciones siguientes se exportan desde `orbital_optimizer`.
El estado es un array real finito de forma `(6,)`: `[x,y,z,vx,vy,vz]`,
en metros y m/s, en un marco inercial centrado en el cuerpo primario.
No incorpora época ni orientación respecto a un catálogo astronómico.
`mu` es un escalar positivo en m³/s²; por defecto `EARTH_MU`.

```python
OrbitalElements(semi_major_axis, eccentricity, inclination,
                raan=0.0, argument_of_periapsis=0.0, true_anomaly=0.0)
elements_to_state(elements, mu=EARTH_MU)
state_to_elements(state, mu=EARTH_MU)
```

El semieje mayor es positivo, en metros; `0 <= eccentricity < 1` y
`0 <= inclination <= pi`. Los ángulos son radianes; la anomalía es verdadera.
Los otros tres ángulos de entrada pueden superar una vuelta; la conversión
inversa los devuelve en `[0,2*pi)`. Las conversiones usan `OrbitalElements`,
una dataclass congelada. `elements_to_state` devuelve un array nuevo `(6,)`.
`state_to_elements` rechaza estados radiales, no ligados y numéricamente
indistinguibles del límite parabólico. Cerca de singularidades se compara
el estado reconstruido, no cada ángulo por separado.

## Propagación

```python
propagate(initial_state, duration, mu=EARTH_MU, *, method="RK45",
          sample_times=None, max_step=60.0, rtol=1e-10, atol=None)
```

| Entrada | Contrato |
| --- | --- |
| `duration` | Duración positiva finita en segundos; el reloj comienza en cero |
| `method` | Exactamente `"RK45"` o `"RK4"` |
| `sample_times` | `None` o vector no vacío, finito, estrictamente creciente dentro de `[0,duration]`; añade extremos ausentes |
| `max_step` | Techo de paso positivo finito, en segundos |
| `rtol` | Escalar positivo finito, al menos `100*eps_float64` |
| `atol` | Escalar o vector `(6,)` positivo finito; `None` equivale a `[1e-3]*3+[1e-6]*3`, en m y m/s |

La entrada se valida como elipse no radial. RK4 acorta pasos para alcanzar
las muestras y el extremo final; `rtol/atol` se validan pero no controlan su
error. RK45 adapta pasos e interpola las muestras solicitadas. Sin muestreo
explícito se devuelven los extremos de los pasos aceptados. Las tolerancias
locales no garantizan una precisión global determinada.

El resultado `PropagationResult` contiene:

| Campo | Significado |
| --- | --- |
| `times` | `(N,)`, segundos crecientes, incluyendo cero y duración |
| `states` | `(N,6)`, estados SI correspondientes |
| `method`, `mu` | Método utilizado y parámetro gravitacional |
| `nfev` | Evaluaciones de la derivada; no es número de muestras |
| `final_state` | Copia independiente `(6,)` del último estado |

Los arrays de resultado se copian y se marcan como solo lectura. No hay
eventos de colisión ni aplicación de impulsos. Entradas inválidas producen
`ValueError`; los fallos de integración detectados producen `RuntimeError`.
Las conversiones y cálculos orbitales también pueden producir `RuntimeError`
cuando el resultado excede la representación numérica.

## Dinámica, invariantes y costes

| Función | Resultado |
| --- | --- |
| `two_body_dynamics(t, state, mu=EARTH_MU)` | Derivada `(6,)`: m/s y m/s²; dinámica autónoma |
| `specific_energy(state, mu=EARTH_MU)` | Escalar, J/kg |
| `specific_angular_momentum(state)` | Vector `(3,)`, m²/s |
| `eccentricity_vector(state, mu=EARTH_MU)` | Vector `(3,)`, adimensional |
| `node_vector(state)` | Vector `(3,)`, m²/s; nulo en órbita ecuatorial |
| `orbital_period(semi_major_axis, mu=EARTH_MU)` | Escalar en segundos |
| `delta_v(velocity_before, velocity_after)` | Norma de la diferencia de velocidades `(3,)`, m/s |
| `total_delta_v(burn_vectors)` | Suma de normas de vectores `(N,3)`, m/s; `[]` cuesta cero |

`validate_state` y `validate_mu` se encuentran en
`orbital_optimizer.dynamics.two_body`; no se reexportan desde la raíz.
Los cálculos de invariantes no imponen por sí mismos todo el dominio elíptico
que sí comprueba `propagate`. Calcular un coste Δv no modifica una trayectoria.
`src/physics` conserva imports de compatibilidad; para código nuevo se usa
`orbital_optimizer`. No todos los símbolos nuevos se reexportan en `physics`.

## Ejemplo ejecutable

```python
import numpy as np
from orbital_optimizer import (
    OrbitalElements, elements_to_state, orbital_period, propagate,
    specific_energy, total_delta_v,
)

orbit = OrbitalElements(7_000_000.0, 0.0, np.deg2rad(28.5))
initial = elements_to_state(orbit)
period = orbital_period(orbit.semi_major_axis)
result = propagate(initial, period, sample_times=np.linspace(0, period, 101))
print(result.states.shape)  # (101, 6)
print(specific_energy(result.final_state))  # J/kg
print(total_delta_v([[3.0, 4.0, 0.0], [-3.0, -4.0, 0.0]]))  # 10.0 m/s
```

## Experimento y archivos

Desde la raíz del repositorio:

```powershell
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/earth_orbit.json --output outputs/earth_orbit --plots
```

`run_experiment(config)` devuelve `(report, trajectory, reference)`:
un diccionario JSON, un `PropagationResult` y los estados circulares analíticos
`(N,6)`. El experimento admite únicamente órbitas circulares. La API analítica
`propagate_kepler`, documentada abajo, no cambia el alcance de este CLI.

La configuración exige exactamente las claves de
[earth_orbit.json](../configs/earth_orbit.json): `radius_m`, `inclination_rad`,
`mu_m3_s2`, `periods` (entero >=1), `samples` (entero >=2), `method`,
`max_step_s`, `rtol`, `atol` y `acceptance`. Esta última exige los cinco
nombres de métricas del archivo, con límites escalares positivos finitos.
El radio se mide desde el centro, no desde la superficie.

`report.json` tiene `schema_version=1`, configuración efectiva, fecha UTC,
versiones, modelo, periodo, estados inicial/final, tiempo de una propagación,
evaluaciones, `metrics`, `checks`, `validation_passed` y limitaciones.
Cada comprobación acepta `metric <= límite`. Las derivas de energía y vector
de momento angular son relativas al valor inicial; la del vector de
excentricidad es absoluta. Los errores de posición y velocidad son máximos
de normas euclídeas sobre las muestras guardadas.

`trajectory.csv` tiene cabecera
`time_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s`. `--plots` añade `trajectory.png`
y `consistency.png` con Matplotlib opcional; los gráficos convierten posición
a km y tiempo a minutos solo en presentación. `outputs/` se ignora en Git.
El CLI devuelve 0 si acepta, 1 si incumple umbrales y 2 para los errores
capturados de configuración, lectura o ejecución. Conviene usar una carpeta
nueva por ejecución para conservar evidencia anterior.

## Ampliación local: ecuación y propagación de Kepler

```python
solve_kepler(mean_anomaly, eccentricity)
propagate_kepler(initial_state, duration, mu=EARTH_MU, *, sample_times=None)
```

Ambas funciones se exportan desde `orbital_optimizer`. `solve_kepler` acepta
dos escalares reales finitos: anomalía media en radianes y `0 <= e < 1`.
Devuelve la anomalía excéntrica E que satisface `E-e*sin(E)=M`, conservando
las vueltas de M. Reduce internamente M a `[-pi,pi]`. Usa Newton con
salvaguarda por bisección; no expone tolerancia ni límite de iteraciones.
El mal condicionamiento aumenta al acercarse a e=1 y al periapsis.

`propagate_kepler` devuelve `PropagationResult`, con `method="Kepler"` y
`nfev=0` (sin evaluaciones de la derivada). Usa estados SI y el dominio
elíptico no radial de `propagate`. Sin `sample_times` devuelve únicamente
`[0,duration]`; con muestras aplica las mismas reglas de extremos y orden.
No recibe `max_step`, `rtol` ni `atol`.

La evolución usa coeficientes de Lagrange f/g directamente sobre el estado
cartesiano. No construye ángulos orbitales para avanzar el estado; sí utiliza
`state_to_elements` para validar el dominio. Es una solución del modelo de
dos cuerpos evaluada mediante una raíz numérica y aritmética float64, sin
error de truncamiento de integración, pero con errores de raíz y redondeo.
No admite perturbaciones ni burns. Entradas inválidas producen `ValueError`;
la falta de convergencia o estados finales no finitos producen `RuntimeError`.

Implementación: `orbital/kepler.py`, `simulation/kepler_propagator.py`.
Estado de sus pruebas: [registro de validación](phase1-validation.md).

```python
from orbital_optimizer import solve_kepler, propagate_kepler

eccentric_anomaly = solve_kepler(1.0, 0.2)  # 1.1853242038613385 rad
analytic = propagate_kepler(initial, period, sample_times=result.times)
print(analytic.method, analytic.nfev)  # Kepler 0
```

Este fragmento continúa el ejemplo anterior; `initial`, `period` y `result`
proceden de él. No requiere instalar dependencias adicionales.
