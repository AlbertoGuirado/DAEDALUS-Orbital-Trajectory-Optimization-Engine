# Fase 3: API del motor de optimización

Documento parcial, en construcción, de la rama `feature/phase-3-optimization`.
Recoge los incrementos (a) y (b): optimizador genérico, ayudas de Kepler y el
problema de transferencia coplanaria. Quedan para el incremento (c) el
benchmark multisemilla y los gráficos. Unidades SI en todas las magnitudes físicas.

## Capas y dependencias en tiempo de importación (Layering and import-time dependencies)

**Implemented:** el subpaquete `orbital_optimizer.optimization` (`problem.py`,
`result.py`, `optimizer.py`, `differential_evolution.py`, `benchmarks.py`)
importa directamente solo NumPy, la biblioteca estándar y sus propios módulos.
Lo comprueba `tests/test_layering.py`, que analiza el código con `ast` sin
ejecutarlo.

La regla física se aplica a `dynamics/`, `orbital/`, `simulation/`,
`visualization/`, `src/physics/` y `orbital_optimizer/__init__.py`. Esos
archivos no pueden:

- importar `optimization` ni `problems`;
- importar `importlib`, `runpy` ni `pkgutil`;
- usar `exec`, `eval` ni `__import__`.

El adaptador `problems/` queda fuera de ambas reglas, porque puede importar la
física y la optimización.

El paquete raíz no reexporta la API de optimización ni la de `problems`. Ambas
se importan desde `orbital_optimizer.optimization` y `orbital_optimizer.problems`.

**Alcance de la garantía:** la independencia es **estática y de imports
directos**. No es un aislamiento en tiempo de importación. Python ejecuta
primero el paquete padre: `import orbital_optimizer.optimization` ejecuta
`orbital_optimizer/__init__.py`, que importa la API de las fases 1–2.

Comprobación ejecutada el 2026-11-03 en el entorno del proyecto (Python 3.12.1,
NumPy 2.3.3, SciPy 1.16.3), tras importar solo `orbital_optimizer.optimization`:

- se cargaron `orbital_optimizer.dynamics.two_body`, `orbital.elements`,
  `orbital.kepler`, `orbital.mechanics`, `orbital.transfers`,
  `simulation.impulsive`, `simulation.kepler_propagator` y
  `simulation.propagator`;
- se cargaron `scipy` y `scipy.integrate`, este último a través de `simulation.propagator`;
- no se cargó `matplotlib`.

**Consecuencias:**

- SciPy debe estar instalado, y su coste de importación se paga aunque solo
  se use la optimización genérica.
- No hay dependencia circular, porque ninguna capa física importa este subpaquete.

**Planned (no implementado):** carga perezosa en el paquete raíz, por ejemplo
con un `__getattr__` de módulo (PEP 562), para que importar el optimizador no
cargue la física ni SciPy. Exigiría cambiar `orbital_optimizer/__init__.py`,
que forma parte de la API pública de las fases 1–2, y por tanto requiere
aprobación previa.

## Optimizador genérico (`orbital_optimizer.optimization`)

**Implemented:**

- **Contratos del problema:** `Bounds`, `EvaluationStatus`, `Evaluation`, el
  protocolo `Problem` y `FunctionProblem`.
- **Terminación y resultados:** `TerminationCriteria`, `TerminationReason`,
  `OptimizationHistory` y `OptimizationResult`, con `to_dict()` en JSON estricto.
- **Optimizador:** `EvaluationBudget` y `DifferentialEvolution`, variante
  `"DE/rand/1/bin@1"`, con NP = 10·D, F = 0,5 y CR = 0,9.

Las ecuaciones, el orden de las extracciones aleatorias y la precedencia de
terminación están en los docstrings de `differential_evolution.py` y `result.py`.

Contrato de aislamiento de `EvaluationBudget`:

- el problema recibe una copia de solo lectura de `x`;
- `evaluation.x` debe coincidir exactamente con el punto evaluado; si no, se
  lanza `ValueError`;
- la mejor evaluación se guarda reconstruida a partir de la copia privada de
  `x`, con los metadatos copiados mediante `copy.deepcopy`.

Por tanto, modificar después los objetos del problema no altera el resultado.

**Referencia fijada (*pin*):** `tests/test_differential_evolution.py` fija
`f_best` y un hash SHA-256 de una ejecución de referencia: Rosenbrock-2D,
600 evaluaciones y semilla 2026.

- Es un **detector de deriva**, no un oráculo de corrección.
- Solo es válido con **NumPy 2.3.x**; con otra versión, la prueba se omite
  (`skipTest`) indicando el motivo.
- El trabajo `pinned-runtime` de `.github/workflows/tests.yml` instala
  `requirements-reproducible.txt` (NumPy 2.3.3, SciPy 1.16.3) y el proyecto con
  `--no-deps`, en Ubuntu y Windows con Python 3.12, para que la referencia se
  ejecute en CI y no se omita.
- A 2026-11-03 ese trabajo de CI aún no se ha ejecutado. La referencia solo se
  ha comprobado en Windows 11.

## Ayudas de Kepler (aditivas, `orbital_optimizer.orbital.kepler`)

**Implemented**, sin cambiar las funciones existentes ni reexportarlas desde la raíz:

| Función | Resultado |
|---|---|
| `true_to_eccentric_anomaly(nu, e)` | E = 2·atan2(√(1−e)·sin(ν/2), √(1+e)·cos(ν/2)), en rad; conserva las vueltas de ν |
| `true_to_mean_anomaly(nu, e)` | M = E − e·sin E, en rad; conserva las vueltas de ν |
| `elliptic_time_of_flight(nu_start, nu_end, a, e, mu=EARTH_MU)` | ((M(ν_end) − M(ν_start)) mod 2π)/n, en s, con n = √(μ/a³); primera llegada, en [0, T) |

Con 0 ≤ e < 1, a > 0 en m y μ en m³/s². Las entradas inválidas producen
`ValueError`. En `tests/test_kepler_anomalies.py` se comprueban:

- el caso circular;
- las identidades geométricas de E;
- la ida y vuelta con `solve_kepler`;
- la semielipse de Hohmann frente a π·√(a³/μ);
- la llegada al ν pedido mediante `propagate_kepler`.

## Problema de transferencia coplanaria (`orbital_optimizer.problems`)

**Implemented:** `CoplanarTwoImpulseProblem(r1, r2, mu=EARTH_MU)`, solo para
transferencias de subida (0 < r1 < r2).

- Si r1 == r2 o r2 < r1, lanza `ValueError`.
- Si la referencia de Hohmann o los límites no son finitos, lanza `RuntimeError`.

**Variables de decisión:**

| Variable | Significado | Unidad | Dominio |
|---|---|---|---|
| Δv1 | Módulo del primer impulso | m/s | [0, (√2−1)·v_c1] |
| α1 | Ángulo respecto a la velocidad local; positivo hacia fuera en la dirección radial | rad | [−π, π), periódica |

Dentro de esa caja, todo candidato es ligado y prógrado. El segundo impulso
circulariza en el primer cruce hacia fuera de r2, y τ se deriva de él.

**Objetivo:** J = |Δv1| + |Δv2|, en m/s.

**Cadena de evaluación:**

1. Predicados analíticos, antes de llamar a la física.
2. `state_to_elements`.
3. ν2 = arccos((p/r2 − 1)/e).
4. `elliptic_time_of_flight`.
5. `propagate_kepler` hasta τ.
6. Δv2 y J.

No se usa RK45 dentro del objetivo.

**Metadatos, con unidades en el nombre:**

- `delta_v_total_m_s`, `delta_v1_m_s`, `delta_v2_m_s`;
- `first_burn_m_s`, `second_burn_m_s`;
- `transfer_time_s`;
- `semi_major_axis_m`, `eccentricity`, `periapsis_radius_m`, `apoapsis_radius_m`;
- `departure_true_anomaly_rad`, `arrival_true_anomaly_rad`;
- `arrival_radius_residual_m`, `tangency_clip`, `reach_margin_relative`;
- `reference_delta_v_m_s`, `reference_transfer_time_s`.

Los contadores `captured_value_errors` y `captured_runtime_errors` exponen
las excepciones de la física capturadas. En las pruebas y en las ejecuciones
de DE han sido cero.

**Tolerancia de tangencia.** El óptimo de Hohmann está sobre la frontera
r_a = r2. El criterio de alcance es

Q = 2ε·r2² + 2μ·r2 − h² ≥ −K·u·S,  con S = 2r2²(v²/2 + μ/r1) + 2μ·r2 + h²,
u = 2⁻⁵³ y K = 32.

Derivación:

- unos 10 redondeos en el camino más largo dan ≈ 10u·S;
- el error de representación del Δv1 de Hohmann calculado (≈ 4,5u·v)
  contribuye ≤ 9u·S;
- la suma es ≈ 19u·S, y K es la siguiente potencia de dos.

K no depende de r2/r1 ni de la escala. Se fijó con este argumento, no
ajustándolo a las pruebas.

Barrido de confirmación (2026-11-03; 3 valores de μ × 3 escalas × 40
cocientes r2/r1 en [1,0001; 100]):

- en el Δv1 de Hohmann, |Q|/(u·S) ≤ 1,56;
- Q resultó negativo en 153 de 360 casos, que un test estricto habría rechazado.

Banda infactible admitida, medida por bisección con α1 = 0:

- **LEO→GEO:** 93 ulps de Δv1 (4,2e-11 m/s), equivalentes a un déficit de
  apoapsis de 6,3e-14 relativo (2,7e-6 m). En el borde de la banda,
  J = J_H − 1,4e-10 m/s.
- **Radios casi iguales:** la banda en Δv1 crece (relativa 2,1e-6 para
  r2/r1 = 1,0001), porque Q apenas depende de Δv1 cuando r2 ≈ r1. Un umbral
  "J ≥ J_H − ε" debe tenerlo en cuenta en ese régimen.

**Verificación (RV-3.04):** `verify_two_impulse_plan(plan)` hace lo siguiente:

1. propaga el plan con `propagate_with_impulses` (RK45 por defecto);
2. comprueba explícitamente la finitud y llama a `state_to_elements` sobre el
   estado final;
3. exige |a − r2| ≤ 1 m, e < 1e-6 y |‖r‖ − r2| ≤ 1 m, los criterios de llegada
   de la Fase 2;
4. propaga una costa de un periodo de la órbita objetivo con `propagate` y
   exige las mismas tolerancias de radio y excentricidad.

`hohmann_two_impulse_plan(r1, r2)` construye el plan de Hohmann orientado en
subida o en bajada. En la Fase 3, la bajada se usa solo para verificación.

**Evidencia de las pruebas** (`tests/test_coplanar_transfer_problem.py`,
2026-11-03; Windows 11, Python 3.12.1, NumPy 2.3.3):

- **Decisión de Hohmann:** |J − J_H| ≤ 1e-6 m/s y τ = t_H con `places=6`.
- **Oráculo cerrado independiente:** coincide con el evaluador en 152 puntos
  factibles de una malla de 59 × 60. El cociente error observado / tolerancia
  derivada fue como máximo 0,021 para J y 0,011 para τ (146 puntos bien condicionados).
- **Muestreo factible:** ninguna de 1000 muestras factibles mejora J_H − 1e-6 m/s.
- **Verificación RK45:**
  - el plan de Hohmann pasa, con un error de semieje de 0,064 m;
  - la bajada GEO→LEO pasa;
  - un plan sin segundo impulso falla (e = 0,727).
- **DE de extremo a extremo, 3 semillas × 2000 evaluaciones:** gaps de
  4,1e-4, 4,4e-5 y 8,1e-3 m/s; los tres pasan la verificación.

Estas tres ejecuciones son una guarda de regresión, no el benchmark
multisemilla, que corresponde al incremento (c).

**Known Limitations:**

- **Clase de transferencias:** solo dos impulsos, coplanaria, entre órbitas
  circulares, con tiempo libre; sin *phasing*, superficie, masa ni perturbaciones.
- **Trayectorias no físicas:** se aceptan candidatos con r_p < R⊕ (trayectorias
  subterráneas), porque la altitud mínima es una restricción de la Fase 4.
- **Advertencia heredada de la Fase 2:** `hohmann_transfer` emite una
  `RuntimeWarning` y devuelve NaN a escalas extremas (inf − inf fuera de su
  bloque `errstate`). El adaptador la suprime solo en sus propias llamadas y
  rechaza el NaN; la función de la Fase 2 no se ha modificado.
