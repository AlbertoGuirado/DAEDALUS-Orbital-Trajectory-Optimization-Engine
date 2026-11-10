# Fase 3: API del motor de optimización

Estado a **2026-11-10**, rama `feature/phase-3-optimization`: incrementos
(a) motor DE e interfaces, (b) ayudas de Kepler y problema coplanario de dos
impulsos, (c) benchmark multisemilla y (d) optimizadores de comparación,
objetivo ponderado y problema de N impulsos. Unidades SI en todas las
magnitudes físicas; ángulos en radianes.

Resultados medidos y trazabilidad de requisitos: [phase3-validation.md](phase3-validation.md).

## Capas y dependencias en tiempo de importación (Layering and import-time dependencies)

**Implemented:** el subpaquete `orbital_optimizer.optimization` importa
directamente solo NumPy, la biblioteca estándar y sus propios módulos. Lo
comprueba `tests/test_layering.py`, que analiza el código con `ast` sin
ejecutarlo.

La regla física se aplica a `dynamics/`, `orbital/`, `simulation/`,
`visualization/`, `src/physics/`, `orbital_optimizer/__init__.py` y
`orbital_optimizer/constants.py`. Esos archivos no pueden:

- importar `optimization`, `problems`, `benchmarking` ni `mission`;
- importar `importlib`, `runpy` ni `pkgutil`;
- usar `exec`, `eval` ni `__import__`.

Por encima de la física las capas se ordenan así, y cada una puede importar
las anteriores pero nunca las posteriores:

| Capa | Contenido | Puede importar |
|---|---|---|
| `optimization/` | Contratos, optimizadores, objetivo ponderado | NumPy y biblioteca estándar |
| `problems/` | Adaptadores de problemas orbitales | Física y `optimization` |
| `benchmarking/` | Ejecución multisemilla, estadística, informe y figuras | Lo anterior y `scipy.stats` |
| `mission/` | Reservada para la Fase 6B | Todo lo anterior |

Ningún archivo de `src/` importa `streamlit`, `plotly`, `pandas` ni `altair`:
la Verification Console depende del núcleo, nunca al revés.
`constants.py` no importa nada (RD-01).

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
que forma parte de la API pública de las fases 1–2.

## Optimizador genérico (`orbital_optimizer.optimization`)

**Implemented:**

- **Contratos del problema:** `Bounds`, `EvaluationStatus`, `Evaluation`, el
  protocolo `Problem` y `FunctionProblem`.
- **Terminación y resultados:** `TerminationCriteria`, `TerminationReason`,
  `OptimizationHistory` y `OptimizationResult`, con `to_dict()` en JSON estricto
  (`jsonable` convierte los no finitos en `"inf"`, `"-inf"` y `"nan"`).
- **Infraestructura común:** `EvaluationBudget` (presupuesto duro, aislamiento
  y clave de comparación), `check_termination`, `GenerationLog` y `finish_run`.
- **Optimizadores:** los seis de la tabla siguiente. Todos implementan el
  protocolo `Optimizer`: `minimize(problem, termination, *, seed)`.

Todas las variables se buscan en el hipercubo unidad mediante el mapa afín de
`Bounds` (RF-3.08). Las variables periódicas (ángulos) ocupan el intervalo
semiabierto `[lower, upper)`.

### Optimizadores

| Clase | Identificador | Tipo | Valores por defecto y su origen |
|---|---|---|---|
| `DifferentialEvolution` | `DE/rand/1/bin@1` | Global, poblacional | NP = 10·D, F = 0,5, CR = 0,9 (Storn y Price, 1997) |
| `ParticleSwarm` | `PSO/gbest/inertia@1` | Global, poblacional | NP = 10·D, w = 0,7298, c1 = c2 = 1,49618 (Clerc y Kennedy, 2002), v_max = 0,5 |
| `GeneticAlgorithm` | `GA/tournament2/SBX/polynomial/elitist@1` | Global, poblacional | NP = 10·D, p_c = 0,9, η_c = 15, p_m = 1/D, η_m = 20 (Deb et al., 2002) |
| `RandomSearch` | `RandomSearch/uniform@1` | Línea base nula | Lotes de 10·D muestras uniformes independientes |
| `NelderMead` | `NelderMead/Lagarias1998@1` | Local, sin derivadas | Reflexión 1, expansión 2, contracción ½, reducción ½; arista inicial 0,05 |
| `HybridOptimizer` | `Hybrid(<global> -> NelderMead/Lagarias1998@1)@1` | Global + refinamiento final | 80 % del presupuesto a la fase global |

Ningún valor por defecto se ajustó sobre los problemas de DAEDALUS. Ningún
optimizador garantiza un óptimo global.

**Contrato común** (comprobado para todos en `tests/test_optimizers.py` y
`tests/test_differential_evolution.py`):

- **Semilla (RF-3.16, RR-01):** un único `Generator(PCG64(seed))` por
  ejecución; nunca el estado global. Misma semilla, configuración, problema y
  entorno ⇒ resultado idéntico bit a bit, salvo `wall_time_s`.
- **Presupuesto (RF-3.17):** `max_evaluations` es un límite duro que cuenta
  también las evaluaciones no válidas. Nunca se supera.
- **Historia (RF-3.18):** mejor valor válido tras cada evaluación y un
  registro por generación.
- **Comparación de candidatos:** clave lexicográfica `(0, f)` si la evaluación
  es `VALID` y `(1, 0)` si no; gana la menor. Un punto no válido nunca
  sustituye a uno válido; dos no válidos empatan.
- **Terminación:** el primer criterio cumplido según `TERMINATION_PRECEDENCE`.

El orden de las extracciones aleatorias de cada algoritmo es parte de su
versión (`@1`) y está en el docstring de su módulo. Los tests lo reproducen
con una reimplementación escalar independiente para DE, PSO, GA y búsqueda
aleatoria.

**Operadores con variables periódicas:**

| Algoritmo | Variable periódica | Variable acotada |
|---|---|---|
| DE | Envuelve módulo 1; la mutación **no** usa el arco más corto | Punto medio entre el progenitor y el límite violado |
| PSO | Diferencias por el arco más corto; envuelve módulo 1 | Recorta al límite y anula esa componente de la velocidad |
| GA | Cruza con la imagen más cercana del segundo progenitor; envuelve | Recorta a [0, 1] |
| Nelder–Mead | Símplex sin envolver; módulo 1 solo al evaluar | Refleja en la cara violada |

**Refinamiento local y esquema híbrido (RF-3.13, RF-3.14):**

- Nelder–Mead refleja en los límites en lugar de recortar. Con recorte el
  símplex colapsa sobre una cara: el 2026-11-10 convergió a un punto no
  estacionario en una cuadrática convexa con mínimo interior.
- El híbrido usa acoplamiento de refinamiento final ([03-fundamentos](03-fundamentos-cientificos-y-algoritmicos.md) § 12.2).
  La fase global recibe ⌊0,8·B⌋ evaluaciones y la local exactamente las
  restantes, de modo que el total no supera B.
- La arista inicial del símplex es la dispersión final de la población
  global, limitada a [1e-6; 0,05].
- `x_tol`, `f_tol`, estancamiento y `max_generations` se aplican solo a la
  fase global; `target_objective`, a ambas.

### Aislamiento de `EvaluationBudget`

- el problema recibe una copia de solo lectura de `x`;
- `evaluation.x` debe coincidir exactamente con el punto evaluado; si no, se
  lanza `ValueError`;
- la mejor evaluación se guarda reconstruida a partir de la copia privada de
  `x`, con los metadatos copiados mediante `copy.deepcopy`.

Por tanto, modificar después los objetos del problema no altera el resultado.

### Referencia fijada (*pin*) de DE

`tests/test_differential_evolution.py` fija `f_best` y un hash SHA-256 de una
ejecución de referencia: Rosenbrock-2D, 600 evaluaciones y semilla 2026.

- Es un **detector de deriva**, no un oráculo de corrección.
- Solo es válido con **NumPy 2.3.x**; con otra versión, la prueba se omite
  (`skipTest`) indicando el motivo.
- El trabajo `pinned-runtime` de `.github/workflows/tests.yml` instala
  `requirements-reproducible.txt` (NumPy 2.3.3, SciPy 1.16.3) y el proyecto con
  `--no-deps`, en Ubuntu y Windows con Python 3.12, para que la referencia se
  ejecute en CI y no se omita.
- A 2026-11-10 ese trabajo de CI aún no se ha ejecutado. La referencia solo se
  ha comprobado en Windows 11.
- El hash sigue siendo el del 2026-11-03: los incrementos (c) y (d) no
  cambiaron el registro de una ejecución de DE.

## Objetivo ponderado (`optimization.objective`, `problems.objectives`)

**Implemented (RF-3.19 a RF-3.21):**

    J(x) = Σ_i w_i · m_i(x) / m_ref,i

Cada término divide una métrica entre su escala de referencia en la misma
unidad, así que todos los sumandos son adimensionales y ningún peso esconde
una conversión de unidades.

- `ObjectiveTerm(name, metric, weight, reference, unit, rationale)`: el peso
  es adimensional y no negativo; la referencia, positiva. `unit` y `rationale`
  son obligatorios: un término sin ellos se rechaza (RF-3.20).
- `WeightedObjective(terms)`: exige nombres únicos y algún peso positivo. Los
  términos con peso cero no se leen.
- `ScalarizedProblem(problem, objective)`: envuelve cualquier `Problem` cuyas
  evaluaciones válidas publiquen las métricas en `metadata`. Conserva el
  objetivo original en `base_objective` y la contribución de cada término en
  `objective_terms`.

`transfer_objective(reference, *, delta_v, time, energy, risk, rationale)`
declara los cuatro términos para transferencias, referidos a la transferencia
de Hohmann entre las mismas órbitas, de modo que cada uno vale 1 en ella:

| Término | Métrica | Unidad | Referencia |
|---|---|---|---|
| `delta_v` | Σ‖Δv_k‖ | m/s | Δv total de Hohmann |
| `time` | Instante del último impulso | s | Tiempo de Hohmann |
| `energy` | E = Σ‖Δv_k‖²/2 | J/kg | (Δv₁² + Δv₂²)/2 de Hohmann |
| `risk` | R = min(r1, r2) / min_t r(t) | adimensional | 1 |

Energía y riesgo son **indicadores declarados, no modelos validados**:

- E es la energía cinética por unidad de masa que cada impulso imparte en el
  marco que acompaña al vehículo. Para un Δv total fijo es mínima con impulsos
  iguales. No es la energía química del propelente.
- R compara el máximo acercamiento al cuerpo central con la órbita más baja
  de la misión: vale 1 si la transferencia nunca desciende por debajo de ella.
  No es una probabilidad.

Los pesos por defecto son `delta_v = 1` y cero el resto. El llamante debe
justificarlos en `rationale`. Una suma ponderada solo alcanza la parte convexa
de un frente de Pareto ([03-fundamentos](03-fundamentos-cientificos-y-algoritmicos.md) § 14.3).

## Ayudas de Kepler (aditivas, `orbital_optimizer.orbital.kepler`)

**Implemented**, sin cambiar las funciones existentes ni reexportarlas desde la raíz:

| Función | Resultado |
|---|---|
| `true_to_eccentric_anomaly(nu, e)` | E = 2·atan2(√(1−e)·sin(ν/2), √(1+e)·cos(ν/2)), en rad; conserva las vueltas de ν |
| `true_to_mean_anomaly(nu, e)` | M = E − e·sin E, en rad; conserva las vueltas de ν |
| `elliptic_time_of_flight(nu_start, nu_end, a, e, mu=EARTH_MU)` | ((M(ν_end) − M(ν_start)) mod 2π)/n, en s, con n = √(μ/a³); primera llegada, en [0, T) |
| `radius_extrema_on_arc(nu_start, duration, a, e, mu=EARTH_MU)` | (r_min, r_max) del arco recorrido durante `duration`, en m |

Con 0 ≤ e < 1, a > 0 en m y μ en m³/s². Las entradas inválidas producen
`ValueError`.

`radius_extrema_on_arc` es exacta para el movimiento de dos cuerpos: r(ν)
solo cambia de sentido en los ápsides, así que los extremos están en los
puntos finales del arco o en un ápside interior. Con M₀ = M(ν_start) reducida
a [0, 2π) y barrido s = n·`duration`:

- el periapsis está dentro si M₀ + s ≥ 2π;
- el apoapsis está dentro si M₀ ≤ π ≤ M₀ + s, o si M₀ + s ≥ 3π;
- si s ≥ 2π, ambos.

No muestrea el arco, de modo que un periapsis situado entre dos muestras no
puede pasar inadvertido (base de RF-4.05).

Comprobaciones (`tests/test_kepler_anomalies.py`, `tests/test_arc_extrema_and_bielliptic.py`):

- el caso circular y las identidades geométricas de E;
- la ida y vuelta con `solve_kepler`;
- la semielipse de Hohmann frente a π·√(a³/μ);
- la llegada al ν pedido mediante `propagate_kepler`;
- los extremos frente a la ecuación de la cónica en anomalías elegidas a mano
  y frente a un muestreo denso de 2001 puntos en 40 arcos aleatorios.

La transferencia bielíptica (`bielliptic_transfer`, RF-2.05) está documentada
en [phase2-api.md](phase2-api.md#transferencia-bielíptica).

## Problema de transferencia coplanaria (`orbital_optimizer.problems`)

**Implemented:** `CoplanarTwoImpulseProblem(r1, r2, mu=EARTH_MU)`, solo para
transferencias de subida (0 < r1 < r2). Es el problema de referencia validado
de la Fase 3.

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
subida o en bajada.

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

Estas tres ejecuciones son una guarda de regresión. El benchmark multisemilla
está en [phase3-validation.md](phase3-validation.md).

## Problema coplanario generalizado de N impulsos

**Implemented (RF-3.05, RF-3.07):** `CoplanarMultiImpulseProblem(r1, r2, mu=EARTH_MU, *, impulses=2)`.
Admite subida o bajada (r1 ≠ r2) y N ≥ 2 impulsos. Para N = 2 y r2 > r1
define el mismo objetivo que el problema de referencia.

**Marco local de un impulso:** r̂ (radial, hacia fuera) y θ̂ = ẑ × r̂
(transversal, prógrado). Un impulso libre k es

    Δv_k = Δv_k·(sin α_k · r̂ + cos α_k · θ̂)

α = 0 es prógrado, α = ±π retrógrado y 0 < α < π apunta hacia fuera.

**Variables de decisión**, D = 2 + 3·(N − 2):

    x = (Δv_1, α_1, [M_k, Δv_k, α_k  para k = 2 … N−1])

| Variable | Significado | Unidad | Dominio |
|---|---|---|---|
| Δv_k | Módulo del impulso libre k | m/s | [0, (√2−1)·√(μ/min(r1, r2))] |
| α_k | Dirección del impulso en el marco local | rad | [−π, π) en subida; [0, 2π) en bajada; periódica |
| M_k | Barrido de anomalía media de la costa entre los impulsos k−1 y k | rad | [0, 2π), periódica |

**Instantes de los impulsos.** Se derivan de la decisión: t_1 = 0 y
t_k = t_{k−1} + M_k/n, con n = √(μ/a³) de la órbita que se recorre. El
primero es cero por la simetría de rotación de la órbita circular inicial.

**Último impulso.** No es variable de decisión: circulariza en el primer
cruce de r2 posterior al impulso N−1 (hacia fuera si el vehículo está por
debajo de r2, hacia dentro si está por encima). Así se elimina por
construcción la restricción de igualdad "llegar a la órbita objetivo".

**Corte del ángulo.** Los dos intervalos de α describen las mismas
direcciones. El corte se sitúa en el sentido opuesto al del impulso de las
soluciones de tipo Hohmann, para que el óptimo esperado no caiga en la
costura de la variable periódica. Motivo: la mutación de `DE/rand/1/bin@1`
envuelve las variables periódicas pero no usa el arco más corto. Con el corte
en el óptimo (α = ±π en una bajada), sus gaps finales en cinco semillas con
2000 evaluaciones fueron de 0,04 a 12,9 m/s (GEO→LEO, 2026-11-10).

**Predicados de dominio**, tras cada impulso libre: órbita ligada y no
numéricamente parabólica; ni radial ni retrógrada. Después, el criterio de
alcance de r2. Una costa nula entre dos impulsos (impulsos coincidentes) es
`INVALID_DOMAIN`, porque los instantes deben ser estrictamente crecientes (DT-13).

**Tolerancia de tangencia tras arcos de costa.** Un impulso libre puede
seguir a un arco propagado, y la transferencia bielíptica es tangente al
objetivo tras uno. El estado propagado arrastra el redondeo de la propagación
analítica, cuyo modelo de error acota el semieje recalculado en
4·eps·(3+e)/(1−e) = 8u·κ, con κ = (3+e)/(1−e) del arco:

- término de energía: |2r2²·δε| ≤ 8u·κ·S;
- término de momento angular: |2h·δh| ≤ 8u·κ·S.

Cada arco de costa previo añade 16·κ al factor:

    aceptar si Q ≥ −(32 + 16·Σ_arcos κ_arco)·u·S

La constante 16 sale de este argumento, no de un ajuste. Confirmación
(2026-11-10; decisiones bielípticas con r2/r1 de 0,05 a 20 y apoapsis de
hasta 12 000·r1, κ hasta 2,4e4):

- |Q|/(u·S·κ) fue como máximo 0,61;
- Q fue negativo en 4 de 7 casos, que la tolerancia sin el término de arco rechazaba.

La banda infactible admitida crece en consecuencia: un déficit relativo de
unos 16·κ·u en el radio del ápside (2e-14 para un arco con e = 0,7; 4e-11
para e = 0,9998). El metadato `tangency_factor` registra el factor usado.

**Condicionamiento en la tangencia.** La anomalía de cruce es arccos(c) con
c = ±1 en la tangencia, donde un error d en c la desplaza √(2d). El tiempo de
una llegada tangente arrastra por ello un error de hasta unos 4·√(F·u) rad de
anomalía verdadera (F es el factor de tangencia): del orden de 1e-3 s en las
transferencias probadas, frente a 1e-6 s en un Hohmann de un solo arco. El
objetivo no se ve afectado a primer orden.

**Metadatos:**

- `delta_v_total_m_s`, `burn_delta_v_m_s` (N), `burn_times_s` (N), `burn_vectors_m_s` (N, 3);
- `transfer_time_s`: instante del último impulso;
- `minimum_radius_m`, `maximum_radius_m`: extremos exactos de ‖r‖ en toda la
  transferencia, obtenidos con `radius_extrema_on_arc` en cada arco;
- `impulse_kinetic_energy_j_kg`, `closest_approach_index`;
- `final_arc_semi_major_axis_m`, `final_arc_eccentricity`;
- `arrival_radius_residual_m`, `tangency_clip`, `reach_margin_relative`, `tangency_factor`;
- `reference_delta_v_m_s`, `reference_transfer_time_s`, `impulses`.

**Decisiones de referencia:** `hohmann_decision` (N = 2; α = 0 en subida y π
en bajada) y `bielliptic_decision(rb)` (N = 3).

**Planes y verificación (RV-3.04):** `ImpulsivePlan(initial_state, burn_times, burn_vectors, final_radius, mu)`
y `verify_impulsive_plan(plan, *, radius_tolerance_m=1.0, eccentricity_tolerance=1e-6, method, max_step, rtol, atol)`.
Los criterios son los de `verify_two_impulse_plan`. Las tolerancias y las
opciones del integrador son explícitas y quedan registradas. Para
transferencias mucho más largas que LEO–GEO se debe **endurecer el
integrador**, no relajar la tolerancia de llegada: una bielíptica de 13 días
(r2 = 20·r1, rb = 60·r1) falla con los valores por defecto (error de radio
2,1 m) y pasa con `rtol=1e-12` y `atol=[1e-6]*3+[1e-9]*3` (0,005 m). La
medición y su reproducción están en [phase3-validation.md](phase3-validation.md).

## Propagación impulsiva con arcos (`simulation.impulsive`)

**Implemented (DT-11):** `propagate_impulsive_arcs(initial_state, impulses, duration, mu, *, method, samples_per_arc, max_step, rtol, atol)`
devuelve un `ImpulsiveTrajectory` con un `PropagationResult` por arco de vuelo libre.

- Mismo contrato de entrada que `propagate_with_impulses`. Con RK45 y RK4 el
  estado final es idéntico al de esa función.
- `method="Kepler"` usa la propagación analítica.
- `samples_per_arc` pide muestras equiespaciadas en cada arco.
- En cada unión la posición es continua y la velocidad salta exactamente el
  vector del impulso.
- `sampled_radius_extrema()` devuelve extremos **muestreados**, no cotas en
  tiempo continuo: con 7 muestras, un arco que pasa por el periapsis lo
  sobrestimó en más de 100 km en el test. Para una cota exacta se usa
  `radius_extrema_on_arc`.

`propagate_with_impulses` no cambia.

## Benchmark multisemilla (`orbital_optimizer.benchmarking`)

**Implemented (incremento c; base de la Fase 6A):**

    python -m experiments.benchmark_optimizers --config configs/benchmarks/<nombre>.json \
        --output outputs/<carpeta-nueva> [--plots] [--jobs N]

**Configuración** (JSON versionado, esquema `daedalus.benchmark/1`): semillas
(`first`, `count`), algoritmos (`id`, `type`, `settings`), escenarios y
criterios de aceptación. Las claves o tipos desconocidos se rechazan, para
que una errata no caiga en un valor por defecto.

| Tipo de escenario | Parámetros | Referencia |
|---|---|---|
| `function` | `function` (sphere, rosenbrock, rastrigin), `dimension`, `lower`, `upper`, `success_threshold` | Mínimo conocido |
| `coplanar_transfer` | `initial_radius_m`, `final_radius_m`, `impulses`, `success_threshold_m_s`, `lower_threshold_m_s`, `verification` | Hohmann, o bielíptica con `reference` |

**Reglas que impone el ejecutor:**

- todos los algoritmos reciben el presupuesto exacto del escenario, sin otro
  criterio de parada (RF-6A.03);
- las mismas semillas para cada pareja escenario–algoritmo (RF-6A.02);
- el tiempo de pared se registra y nunca detiene una ejecución;
- cada mejor plan de transferencia se verifica por propagación numérica (RV-3.04);
- la primera semilla de cada pareja se repite y se compara bit a bit (RV-3.03).

`--jobs` reparte las ejecuciones entre procesos. Los registros no dependen de
`--jobs`, salvo el tiempo de pared.

**Criterios por ejecución:**

- **éxito:** punto válido con f − f_ref ≤ umbral y, si hay verificador, plan verificado;
- **anomalía:** f − f_ref < −umbral inferior, es decir, mejor que un óptimo
  demostrado (RV-3.01).

**Estadística** (`benchmarking.statistics`):

| Función | Uso |
|---|---|
| `quantile`, `summarize` | Estadísticos de orden por CDF empírica invertida: todo cuantil es un valor observado y +inf no requiere aritmética |
| `wilson_interval` | Intervalo de Wilson al 95 % de la tasa de éxito |
| `evaluations_to_target`, `convergence_band` | Evaluaciones hasta el umbral; cuantiles entre semillas del mejor valor |
| `mann_whitney`, `vargha_delaney_a12` | Comparación de dos muestras independientes y tamaño del efecto |
| `holm_correction` | Corrección de comparaciones múltiples de Holm |
| `friedman`, `wilcoxon_signed_rank` | Comparación entre escenarios con tres o más, o dos, algoritmos |

Se usan mediana y cuartiles en lugar de media y desviación típica porque las
distribuciones son asimétricas y pueden contener +inf (ejecuciones sin punto
válido). Los tests de significación proceden de `scipy.stats`; la corrección
de Holm y A12 están implementadas aquí y se comprueban a mano.

**Aceptación** (clave `acceptance`, opcional): `max_anomalies`,
`require_verification` y `minimum_success_rate` por escenario y algoritmo. El
proceso devuelve 0 si se cumplen y la reproducibilidad es total, 1 si no, y 2
ante un error de configuración o de ejecución.

**Salidas** (SI):

| Archivo | Contenido |
|---|---|
| `report.json` | Configuración, entorno, versión del código, resumen estadístico, reproducibilidad y aceptación |
| `runs.csv` | Una fila por ejecución |
| `convergence.csv` | Cuartiles entre semillas del mejor valor |
| `*.png` (con `--plots`) | Convergencia, distribución final, tasa de éxito, coste frente a tiempo de pared y mejor transferencia verificada |

Todas las figuras llevan ejes con unidades y declaran el estadístico agregado
(RF-6A.08). Sus valores están en los CSV.

**Evidencia versionable:** `outputs/` no se versiona (RR-07).
`python -m experiments.export_evidence --report <report.json> --output docs/evidence/<nombre>.json`
deriva del informe un registro compacto con el SHA-256 del informe de origen.
Nada se escribe a mano.

## Known Limitations

- **Clase de transferencias:** coplanaria, entre órbitas circulares, con
  tiempo libre; sin *phasing*, cambio de plano, superficie, masa ni perturbaciones.
- **Trayectorias no físicas:** se aceptan candidatos con r_p < R⊕
  (trayectorias subterráneas), porque la altitud mínima es una restricción de
  la Fase 4. `minimum_radius_m` permite detectarlas.
- **Costura periódica en DE:** el operador de mutación de `DE/rand/1/bin@1`
  no usa el arco más corto. Una bielíptica descendente, cuyo primer impulso es
  prógrado, sigue teniéndolo en la costura.
- **Nelder–Mead junto a una frontera de evaluaciones no válidas:** el símplex
  puede colapsar contra ella y estancarse (observado el 2026-11-10: objetivo
  lineal junto a un semiplano no válido, detenido en f = 0,033 con mínimo 0).
  El óptimo de Hohmann está sobre una frontera de ese tipo.
- **No implementado:** control óptimo (RF-3.15), duración del encendido como
  variable (RF-3.06) y optimización multiobjetivo por dominancia (RF-3.22, RF-3.23).
- **Advertencia heredada de la Fase 2:** `hohmann_transfer` emite una
  `RuntimeWarning` y devuelve NaN a escalas extremas (inf − inf fuera de su
  bloque `errstate`). Los adaptadores la suprimen solo en sus propias llamadas
  y rechazan el NaN; la función de la Fase 2 no se ha modificado.
- **Coste por evaluación:** una evaluación válida cuesta alrededor de 1 ms
  con dos impulsos y 1,8 ms con tres (medido el 2026-11-10 en un solo
  proceso), dominada por la sobrecarga de NumPy en arrays pequeños de las
  funciones validadas de la Fase 1, que no se han reescrito.

## Referencias

- Storn, R. y Price, K. (1997). *Differential Evolution: A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces*. Journal of Global Optimization 11.
- Kennedy, J. y Eberhart, R. (1995). *Particle Swarm Optimization*. Proc. IEEE ICNN. Shi, Y. y Eberhart, R. (1998). *A Modified Particle Swarm Optimizer*.
- Clerc, M. y Kennedy, J. (2002). *The Particle Swarm: Explosion, Stability, and Convergence in a Multidimensional Complex Space*. IEEE Trans. Evolutionary Computation 6.
- Deb, K. y Agrawal, R. B. (1995). *Simulated Binary Crossover for Continuous Search Space*. Complex Systems 9. Deb, K. et al. (2002). *A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II*.
- Nelder, J. A. y Mead, R. (1965). *A Simplex Method for Function Minimization*. Lagarias, J. C. et al. (1998). *Convergence Properties of the Nelder–Mead Simplex Method in Low Dimensions*. SIAM J. Optimization 9. McKinnon, K. I. M. (1998). *Convergence of the Nelder–Mead Simplex Method to a Nonstationary Point*.
- Edelbaum, T. N. (1967). *How Many Impulses?* Astronautics & Aeronautics 5.
- Vargha, A. y Delaney, H. D. (2000). *A Critique and Improvement of the CL Common Language Effect Size Statistics*. Holm, S. (1979). *A Simple Sequentially Rejective Multiple Test Procedure*. Wilson, E. B. (1927). *Probable Inference, the Law of Succession, and Statistical Inference*.
- Derrac, J. et al. (2011). *A Practical Tutorial on the Use of Nonparametric Statistical Tests…*. Swarm and Evolutionary Computation 1.

Las referencias se citan de memoria como fuente de los métodos; no se
consultaron en línea durante esta fase. La evidencia de implementación procede
del código y los tests del repositorio.
