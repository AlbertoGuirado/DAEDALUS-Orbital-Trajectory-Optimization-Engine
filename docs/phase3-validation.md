# Fase 3: validación ejecutada

Registro de las ejecuciones del **2026-11-10** sobre la rama
`feature/phase-3-optimization`. Entorno: Windows 11, Python 3.12.1,
NumPy 2.3.3, SciPy 1.16.3. API y formulación: [phase3-api.md](phase3-api.md).

Las cifras de este documento proceden de ejecuciones reales. Son resultados
de las configuraciones indicadas, no cotas universales ni clasificaciones
generales de algoritmos.

## Suite de tests

```powershell
.\.venv\Scripts\python.exe -B -m unittest discover -s tests
```

**266 tests OK, salida 0** (82 s). Corresponden a la Fase 3 los 215 tests de
estos módulos; los 51 restantes son de las fases 1 y 2 y siguen pasando sin
cambios:

| Módulo | Tests | Qué comprueba |
|---|---:|---|
| `test_optimization_problem` | 15 | Contratos de `Bounds`, `Evaluation` y `FunctionProblem` |
| `test_differential_evolution` | 28 | Operadores, orden de extracciones, presupuesto, terminación, reproducibilidad y referencia fijada |
| `test_benchmark_functions` | 8 | Funciones de prueba y guardas de DE (RV-3.02) |
| `test_optimizers` | 32 | Contrato común y operadores de PSO, GA, búsqueda aleatoria, Nelder–Mead e híbrido |
| `test_weighted_objective` | 9 | Objetivo ponderado adimensional (RF-3.19 a RF-3.21) |
| `test_kepler_anomalies` | 10 | Anomalías y tiempo de vuelo |
| `test_arc_extrema_and_bielliptic` | 13 | Extremos de radio por arco y transferencia bielíptica |
| `test_coplanar_transfer_problem` | 15 | Problema de referencia de dos impulsos (RV-3.01, RV-3.04) |
| `test_coplanar_multi_impulse` | 24 | Problema de N impulsos, planes y verificación |
| `test_impulsive_arcs` | 11 | Propagación con arcos, proyección, colores y procedencia |
| `test_benchmark_statistics` | 15 | Estadística frente a cálculo a mano y enumeración exacta |
| `test_benchmarking` | 22 | Registro, ejecutor, informe y experimento |
| `test_export_evidence` | 4 | Exportación de evidencia |
| `test_layering` | 13 | Reglas de capas |

Referencias independientes usadas por los tests (RNF-14):

- **Óptimo de dos impulsos:** forma cerrada de Hohmann y un oráculo cerrado
  escrito en el propio test a partir de vis-viva y la ecuación de Kepler.
- **Tres impulsos:** forma cerrada bielíptica; invariancia del coste ante un
  impulso intermedio nulo; propagación RK45 del plan.
- **Umbrales clásicos de la bielíptica:** el límite biparabólico iguala a
  Hohmann en r2/r1 = 11,94 (localizado por bisección sobre las formas
  cerradas) y la bielíptica gana para todo rb a partir de 15,58.
- **Optimalidad:** por debajo de 11,94 ninguna muestra factible de dos o tres
  impulsos mejora a Hohmann en más de 1e-6 m/s (6000 y 4000 muestras).
- **Extremos de radio:** ecuación de la cónica y muestreo denso de 2001
  puntos en 40 arcos aleatorios.
- **Orden de extracciones:** reimplementación escalar independiente de DE,
  PSO, GA y búsqueda aleatoria.
- **Estadística:** cuantiles y corrección de Holm calculados a mano; A12 y U
  recontados par a par; p de Mann–Whitney por enumeración de todas las
  asignaciones; el intervalo de Wilson, por su ecuación cuadrática.

## Benchmark de cierre de la Fase 3

```powershell
.\.venv\Scripts\python.exe -B -m experiments.benchmark_optimizers `
    --config configs/benchmarks/phase3_differential_evolution.json `
    --output outputs/phase3-benchmark-2026-11-10 --jobs 12 --plots
```

- **Código:** commit `cef41a0`, árbol limpio.
- **Algoritmo:** `DE/rand/1/bin@1` con sus valores de la literatura (NP = 10·D, F = 0,5, CR = 0,9).
- **Semillas:** 0 a 29 (30 ejecuciones por escenario; 210 en total).
- **Duración:** 131 s con 12 procesos.
- **Evidencia versionada:** [evidence/phase3_benchmark.json](evidence/phase3_benchmark.json),
  derivada de `report.json` (SHA-256 `b4dde211…450cb`).

Escenarios, presupuestos, semillas, umbrales y criterios de aceptación se
fijaron en la configuración **antes** de la primera ejecución:

- **Funciones de prueba:** formas estándar en cajas asimétricas, para que el
  óptimo no sea el centro. Presupuesto de 10⁴·D evaluaciones y nivel de
  precisión 1e-8, la convención del benchmark CEC 2005. Las funciones no
  están desplazadas ni rotadas, así que los resultados no son comparables con
  esa colección.
- **Transferencias:** presupuesto de 10³·D evaluaciones, una décima parte,
  porque cada evaluación propaga una trayectoria. Umbral de éxito de 1 m/s (el
  de convergencia práctica del diseño de la Fase 3). Una ejecución mejor que
  la referencia en más de 1e-6 m/s es una anomalía.

### Resultados

Gap = mejor objetivo − referencia. Éxito = gap ≤ umbral y, en las
transferencias, plan verificado por propagación.

| Escenario | D | Presupuesto | Gap mediano | Gap Q1 – Q3 | Peor gap | Éxitos | Intervalo de Wilson 95 % |
|---|---:|---:|---:|---|---:|---:|---|
| `sphere-5d` | 5 | 50 000 | 0 | 0 – 0 | 0 | 30/30 | 88,6 % – 100 % |
| `rosenbrock-5d` | 5 | 50 000 | 4,08e-2 | 4,37e-4 – 2,15e-1 | 1,09 | 6/30 | 9,5 % – 37,3 % |
| `rastrigin-5d` | 5 | 50 000 | 0 | 0 – 0 | 0 | 30/30 | 88,6 % – 100 % |
| `leo-geo` | 2 | 2 000 | 7,91e-4 m/s | 1,22e-4 – 8,08e-3 m/s | 11,5 m/s | 28/30 | 78,7 % – 98,2 % |
| `leo-meo` | 2 | 2 000 | 1,49e-4 m/s | 3,95e-5 – 1,83e-3 m/s | 13,8 m/s | 29/30 | 83,3 % – 99,4 % |
| `geo-leo` | 2 | 2 000 | 1,65e-1 m/s | 4,43e-2 – 1,74 m/s | 35,6 m/s | 21/30 | 52,1 % – 83,3 % |
| `leo-geo-3-impulse` | 5 | 5 000 | 1,44 m/s | 8,62e-1 – 10,5 m/s | 72,1 m/s | 11/30 | 21,9 % – 54,5 % |

Referencias de Hohmann: 3 892,608 m/s (`leo-geo`, `geo-leo`,
`leo-geo-3-impulse`) y 3 459,328 m/s (`leo-meo`). Radios: LEO 6 678 km,
MEO 26 560 km, GEO 42 164 km.

Datos complementarios:

| Escenario | Evaluaciones hasta el éxito (mediana, solo éxitos) | Fracción de evaluaciones no válidas (mediana) | Mejor gap |
|---|---:|---:|---:|
| `sphere-5d` | 4 841 | 0 | 0 |
| `rosenbrock-5d` | 10 120 | 0 | 0 |
| `rastrigin-5d` | 20 280 | 0 | 0 |
| `leo-geo` | 923 | 42,6 % | 7,1e-6 m/s |
| `leo-meo` | 952 | 39,7 % | 1,1e-5 m/s |
| `geo-leo` | 1 266 | 38,4 % | 6,5e-5 m/s |
| `leo-geo-3-impulse` | 4 482 | 51,6 % | 0,29 m/s |

### Criterios de aceptación

| Criterio (fijado en la configuración) | Observado | Resultado |
|---|---|---|
| Ninguna ejecución mejora la referencia en más de 1e-6 m/s (RV-3.01) | 0 anomalías en 210 ejecuciones | Cumple |
| Todo mejor plan pasa la verificación por propagación (RV-3.04) | 120 de 120 planes | Cumple |
| Tasa de éxito de DE en `leo-geo` ≥ 90 % | 93,3 % (28/30) | Cumple |
| Reproducibilidad bit a bit al repetir la primera semilla (RV-3.03) | 7 de 7 parejas idénticas | Cumple |

El proceso terminó con salida 0.

### Lectura de los resultados

- **RV-3.01 se cumple en el sentido del requisito:** en `leo-geo`, DE
  converge al óptimo de Hohmann sin mejorarlo. La mediana queda a 0,8 mm/s y
  28 de 30 semillas a menos de 1 m/s.
- **La tasa de éxito observada tiene incertidumbre:** con 30 semillas, el
  93,3 % de `leo-geo` es compatible con una tasa real de entre el 78,7 % y el
  98,2 %. El criterio del 90 % se cumple en la muestra, no queda demostrado
  para la población de semillas.
- **Rosenbrock-5D no se resuelve de forma fiable:** 6 de 30 semillas alcanzan
  1e-8 con 50 000 evaluaciones. Las demás se detienen en el valle curvo
  (mediana 0,04). Es un resultado de la configuración por defecto, no ajustada.
- **La bajada es más difícil que la subida con el mismo presupuesto:**
  `geo-leo` tiene 21 de 30 éxitos frente a 28 de `leo-geo`, aunque el óptimo
  es el mismo por simetría temporal. La caja de decisión de la bajada es más
  ancha respecto a la velocidad circular de partida.
- **Tres impulsos con 5 000 evaluaciones no bastan:** 11 de 30 éxitos y
  mediana de 1,44 m/s. El problema tiene 5 variables y un conjunto degenerado
  de óptimos (cualquier división del impulso de circularización da el coste
  de Hohmann). Ninguna ejecución mejoró a Hohmann, como exige la teoría para
  r2/r1 = 6,31 < 11,94.
- **Ninguna excepción de la física capturada** (`captured_value_errors` y
  `captured_runtime_errors` a cero) en las 120 ejecuciones de transferencia.

## Verificación de una transferencia larga

Medición del 2026-11-10 que sustenta la regla "endurecer el integrador, no
relajar la tolerancia" de `verify_impulsive_plan`. Bielíptica con r1 = 6 678 km,
r2 = 20·r1 y rb = 60·r1; duración 13,25 días; Δv = 4 023,147 m/s (Hohmann:
4 131,249 m/s).

```python
from orbital_optimizer.problems import CoplanarMultiImpulseProblem, verify_impulsive_plan

LEO = 6_678_000.0
problem = CoplanarMultiImpulseProblem(LEO, 20 * LEO, impulses=3)
plan = problem.plan(problem.bielliptic_decision(60 * LEO))
print(verify_impulsive_plan(plan).metrics)
print(verify_impulsive_plan(plan, rtol=1e-12, atol=[1e-6] * 3 + [1e-9] * 3).metrics)
```

| Integrador | Error de semieje | Error de radio final | Error máximo de radio en la costa | Resultado (tolerancia 1 m) |
|---|---:|---:|---:|---|
| RK45 por defecto (`rtol=1e-10`, `atol` 1e-3 m y 1e-6 m/s) | 1,37 m | 2,10 m | 10,7 m | No pasa |
| RK45 con `rtol=1e-12`, `atol` 1e-6 m y 1e-9 m/s | 0,0035 m | 0,0054 m | 0,028 m | Pasa |

El plan es el mismo: lo que cambia es el error global del integrador de
verificación, que crece con la duración.

## Trazabilidad de requisitos de la Fase 3

| Requisito | Prioridad | Estado | Evidencia |
|---|---|---|---|
| RF-3.01 Variables de decisión con significado, unidades y límites | Obligatorio | Implemented | `Bounds`; tablas de [phase3-api.md](phase3-api.md) |
| RF-3.02 Función objetivo, J = Δv total | Obligatorio | Implemented | Problemas coplanarios |
| RF-3.03 Restricciones, igualdades y desigualdades | Obligatorio | Implemented en parte | La igualdad de llegada se elimina por construcción; los predicados de dominio son desigualdades. Las restricciones de misión son de la Fase 4 |
| RF-3.04 Método de evaluación | Obligatorio | Implemented | Cadena de evaluación documentada |
| RF-3.05 Instantes, magnitud y dirección de los impulsos | Obligatorio | Implemented | `CoplanarMultiImpulseProblem`: los instantes se derivan de las costas M_k; t_1 = 0 por simetría |
| RF-3.06 Duración de los impulsos como variable | Recomendado | No implementado | — |
| RF-3.07 Número de maniobras como parámetro | Recomendado | Implemented | `impulses` ≥ 2; probado con 2 y 3 |
| RF-3.08 Normalización al hipercubo unidad | Obligatorio | Implemented | `Bounds.to_unit` / `from_unit` |
| RF-3.09 Differential Evolution | Obligatorio | Implemented | `DE/rand/1/bin@1` |
| RF-3.10 Interfaz común de optimizador | Obligatorio | Implemented | Protocolo `Optimizer`; seis implementaciones |
| RF-3.11 Algoritmo genético | Recomendado | Implemented | `GeneticAlgorithm` |
| RF-3.12 Enjambre de partículas | Recomendado | Implemented | `ParticleSwarm` |
| RF-3.13 Refinamiento local | Recomendado | Implemented | `NelderMead` |
| RF-3.14 Esquema híbrido | Recomendado | Implemented | `HybridOptimizer`, refinamiento final |
| RF-3.15 Control óptimo | Opcional | No implementado | — |
| RF-3.16 Semilla y resultados idénticos | Obligatorio | Implemented | Tests de reproducibilidad; 7 de 7 en el benchmark |
| RF-3.17 Presupuesto máximo respetado | Obligatorio | Implemented | `EvaluationBudget`; tests de presupuesto exacto |
| RF-3.18 Historia de convergencia | Obligatorio | Implemented | `OptimizationHistory` |
| RF-3.19 Forma ponderada con Δv, tiempo, energía y riesgo | Obligatorio | Implemented | `transfer_objective`; energía y riesgo son indicadores declarados |
| RF-3.20 Todo peso y parámetro documentado | Obligatorio | Implemented | `ObjectiveTerm` exige unidad y justificación |
| RF-3.21 Términos adimensionalizados | Obligatorio | Implemented | Referencias de Hohmann |
| RF-3.22 Optimización por dominancia de Pareto | Recomendado | No implementado | — |
| RF-3.23 Métricas del frente de Pareto | Recomendado | No implementado | — |
| RV-3.01 No significativamente mejor que Hohmann | Obligatorio | Validated | 0 anomalías en 210 ejecuciones; tests de muestreo |
| RV-3.02 Funciones de prueba con óptimo conocido | Obligatorio | Validated | Sphere y Rastrigin 30/30; Rosenbrock 6/30 |
| RV-3.03 Reproducibilidad por semilla | Obligatorio | Validated | Bit a bit, en proceso y entre procesos |
| RV-3.04 Verificación por propagación | Obligatorio | Validated | 120 de 120 planes |

**Criterio de cierre de la Fase 3** ([02-requerimientos.md](02-requerimientos.md) § 12):
"Differential Evolution reproduce el óptimo de Hohmann en el caso donde es el
óptimo conocido". Se cumple en `leo-geo` y `leo-meo` con los criterios
declarados; en `geo-leo` la tasa de éxito es del 70 %.

## Deuda técnica

| ID | Estado a 2026-11-10 | Evidencia |
|---|---|---|
| DT-10 | Resuelto | Estados de evaluación `INVALID_DOMAIN` y `NUMERICAL_FAILURE`; predicados analíticos antes de la física |
| DT-11 | Resuelto | `propagate_impulsive_arcs` y `radius_extrema_on_arc` |
| DT-13 | Resuelto | Los impulsos coincidentes son `INVALID_DOMAIN` |
| DT-14 | Medido | El objetivo usa `propagate_kepler`, no RK45. En un solo proceso, una evaluación válida cuesta 0,85 ms (problema de referencia), 0,99 ms (general, dos impulsos) y 1,76 ms (tres impulsos); una rechazada por un predicado analítico, entre 19 y 29 µs |

## Known Limitations de la evidencia

- Los resultados corresponden a un único entorno (Windows 11, NumPy 2.3.3).
  La CI remota no se ha ejecutado en esta rama.
- Los hiperparámetros son los de la literatura; no se ha estudiado su sensibilidad.
- El benchmark de cierre evalúa un solo algoritmo. La comparación entre
  algoritmos corresponde a la Fase 6A.
- Los tiempos de pared dependen de la máquina y de la carga de los 12 procesos.
- Un éxito es una ejecución dentro del umbral declarado; no demuestra
  optimalidad global.
