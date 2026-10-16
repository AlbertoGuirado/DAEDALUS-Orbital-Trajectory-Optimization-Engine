# Fase 1: auditoría documental y cierre

Revisión local del **2026-10-16**, sobre `fd9d785`, tag `v0.1.0`, rama
`develop` como baseline. Durante la revisión aparecieron cambios de ingeniería
en `feature/phase-1-kepler` aún sin commit, que se documentan por separado.
La revisión de Codex cubre documentación y reproducción de verificaciones;
Codex no modifica algoritmos ni tests. Claude conserva la responsabilidad de ingeniería.

## Dictamen y alcance

Se confirma la validación del baseline de dos cuerpos según el contrato de
[phase1.md](phase1.md): paquete importable, 20 tests superados y experimento
circular aceptado con gráficos. La cobertura está limitada a los casos
enumerados abajo. La ampliación local también cumple la suite de 32 tests.
RV-1.10 compara ahora los métodos en 201 muestras de una órbita circular
y una elipse inclinada, con límites de 1 m y 0.001 m/s.

La ampliación local añade `solve_kepler` y `propagate_kepler`.
**Se confirma su cierre local: 32 tests OK**, tras incorporar al criterio
analítico el acondicionamiento del estado inicial redondeado. El primer
fallo y el cambio justificado de criterio se conservan en el registro de
validación; no se atribuye este resultado al tag anterior.
Hohmann, aplicación de burns, optimizadores, restricciones de
misión y perturbaciones siguen sin implementar. Los requisitos recomendados y opcionales no se convierten
en obligatorios para declarar cerrada esta entrega.

## Inventario contrastado

| Área | Evidencia en el checkout |
| --- | --- |
| Física | `src/orbital_optimizer/dynamics/two_body.py`: campo gravitatorio, constante y validación |
| Geometría e invariantes | `src/orbital_optimizer/orbital/{elements,mechanics}.py`: conversiones, energía, vectores, periodo y costes |
| Kepler local | `orbital/kepler.py`, `simulation/kepler_propagator.py`: solver y propagación f/g, probados localmente |
| Simulación | `src/orbital_optimizer/simulation/propagator.py`: RK4/RK45 y resultado |
| Visualización científica | `src/orbital_optimizer/visualization/plots.py`: dos PNG opcionales |
| Compatibilidad | `src/physics/`: adaptadores de imports, sin duplicar algoritmos |
| Reproducción | `experiments/validate_two_body.py`, `configs/earth_orbit.json`, `requirements-reproducible.txt` |
| Pruebas | Baseline: 20 tests. Ampliación: 32 tests en cinco módulos (10 nuevos de Kepler y 2 adicionales de propagación) |
| CI | `.github/workflows/tests.yml`: matriz de tests y experimento; ejecuciones remotas no consultadas |
| Reservas | `optimization/`, `missions/`, `agent/`: archivos vacíos de fases posteriores |
| Conocimiento | `../DAEDALUS-Vault/`: inicio, cinco áreas temáticas y plantilla; notas en español con wikilinks |

Al inicio había archivos ajenos sin seguimiento: `docs/informe-visual.html`,
`vault/`, `visualization/` y dos imágenes `world map*.jpg`. Se preservan como
trabajo local fuera del cierre del paquete científico. El vault de conocimiento
usado en esta revisión es el directorio hermano `DAEDALUS-Vault`, conservando
su estructura, sin fusionarlo con `orbital-optimizer/vault/`.

## Trazabilidad de requisitos funcionales

Identificadores originales: [02-requerimientos.md](02-requerimientos.md).

| Requisitos | Implementación y evidencia |
| --- | --- |
| RF-1.01–04 | `two_body.py`; `DynamicsTests`: derivada circular, μ personalizado, entradas inválidas y copia del estado |
| RF-1.05–07 | `propagate`; `PropagatorTests`: circular RK4/RK45, muestreo, intervalo pequeño y μ personalizado |
| RF-1.08–09, RF-1.14 | `elements.py`; `ElementsTests`: apsides, geometría general y degenerada, rechazo de estados no soportados |
| RF-1.10–13, RF-1.15 | `mechanics.py`; `DynamicsTests`: cantidades circulares, periodo, μ y costes sin cancelación |
| RF-1.16 (recomendado) | `solve_kepler` local: residuos, valores especiales y bisección independiente; tests propios superados en la ejecución de 32 tests |
| RF-1.17 (opcional) | `propagate_kepler` local: f/g cartesianos; tests propios superados, incluido e=0.95 con modelo de error explícito |
| RF-1.18–20 | Dominio, ausencia de colisiones y efecto del muestreo RK4 documentados en `phase1.md`, API y docstrings |

## Trazabilidad de validación

Los nombres de funciones siguientes están en los cinco archivos
`tests/test_{dynamics,elements,propagator,experiment,kepler}.py`.

| Requisitos | Test y alcance comprobado |
| --- | --- |
| RV-1.01–03 | `test_ten_period_invariants_in_inclined_ellipse`: 10 periodos, `a=9e6 m`, `e=0.2`, `i=0.7 rad`; tres derivas <1e-7 con RK45 |
| RV-1.04–05 | `test_circular_analytic_motion_and_period_closure`: RK4 y RK45, 101 muestras, error <1 m y <0.001 m/s; incluye extremos |
| RV-1.06, RV-1.12 | `test_roundtrip_general_and_degenerate_geometry`: e=0, 1e-9, 0.2, 0.8; seis inclinaciones, incluida retrógrada ecuatorial; estado reconstruido y ángulos generales |
| RV-1.07 | `test_known_periapsis_and_apsides`: vis-viva en ambos apsides y orientación conocida |
| RV-1.08, RV-1.11 (baseline) | `test_eccentric_motion_against_independent_kepler_solution`: `a=18e6 m`, `e=0.6`, 121 muestras; bisección independiente, ambos apsides, posición <1 m y velocidad <0.001 m/s |
| RV-1.09 | `test_rk4_fourth_order_convergence`: pasos T/100, T/200 y T/400 sin muestreo; orden estrictamente entre 3.7 y 4.5 |
| RV-1.10 | `test_agreement_between_integrators`: máximos de diferencia RK4/RK45 en 201 muestras comunes, circular y elipse inclinada, <1 m y <0.001 m/s |
| RV-1.13 | `test_invalid_*` en dinámica, elementos y propagación |
| RV-1.14 | `test_sampling_copies_and_small_interval`: arrays de solo lectura y copia de `final_state` |
| Fallo explícito | `test_solver_failure_is_explicit`: fallo simulado de `solve_ivp` produce `RuntimeError` |
| Experimento | `test_reference_experiment_and_strict_failure`, `test_bad_configuration_rejected`: aceptación, incumplimiento y rechazo de configuración |

Ingeniería añadió `test_propagated_apsides_in_inclined_ellipse` y
`test_numerical_against_analytic_kepler_propagation` (ambos pasaron),
y `tests/test_kepler.py` con 10 tests. La ejecución final pasó los 32.
Ingeniería también amplió la comparación directa RK4/RK45; Codex no alteró
el código ni los criterios de los tests. Que la suite pase no demuestra todos los
posibles fallos de integradores, escalas numéricas o geometrías orbitales.

## Brechas y correcciones documentales

| Prioridad | Hallazgo inicial | Acción de esta revisión |
| --- | --- | --- |
| P0 | Documento 01 afirmaba que el paquete no arrancaba, sin dinámica ni tests | Inventario y estado actualizados; auditoría histórica resumida como resuelta |
| P0 | Documento 02 mantenía DT-01–09 como defectos actuales | Estados y decisiones explícitos, conservando identificadores |
| P0 | Documento 03 atribuía J2 y arrastre a la implementación | Etiquetados como planificación de Fase 5 |
| P1 | Faltaban guía de API, formato de datos y rutas de lectura | `docs/README.md` y `phase1-api.md` |
| P1 | Cierre sin reproducción reciente ni matriz de requisitos | Esta auditoría y actualización de validación/vault |
| P2 | Cierre podía interpretarse como validación de todo RF-1 | Baseline separado de la ampliación local y del historial de su validación |

La estructura mínima conserva los documentos 01–03, el contrato de Fase 1 y
su registro de validación. Añade un índice, una guía de API y esta auditoría;
mantiene la teoría desarrollada en el vault.

## Reproducción y comprobaciones de cierre

Desde la raíz del repositorio, con el entorno existente:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/earth_orbit.json --output outputs/phase1-closeout-2026-10-16 --plots
```

Primero se reprodujo el baseline con **20 tests OK**. Tras los cambios de
ingeniería, la suite ampliada terminó con código 0 y **32 tests OK**. El
experimento se repitió sobre la ampliación con código 0,
`validation_passed=true`, cinco comprobaciones positivas y 788 evaluaciones.
Resultados y entorno: [phase1-validation.md](phase1-validation.md).
También se ejecutó el ejemplo Python de la guía de API y se comprobaron
imports públicos, enlaces Markdown locales y wikilinks del vault.
Resultado: enlaces locales y wikilinks correctos en 8 documentos técnicos
y 72 notas; ejemplo e imports públicos correctos; JSON aceptado, CSV de
1001 filas y cabeceras PNG verificadas.
No se reconstruyó el wheel ni se verificó CI remota en este cierre.

## Entrega a Claude y siguiente fase

El fallo inicial de e=0.95 quedó resuelto en la validación final mediante
un modelo de error de test que incluye el acondicionamiento del estado
inicial redondeado. El código del propagador analítico no cambió entre
ambas ejecuciones. La comparación RV-1.10 también se amplió y pasó.
No queda una brecha documental conocida para cerrar esta ampliación local.

El siguiente hito de ingeniería es Fase 2: definir aplicación de impulsos y
baseline Hohmann, comprobando la órbita final objetivo (RV-2.01). Las
funciones actuales de coste no aplican burns. Claude debe entregar a Codex
API, ecuaciones, unidades, parámetros, tests y resultados reproducibles de
ese hito. Este archivo es una entrega local; no implica envío de mensajes.

El solver y el propagador analítico son ampliaciones de ingeniería locales
aún sin commit; el tag `v0.1.0` no incluye estas capacidades.
La ausencia de CI remota y de rebuild reciente limita la evidencia de esta
revisión a la instalación local existente.

## Referencias consultadas

Se comprobaron los enlaces oficiales de
[SciPy solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html),
[NAIF oscelt](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/oscelt_c.html)
y [IERS TN36](https://iers-conventions.obspm.fr/content/tn36.pdf).
El enlace estable de SciPy actualmente muestra otra versión; los contratos
de este checkout se contrastaron directamente con código, docstrings y la
instalación local SciPy 1.16.3. No se afirma una comparación numérica con NAIF
ni una auditoría bibliográfica completa del material teórico de fases futuras.
