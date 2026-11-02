---
tags:
  - validacion
aliases:
  - Criterios de aceptación Fase 1
  - Tests Fase 1
estado: implementado
fase: 1
codigo: tests/
---

# Batería de validación de la Fase 1

Criterios **fijados antes de ejecutar** (`docs/phase1.md`). Son criterios para escala terrestre, no garantías universales del solver.

| Prueba | Criterio | Referencia | Test |
|---|---|---|---|
| Circular analítica, 1 periodo | < 1 m, < 0.001 m/s | Geometría exacta | `test_circular_analytic_motion_and_period_closure` |
| Elipse frente a Kepler | Mismos límites, incl. apsides e instantes intermedios | [[Ecuación de Kepler]] por bisección | `test_eccentric_motion_against_independent_kepler_solution` |
| 10 periodos, inclinada, $e = 0.2$ | Energía y $\mathbf{h}$ relativos < 10⁻⁷; $\mathbf{e}$ absoluto < 10⁻⁷ | [[Invariantes como oráculo de validación\|Invariantes]] | `test_ten_period_invariants_in_inclined_ellipse` |
| Cierre orbital | < 1 m, < 0.001 m/s | [[Tercera ley de Kepler]] | (en el test circular) |
| Conversión de elementos | `rtol=1e-12`, `atol=1e-7 m` y `1e-10 m/s` | Composición identidad | `test_roundtrip_general_and_degenerate_geometry` |
| Apsides | $r_p = a(1-e)$, $r_a = a(1+e)$ | [[Ecuación orbital de las cónicas]] | `test_known_periapsis_and_apsides` |
| Orden RK4 ($T/100, T/200, T/400$) | Orden observado ∈ [3.7, 4.5] | [[Orden de convergencia]] | `test_rk4_fourth_order_convergence` |
| RK4 fino frente a RK45 | < 1 m, < 0.001 m/s | Consistencia mutua | `test_agreement_between_integrators` |
| Costes impulsivos | Opuestos no se cancelan | [[Maniobra impulsiva]] | `test_impulse_costs_do_not_cancel` |
| Entradas inválidas | `ValueError` | — | `test_invalid_*` |
| Fallo del integrador | `RuntimeError`, nunca éxito parcial | — | `test_solver_failure_is_explicit` |
| Experimento end-to-end | Acepta referencia, rechaza configuración estricta/mala | — | `test_experiment.py` |

Total: **20 tests** en 4 ficheros.

Ese total corresponde al baseline `v0.1.0`. La ampliación local
`feature/phase-1-kepler` ejecuta **32 tests en 5 ficheros**: añade 10 tests
de raíces/propagación analítica y 2 de propagación numérica (apsides y
comparación con Kepler). La primera ejecución falló en e=0.95; tras justificar el término de
acondicionamiento del estado inicial en el test, la suite final pasó los
32 tests. RV-1.10 también se amplió a 201 muestras comunes. Ver
[[Cierre documental de la Fase 1]] para el historial.

> [!info] Por qué el error de $\mathbf{e}$ es absoluto
> Con $e = 0$ el error relativo dividiría por cero.

> [!info] Por qué los límites de aceptación son mayores que `atol`
> `atol` es local por paso; los límites de aceptación miden error **global** acumulado. Ver [[Tolerancias justificadas]].

Resultados: [[Resultados de validación de la Fase 1]].

## Contexto posterior — 2026-11-02

Las menciones anteriores a ampliaciones locales o Phase 2 futura son históricas.
Kepler y Phase 2 están integrados en `develop` (`27b4116`); la suite actual
pasó 47 tests. No se reemplazan las cifras, fechas ni el historial del criterio
de Kepler de esta nota. Ver [[Validación de la Fase 2 e integración]].
