# Auditoría documental de integración — Phase 1 y Phase 2

Revisión del **2026-11-02**, sobre `develop`, código `27b4116`.
Kepler consta en `b11150f`; Phase 2 en `3df6f11`; integración en `16c5ac0`.
El tag `v0.1.0` conserva el baseline `fd9d785`. No se han modificado código,
tests, configuración, constantes, tolerancias ni criterios de aceptación.

## Estado y trazabilidad

| Estado | Capacidad | Código y evidencia |
| --- | --- | --- |
| Implemented | Dos cuerpos, invariantes y elementos | `dynamics/two_body.py`, `orbital/{mechanics,elements}.py`; tests de dinámica y elementos |
| Implemented | RK4/RK45 | `simulation/propagator.py`; 11 tests de propagación |
| Implemented | Kepler y coeficientes f/g | `orbital/kepler.py`, `simulation/kepler_propagator.py`; 10 tests de Kepler |
| Implemented | Hohmann y propelente ideal | `orbital/transfers.py`, `orbital/mechanics.py`; 9 tests de transferencias |
| Implemented | Impulsos inerciales instantáneos | `simulation/impulsive.py`; 6 tests de impulsos |
| Implemented | Experimento circular | `experiments/validate_two_body.py`; 2 tests de aceptación/rechazo |
| Planned | Bielíptica, cambio de plano, Lambert | Requisitos recomendados/opcionales, sin API |
| Planned | Optimización, restricciones, perturbaciones, planificador | Reservas vacías; Phase 3 aún no iniciada |

Las rutas de código son relativas a `src/orbital_optimizer/`, salvo donde
se indica otra raíz. La suite suma 47 tests: 32 de Phase 1 y 15 de Phase 2.

## Comprobaciones ejecutadas

Desde la raíz del repositorio, con el entorno existente:

```powershell
.\.venv\Scripts\python.exe -B -m unittest discover -s tests -v
```

**47 tests OK, salida 0.** Python 3.12.1, NumPy 2.3.3, SciPy 1.16.3,
Windows 11. Se verificó que el import procede de `src/orbital_optimizer`.
No se hizo instalación limpia, reconstrucción de wheel ni consulta de CI remota.

El siguiente ejemplo reproduce las mediciones de ambas fases. Se ejecutó
su equivalente en memoria; no se sobrescribieron los artefactos históricos.

```python
import json
from pathlib import Path
import numpy as np
from orbital_optimizer import (
    EARTH_MU, hohmann_transfer, propellant_mass,
    propagate_with_impulses, state_to_elements,
)
from experiments.validate_two_body import run_experiment

report, _, _ = run_experiment(
    json.loads(Path("configs/earth_orbit.json").read_text())
)
print(report["metrics"], report["checks"], report["function_evaluations"])

r1, r2 = 6_678_000.0, 42_164_000.0
h = hohmann_transfer(r1, r2)
y0 = np.array([r1, 0, 0, 0, np.sqrt(EARTH_MU/r1), 0])
dv1 = np.array([0, h.first_delta_v, 0])
arrival = propagate_with_impulses(y0, [(0.0, dv1)], h.transfer_time)
dv2 = arrival[3:]/np.linalg.norm(arrival[3:])*h.second_delta_v
final = propagate_with_impulses(
    y0, [(0.0, dv1), (h.transfer_time, dv2)], h.transfer_time
)
elements = state_to_elements(final)
print(h, propellant_mass(2000.0, h.total_delta_v, 320.0))
print("position_error_m", np.linalg.norm(final[:3]-[-r2, 0, 0]))
print("velocity_error_m_s", np.linalg.norm(final[3:]-[0, -np.sqrt(EARTH_MU/r2), 0]))
print("eccentricity", elements.eccentricity)
print("semi_major_axis_error_m", elements.semi_major_axis-r2)
print("without_second_impulse_e", state_to_elements(arrival).eccentricity)
```

### Phase 1: caso circular configurado

Radio 7 000 000 m, inclinación 28,5°, un periodo, 1001 muestras,
RK45 y parámetros originales de `configs/earth_orbit.json`.

| Magnitud | Observado | Criterio de la configuración |
| --- | ---: | ---: |
| Máximo error de posición muestreado | 0.015833654949214758 m | 1 m |
| Máximo error de velocidad muestreado | 1.59444751686348e-5 m/s | 0.001 m/s |
| Deriva relativa de energía | 9.429883984432415e-10 | 1e-7 |
| Deriva relativa del vector h | 4.714942585888468e-10 | 1e-7 |
| Deriva absoluta del vector e | 1.002666319192962e-9 | 1e-7 |

`validation_passed=true`, cinco checks verdaderos y 788 evaluaciones.
Las métricas coinciden con el registro anterior; eso no convierte sus fechas
o tiempos de ejecución en mediciones nuevas. No se generaron gráficos nuevos.

### Phase 2: LEO→GEO circular coplanaria ideal

Radios 6 678 000 m y 42 164 000 m; `mu=EARTH_MU`, masa inicial 2000 kg,
Isp 320 s; wrapper RK45 con sus valores por defecto.

| Magnitud | Observado |
| --- | ---: |
| Primer impulso | 2425.7690283068578 m/s |
| Segundo impulso | 1466.8387152844527 m/s |
| Δv total | 3892.6077435913103 m/s |
| Tiempo | 18990.051838481286 s |
| Semieje de transferencia | 24421000 m |
| Propelente ideal | 1421.4766255239515 kg |
| Norma del error final de posición | 0.07199410476581594 m |
| Norma del error final de velocidad | 1.4275188021570504e-5 m/s |
| Excentricidad final | 3.7036434861574957e-9 |
| Error firmado del semieje final | 0.06418003886938095 m |
| Excentricidad sin segundo impulso | 0.7265468247853418 |

Las normas son mediciones informativas, no nuevas aserciones. El test de
llegada mantiene `delta=1 m` para radio y semieje, `e<1e-6`, y
`assert_allclose` de componentes con `atol=1 m` y `1e-3 m/s`; su `rtol`
es el predeterminado de NumPy, no cero. La comparación cerrada de impulsos
usa `assertAlmostEqual(..., places=9)` y la de tiempo `places=6`.
El test de masa final usa `rtol=1e-12` con `atol` predeterminado.
Estos detalles no equivalen a cotas universales de error absoluto.

## Known Limitations y discrepancias pendientes

- Estados elípticos no radiales para los arcos; sin colisiones, superficie,
  perturbaciones, empuje finito, masa dinámica, época ni efemérides.
- Hohmann devuelve magnitudes. Al invertir radios cambia el orden de los
  impulsos. Con radios iguales conserva medio periodo aunque el coste sea cero.
- El wrapper devuelve solo el estado final. El impulso final no recibe una
  nueva validación orbital; las escalas extremas de Hohmann tampoco tienen
  una comprobación explícita de finitud de todos los Δv de salida.
- El docstring de `HohmannTransfer` dice que el resultado es el mismo en
  ambos sentidos. La documentación precisa la simetría de total/tiempo y el
  intercambio de impulsos; se deja el docstring intacto por el alcance autorizado.
- La prueba de inversión comprueba total y tiempo, no una propagación completa
  descendente. No se infiere esa cobertura ni se añaden tests.
- No se verificaron CI remota, ruedas construidas ni la rama de UI. DT-14
  sigue siendo una hipótesis de coste, sin benchmark poblacional reproducido.
- El registro histórico de fallo y revisión del modelo de error de Kepler
  se conserva. No se cambió ninguna tolerancia en esta actualización.

Estas observaciones constituyen un registro para ingeniería, no una corrección
del código ni una afirmación de haber enviado mensajes a otro agente.

## Cambios documentales

Se reconciliaron índice, mapa, informe visual y documentos 01–03; se
precisaron contratos de ambas fases y la condición histórica de sus cierres.
Se actualizaron el flujo de datos y la nota de contribución obsoleta.
El vault conserva sus áreas, frontmatter y wikilinks; amplía derivaciones y
separa teoría, implementación, resultados y trabajo futuro.

Referencias comprobadas: NASA Glenn para el cohete ideal, Orbital Mechanics
& Astrodynamics para Hohmann y movimiento relativo, y documentación oficial
de SciPy para `solve_ivp`. El enlace estable de SciPy muestra una versión más
reciente; el enlace web específico 1.16.3 no fue accesible mediante la herramienta.
El comportamiento de esta instalación se contrastó con su código y tests.
No se realizó una auditoría bibliográfica exhaustiva de la teoría de fases futuras.

## Inventario de entrega y control documental

Archivos existentes actualizados:

- `CONTRIBUTING.md`: elimina la descripción obsoleta de Hohmann ausente.
- `docs/README.md`: estado, navegación y enlace a la copia versionada del vault.
- `docs/mapa-del-proyecto.md`: módulos, tests y documentos de Phase 2.
- `docs/informe-visual.html`: 18 nodos y 29 relaciones; añade Kepler e impulsos,
  actualiza evidencia y distingue Implemented, Planned y Known Limitations.
- `docs/01-caracteristicas-del-proyecto.md`: arquitectura y estado coherentes.
- `docs/02-requerimientos.md`: trazabilidad actual sin cambiar requisitos.
- `docs/03-fundamentos-cientificos-y-algoritmicos.md`: estado, signos y alcance
  de las afirmaciones sobre Hohmann y futuros optimizadores.
- `docs/phase1.md` y `docs/phase1-api.md`: contrato actual y frontera entre fases.
- `docs/phase1-validation.md` y `docs/phase1-closeout.md`: aviso histórico;
  se conservan sus mediciones y el historial de validación.
- `docs/phase2-api.md`: ecuaciones, unidades, supuestos, APIs y limitaciones.
- `docs/data-flow-mission-planning.md` y `docs/mission-planning-interaction.html`:
  contexto de integración actualizado. El SVG mantiene su flujo, que sigue vigente.

Archivo nuevo: `docs/phase1-phase2-audit.md` (este informe).

Se incorpora una **copia documental de 73 notas Markdown** en
`DAEDALUS-Vault/`. El vault original, hermano del repositorio, se conserva;
no se incluyen preferencias `.obsidian`. Los enlaces relativos del cierre
de Phase 1 se adaptan en la copia para apuntar al `docs/` versionado.
La sincronización es manual, sin automatización añadida.

Notas revisadas del vault: `00 Inicio`; `MOC Física orbital`; `Problema de dos
cuerpos`; `Ecuación de Kepler`; `Transferencia de Hohmann`; `Maniobra impulsiva`;
`Ecuación de Tsiolkovsky`; `Newton-Raphson`; `Método de bisección`;
`Differential Evolution`; `Arquitectura del código`; `Roadmap por fases`;
`Cierre documental de la Fase 1`; `Resultados de validación de la Fase 1`;
`Batería de validación de la Fase 1`; `MOC Proyecto DAEDALUS` y
`MOC Validación y benchmark`. Nota nueva: `Validación de la Fase 2 e integración`.
Las demás notas se copian para preservar navegación y contexto; incorporarlas
a Git no implica haber validado exhaustivamente su teoría.

Comprobaciones documentales: 86 archivos Markdown de docs/vault, 75 enlaces
locales y 421 wikilinks sin destinos ausentes; los alias escapados de tablas
se interpretaron como tales. Se ejecutaron los ejemplos de importación/uso
de `phase1-api.md` y el ejemplo de reproducción de este informe. Se comprobó
la sintaxis JavaScript con `node --check`, unicidad y referencias de los
18 nodos/29 relaciones y sus estados. `git diff --check` sin errores.
No se hizo una comprobación visual en navegador ni de cada interacción HTML.
