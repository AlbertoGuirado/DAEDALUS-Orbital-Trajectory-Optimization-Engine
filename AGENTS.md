# DAEDALUS — contexto para agentes

## Proyecto y fuente de verdad

DAEDALUS (Orbital Trajectory Optimization Engine) es un proyecto científico
Python de dinámica orbital y futura optimización matemática de transferencias.
El núcleo es físico y numérico; no hay un agente de IA ni un planificador de
misiones implementado.

Este archivo resume el contexto revisado el **2026-11-03**. Antes de trabajar,
comprueba la rama, `git status --short --branch`, el código y los tests:
esta instantánea no sustituye la evidencia del checkout.
Lee `README.md`, `CONTRIBUTING.md` y los documentos de la fase afectada.
Conserva los cambios del usuario y limita las modificaciones a la tarea pedida.

## Estado de las fases

- **Implemented — Phase 1:** gravedad de dos cuerpos, invariantes, elementos
  orbitales clásicos y conversiones, RK4/RK45, `solve_kepler` y propagación
  analítica cartesiana mediante coeficientes f/g (`propagate_kepler`).
- **Implemented — Phase 2:** Hohmann circular coplanaria (`hohmann_transfer`),
  propelente ideal por Tsiolkovsky (`propellant_mass`) y aplicación de impulsos
  instantáneos (`propagate_with_impulses`).
- **Planned:** bielíptica, cambio de plano y Lambert; Phase 3 (Differential
  Evolution), Phase 4 (restricciones), Phase 5 (perturbaciones), Phase 6A
  (benchmark científico) y Phase 6B (orquestación de misiones).
- **Evidencia histórica reciente:** 47 tests OK y casos circular y LEO→GEO
  reproducidos el 2026-11-02. Ver `docs/phase1-phase2-audit.md`. No atribuyas
  esos resultados a una ejecución nueva ni a configuraciones no ensayadas.

`develop` contiene la integración de fases 1–2. `v0.1.0` identifica el baseline
histórico de Phase 1, no todas las capacidades actuales. La versión declarada
del paquete no demuestra el alcance de un checkout o de un wheel local.

## Mapa del repositorio

| Ruta | Responsabilidad |
| --- | --- |
| `src/orbital_optimizer/dynamics/two_body.py` | Campo gravitatorio, μ terrestre y validaciones básicas |
| `src/orbital_optimizer/orbital/` | Invariantes, elementos, Kepler, Hohmann y cálculo de propelente |
| `src/orbital_optimizer/simulation/` | Propagación RK4/RK45, Kepler y arcos con impulsos |
| `src/orbital_optimizer/visualization/` | Gráficos científicos opcionales con Matplotlib |
| `src/physics/` | Adaptadores de imports históricos; no duplicar física |
| `tests/` | `unittest`, referencias analíticas y contratos de entrada/salida |
| `experiments/validate_two_body.py` | Validación reproducible de una órbita circular |
| `configs/earth_orbit.json` | Caso de referencia y criterios de aceptación |
| `configs/experiments/` | Escenarios exploratorios, distintos del caso de referencia |
| `docs/` | Contratos, API, requisitos, evidencia y mapas del proyecto |
| `DAEDALUS-Vault/` | Copia versionada de las notas científicas Markdown |
| `optimization/`, `missions/`, `agent/` | Reservas históricas vacías, no algoritmos disponibles |

La física no debe depender de simulación u optimización. Las implementaciones
nuevas pertenecen al paquete `orbital_optimizer`, respetando esas capas.
La API pública se expone desde `src/orbital_optimizer/__init__.py`.
La UI Streamlit se ha desarrollado aparte, en `feature/streamlit-ui`;
no asumir que está disponible en la rama actual.

## Convenciones científicas y Known Limitations

- SI: posición/radio en m, velocidad y Δv en m/s, tiempo en s, masa en kg,
  μ en m³/s². Ángulos internos en radianes. Radio no significa altitud.
- Estado `(6,)`: `[x, y, z, vx, vy, vz]`, marco inercial dextrógiro centrado
  en el primario; sin época, efemérides ni orientación de catálogo.
- Los propagadores trabajan con elipses ligadas no radiales. La validación
  básica del campo gravitatorio tiene un dominio más amplio; no confundirlas.
- No existen superficie, colisiones, arrastre, J2, terceros cuerpos, empuje
  finito ni masa integrada en el estado. Calcular propelente no integra masa.
- `propagate` es vuelo libre. El wrapper impulsivo suma vectores de velocidad
  inercial entre arcos, exige tiempos estrictamente crecientes y devuelve solo
  el estado final. No revalida orbitalmente un impulso aplicado al final.
- Hohmann devuelve magnitudes, no vectores orientados. Invertir radios
  intercambia las magnitudes de los impulsos y conserva total y tiempo.
  Con radios iguales mantiene medio periodo aunque Δv sea cero.
- `total_delta_v` suma normas; no calcula la norma de la suma de vectores.
- Tolerancias locales de RK45 no son cotas del error global. Separar errores
  de posición y velocidad. Cerca de singularidades de elementos, contrastar
  estados reconstruidos antes de interpretar ángulos mal condicionados.
- La propagación analítica conserva errores de raíz y redondeo. Los invariantes
  de vuelo libre pueden cambiar en cada impulso. No afirmar precisión universal.
- No relajar tolerancias, retirar aserciones ni alterar parámetros físicos
  para hacer pasar una validación. Toda revisión científica requiere justificación
  explícita y trazabilidad; una tarea documental no autoriza esos cambios.

## Entorno y comprobaciones

Python >=3.11; núcleo NumPy/SciPy; Matplotlib opcional. Usa el entorno del
proyecto. Desde la raíz del repositorio, en PowerShell:

```powershell
# Solo si hace falta preparar el entorno:
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[plots]"

# Suite existente:
.\.venv\Scripts\python.exe -B -m unittest discover -s tests -v

# Experimento; elegir una carpeta nueva para preservar resultados previos:
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/earth_orbit.json --output outputs/<ejecucion-nueva> --plots
```

Sustituye `<ejecucion-nueva>` por un nombre de carpeta antes de ejecutar.
Para el núcleo sin gráficos instala `-e .` y omite `--plots`.
En Linux/macOS usa `.venv/bin/python`.

Ejecuta comprobaciones proporcionales al cambio: tests pertinentes para
ingeniería; enlaces, ejemplos y diff para documentación. Si cambia la física
o propagación, reproduce el experimento. No declares éxito sin observarlo.
Los experimentos actuales son deterministas; no requieren seeds.
`outputs/`, entornos, cachés, `build/` y `dist/` son artefactos regenerables,
no evidencia versionada automáticamente. No sobrescribas registros históricos.

## Documentación y conocimiento

Puntos de entrada: `docs/README.md`, `docs/mapa-del-proyecto.md`,
`docs/phase1.md`, `docs/phase1-api.md`, `docs/phase2-api.md` y
`docs/phase1-phase2-audit.md`. `docs/informe-visual.html` presenta el sistema;
no ejecuta el motor orbital.

Mantén explícitos **Implemented**, **Planned** y **Known Limitations**.
Distingue código existente, tests ejecutados, validación de un caso y teoría.
Documenta ecuaciones, unidades, supuestos, dominio y evidencia. No inventes
resultados, bibliografía, benchmarks ni optimalidad de misiones generales.
Conserva las fechas y el historial de los cierres anteriores.

El vault usa español, frontmatter y wikilinks. En el espacio de trabajo del
usuario existe también `../DAEDALUS-Vault/`, fuera del repositorio. La copia
versionada contiene solo Markdown, sin `.obsidian`; no hay sincronización
automática. Si una tarea afecta ambas copias, actualízalas de forma coherente
y adapta sus enlaces relativos. No asumas que el vault externo existe en CI.

Para tareas documentales, limita los cambios a documentación y registra los
posibles defectos del código como discrepancias para ingeniería. No implementes
funcionalidades ni cambies tests bajo una actualización documental.


## Documentation Delegation

For documentation-related tasks, use specialized agents.

- Internal project documentation:
  delegate to `daedalus_docs`.

- Scientific knowledge and reference documentation:
  delegate to `daedalus_vault`.

- Tasks involving both documentation domains:
  delegate to both agents in parallel when their work
  can be performed independently.

The main agent coordinates the work and reviews
the final output.

Never delegate scientific implementation changes
to documentation agents.

## Git y entrega

Sigue `CONTRIBUTING.md`: `main` estable, `develop` integración y ramas
`feature/*` o `fix/*` para trabajo acotado. No cambies de rama ni reescribas
historial sin autorización; nunca descartes cambios ajenos.
Haz commit/push cuando la tarea o la conversación lo autoricen.

Al cerrar, informa de archivos y comportamiento cambiados, comprobaciones
realmente ejecutadas, resultados y límites pendientes. Si hay commit/push,
indica rama y commit. No confundas configuración de CI con una ejecución aprobada.

