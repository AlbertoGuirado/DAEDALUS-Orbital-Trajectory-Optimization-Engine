# Fase 4: restricciones de misión y factibilidad

Estado a **2026-11-10**, rama `feature/phase-4-constraints`. Cubre RF-4.01 a
RF-4.10, RV-4.01 y RV-4.02 de [02-requerimientos.md](02-requerimientos.md).
Teoría: [03-fundamentos](03-fundamentos-cientificos-y-algoritmicos.md) § 15.
Unidades SI en todas las magnitudes.

Cobertura funcional: [tests/test_spacecraft_and_constraints.py](../tests/test_spacecraft_and_constraints.py)
y [tests/test_constraint_handling.py](../tests/test_constraint_handling.py).
Configuración reproducible del benchmark:
[phase4_constraint_handling.json](../configs/benchmarks/phase4_constraint_handling.json).

## Qué hace esta fase y qué no

**Implemented:**

- un modelo de nave (`Spacecraft`) con masa seca, propelente, empuje máximo e
  impulso específico (RD-03);
- un catálogo de seis restricciones de misión, cada una con su expresión,
  unidad, signo, tolerancia y magnitud de la violación;
- tres estrategias de manejo de restricciones para los optimizadores;
- un adaptador que aplica las restricciones a un problema de transferencia;
- la verificación independiente de las restricciones mediante propagación numérica;
- una solución de referencia independiente para el óptimo con límite de tiempo.

**No pertenece a esta fase:** J2, arrastre, encendidos de duración finita y
masa variable integrada (Fase 5); selección de estrategia de misión (Fase 6B).

## Capas

| Módulo | Capa | Puede importar |
|---|---|---|
| `orbital_optimizer/spacecraft.py` | Física | Constantes y validación básica |
| `orbital_optimizer/constraints/` | Física | Constantes, nave y simulación |
| `orbital_optimizer/optimization/constraint_handling.py` | Optimización genérica | NumPy y biblioteca estándar |
| `orbital_optimizer/problems/constrained.py`, `references.py` | Adaptadores | Física y optimización |

`constraints/` no importa `optimization`, `problems`, `benchmarking` ni
`mission`. Qué hace un optimizador con un candidato infactible es un asunto
aparte del modelo físico de restricciones. `tests/test_layering.py` lo comprueba.

## Nave (`orbital_optimizer.spacecraft`)

```python
Spacecraft(dry_mass_kg, propellant_mass_kg, max_thrust_n, specific_impulse_s)
```

Con velocidad efectiva de escape \(v_e = I_{sp}\,g_0\):

| Magnitud | Expresión | Unidad |
|---|---|---|
| Masa húmeda | \(m_0 = m_{seca} + m_{prop}\) | kg |
| Capacidad ideal de Δv | \(v_e \ln(m_0/m_{seca})\) | m/s |
| Propelente de una secuencia | \(m_0\,(1 - e^{-\sum\Delta v_k / v_e})\) | kg |
| Propelente del impulso k | \(m_{k-1} - m_k\), con \(m_k = m_{k-1}\,e^{-\Delta v_k/v_e}\) | kg |
| Duración de un encendido a empuje F | \(t_k = (m_{k-1} - m_k)\,v_e / F\) | s |

- Se distinguen tres masas: la seca excluye el propelente, la de propelente
  es la de los depósitos y la húmeda es su suma al inicio.
- La duración sale del gasto másico constante \(F/v_e\). Es la de un
  encendido ideal a empuje constante, sin pérdida por gravedad. El modelo
  impulsivo no la usa: indica cuánto dista un impulso de ser instantáneo.
- `propellant_per_burn_kg` y `propellant_required_kg` no limitan la petición
  a lo que hay en los depósitos. Pedir más es una restricción violada, no un
  error de entrada.
- Se supone un único impulso específico para todos los impulsos.
- Los valores físicamente imposibles producen `ValueError` (RD-05).

## Convención de signos y tolerancia (`constraints.model`)

Toda restricción se escribe **g ≤ 0**. El residuo g está en la unidad de la
magnitud restringida y es positivo en la cantidad en que se excede el límite.

| Concepto | Definición | Unidad |
|---|---|---|
| Satisfecha | g ≤ tolerancia; un g no finito nunca lo está | — |
| Violación | max(g, 0) si no está satisfecha; 0 si lo está | La de la magnitud |
| Violación normalizada | violación / escala | Adimensional |
| Margen | −g (positivo: queda holgura) | La de la magnitud |

- La **escala** es una magnitud positiva de referencia (normalmente el propio
  límite). Solo sirve para que violaciones de unidades distintas se puedan sumar.
- La **tolerancia** es absoluta y explícita. Vale
  `relative_tolerance` × escala, con `relative_tolerance = 1e-9` por defecto.
- Cuando una restricción no se satisface, la violación es el residuo
  completo, no el residuo menos la tolerancia.
- Un residuo NaN da violación infinita: nunca hay falsos factibles por
  valores no finitos.

**Por qué 1e-9.** Las magnitudes evaluadas arrastran un redondeo de hasta
1e-10 relativo (condicionamiento de la propagación analítica). Con 1e-9, un
candidato situado exactamente sobre un límite no se rechaza por redondeo.
Equivale a 4e-6 m/s en un presupuesto de 3,9 km/s y a 7 mm en un radio LEO,
muy por debajo de la exactitud del modelo físico. Es un parámetro y queda
registrado en cada informe.

`ConstraintReport` reúne las restricciones evaluadas de un candidato:

- `feasible`: todas las restricciones **duras** satisfechas;
- `total_violation`: suma de las violaciones normalizadas de las duras;
- `soft_violation`: lo mismo para las blandas;
- `violated`, `by_name(name)`, `normalized()`, `to_dict()`;
- `not_applicable`: nombre → motivo, para restricciones que el modelo físico
  no puede evaluar. Se informan; nunca se dan por satisfechas en silencio.

Las restricciones **blandas** (`MissionLimits.soft`) se evalúan e informan
igual, pero no deciden la factibilidad. Por defecto todas son duras.

## Catálogo de restricciones (`constraints.limits`)

```python
MissionLimits(max_wet_mass_kg=None, max_thrust_n=None, max_delta_v_m_s=None,
              max_duration_s=None, min_altitude_m=None,
              body_radius_m=EARTH_EQUATORIAL_RADIUS, relative_tolerance=1e-9, soft=frozenset())
TransferSummary(burn_delta_v_m_s, duration_s, minimum_radius_m, thrust_n=None)
evaluate_constraints(summary, spacecraft, limits) -> ConstraintReport
```

`None` significa que el límite no se impone. Con F el empuje al que se
ejecutan los impulsos y R el radio de referencia del cuerpo central:

| Requisito | Nombre | Residuo g | Unidad | Escala | Origen del límite |
|---|---|---|---|---|---|
| RF-4.01 | `max_wet_mass` | \(m_0 - m_{max}\) | kg | \(m_{max}\) | Misión |
| RF-4.02 | `max_thrust` | \(F - \min(F_{motor}, F_{max})\) | N | Ese límite | Nave y misión |
| RF-4.03 | `max_delta_v` | \(\sum_k \lVert\Delta v_k\rVert - \Delta v_{max}\) | m/s | \(\Delta v_{max}\) | Misión |
| RF-4.04 | `max_duration` | \(t_f - t_{max}\) | s | \(t_{max}\) | Misión |
| RF-4.05 | `min_altitude` | \(h_{min} - (\min_t \lVert r(t)\rVert - R)\) | m | \(R + h_{min}\) | Misión |
| RF-4.06 | `propellant` | \(m_0\,(1 - e^{-\sum\Delta v/v_e}) - m_{prop}\) | kg | \(m_0\) | Nave |

Precisiones físicas:

- **Altitud mínima.** Es una restricción de **trayectoria**: usa el mínimo de
  ‖r‖ en toda la transferencia, no los radios en los impulsos. El resumen
  debe traer un mínimo de trayectoria. Para arcos keplerianos es exacto
  (`radius_extrema_on_arc`), no muestreado.
- **Altitud.** Se mide sobre una esfera de radio R. El cuerpo no tiene
  achatamiento ni relieve.
- **Propelente y Δv máximo son límites distintos.** El primero se deduce de
  la nave con la ecuación del cohete; el segundo es un presupuesto
  independiente. Ninguno se infiere del otro.
- **Empuje.** Necesita un nivel de empuje. El modelo impulsivo no tiene
  ninguno, así que sin él la restricción se informa como **no aplicable**,
  nunca como satisfecha. Un empuje no se infiere de un Δv.
- **Masa y empuje no dependen de la trayectoria.** Rechazan una configuración
  de nave, sea cual sea la transferencia.
- **Duración.** Es el instante del último impulso, con t = 0 en el primer
  suceso del plan. El tiempo es transcurrido, sin época.
- **Propelente.** Se evalúa siempre. Las demás, solo si su límite está fijado.

## Manejo de restricciones (`optimization.constraint_handling`)

Todos los optimizadores del proyecto solo **comparan** candidatos, mediante
una clave lexicográfica (gana la menor). Una estrategia es una regla que
convierte una evaluación —objetivo f, validez, violación total V y
violaciones v_j— en una clave.

| Estrategia | Clave | Parámetros |
|---|---|---|
| `FeasibilityRules` (por defecto) | factible (0, f) < infactible (1, V) < no válido (2, 0) | Ninguno |
| `DeathPenalty` | factible (0, f) < todo lo demás (1, 0) | Ninguno |
| `StaticPenalty(coefficient, rationale)` | válido (0, f + ρ·Σ v_j²) < no válido (1, 0) | ρ, en la unidad del objetivo |

- **`FeasibilityRules`** son las reglas de Deb (2000): un punto factible gana
  a uno infactible; dos factibles se comparan por objetivo; dos infactibles,
  por violación. No tiene parámetros, y por eso es la opción por defecto.
- **`DeathPenalty`** descarta todo punto infactible. No conserva información
  sobre la distancia a la factibilidad, así que no puede avanzar mientras la
  población no contenga ningún punto factible.
- **`StaticPenalty`** es la penalización cuadrática exterior de
  [03-fundamentos](03-fundamentos-cientificos-y-algoritmicos.md) § 15.2,
  escrita con violaciones normalizadas para que ρ tenga la unidad del
  objetivo. Es una penalización **inexacta**: para todo ρ finito, el mínimo
  de la función penalizada viola una restricción activa en una cantidad del
  orden de (sensibilidad del objetivo)/ρ. El coeficiente no tiene valor por
  defecto y debe darse con su justificación (RF-4.09).

**El objetivo nunca se modifica.** `Evaluation.objective` es siempre el
objetivo sin penalizar. La violación viaja aparte:

- `Evaluation.constraint_violation`: violación total normalizada V ≥ 0 de las
  restricciones duras; `None` si el problema no tiene restricciones;
- `Evaluation.constraints`: violación normalizada de cada una, por nombre;
- `Evaluation.is_feasible`: válida y sin violación.

**El punto que se informa se elige por factibilidad, sea cual sea la
estrategia** (`EvaluationBudget`):

| Campo de `OptimizationResult` | Contenido en un problema con restricciones |
|---|---|
| `x_best`, `u_best`, `f_best`, `best_evaluation` | El mejor punto **factible**; `None` y +inf si no hay ninguno |
| `found_feasible` | Si se encontró algún punto factible |
| `found_valid` | Si hubo alguna evaluación válida, factible o no |
| `least_violation`, `x_least_violation`, `u_least_violation`, `least_violation_evaluation` | El punto válido infactible con menor violación. Es un diagnóstico, nunca una solución |
| `constraint_handling` | Descripción de la estrategia usada |
| `history.violation` | Violación total de cada evaluación (0 si es factible, NaN si no es válida) |

Con restricciones, `history.best_so_far`, el criterio de objetivo alcanzado y
el recuento de estancamiento se refieren solo a puntos factibles.

**Compatibilidad.** Un problema sin restricciones produce exactamente el
mismo registro que antes de esta fase: los campos nuevos valen `None` y
`to_dict()` no los incluye. La referencia fijada de DE (hash SHA-256 del
2026-11-03) sigue pasando sin cambios. La selección de DE se aplica ahora
mediante claves; sin restricciones da las mismas decisiones que `_trial_wins`
(comprobado caso a caso en `tests/test_constraint_handling.py`).

La estrategia se toma, por este orden, del argumento `handling` de
`EvaluationBudget`, del atributo `constraint_handling` del problema y, si no
hay ninguno, de `FeasibilityRules`.

## Problema de transferencia con restricciones (`problems.constrained`)

```python
ConstrainedTransferProblem(problem, spacecraft, limits, *, handling=None, commanded_thrust_n=None)
```

Envuelve un problema de transferencia cuyas evaluaciones válidas publiquen
`burn_delta_v_m_s`, `transfer_time_s` y `minimum_radius_m`
(`CoplanarMultiImpulseProblem`). Un candidato queda en uno de tres estados
explícitos (RF-4.07):

| Estado | Significado | Qué se registra |
|---|---|---|
| No válido | No hay transferencia (dominio o fallo numérico) | Estado y mensaje del problema base |
| Válido e infactible | Una transferencia que incumple un límite | Objetivo sin penalizar y cada violación con su magnitud |
| Válido y factible | Cumple todos los límites duros | Objetivo y márgenes |

Metadatos añadidos a cada evaluación válida:

- `constraint_report`: el informe completo (`ConstraintReport.to_dict()`);
- `feasible`;
- `propellant_required_kg`, `propellant_per_burn_kg`, `final_mass_kg`;
- `burn_thrust_n` y `burn_durations_s` (solo con `commanded_thrust_n`).

Otros métodos:

- `report(x)` devuelve el `ConstraintReport` de un vector de decisión;
- `configuration_report()` evalúa las restricciones que no dependen de la
  trayectoria (masa y empuje). Si no es factible, ningún candidato puede serlo;
- `describe()` registra la nave, los límites y la estrategia.

**Empuje comandado.** `commanded_thrust_n` es el empuje al que la misión
piensa ejecutar los impulsos. Con él se evalúa la restricción de empuje y se
informan las duraciones ideales de los encendidos. Esas duraciones son un
diagnóstico: no se restringen y no se simula ninguna dinámica de encendido
finito (Fase 5).

`ScalarizedProblem` conserva las violaciones y la estrategia del problema que
envuelve, de modo que un objetivo ponderado y las restricciones se combinan.

## Verificación independiente (`constraints.verification`)

El optimizador evalúa la restricción de trayectoria con el mínimo exacto de
cada arco kepleriano. `numerical_transfer_summary(initial_state, impulses, mu, ...)`
obtiene las mismas magnitudes por otro camino: la integración numérica de los
impulsos tal como se aplican.

Mínimo de radio de un arco propagado (`arc_minimum_radius`):

1. se integra con RK45 y se toma el paso aceptado de menor ‖r‖ (pasos de
   60 s como máximo);
2. se integra de nuevo el intervalo entre sus dos vecinos con 201 salidas
   equiespaciadas;
3. se ajusta una parábola por la muestra menor y sus vecinas y se toma su vértice.

Con una separación h entre muestras refinadas, el vértice difiere del mínimo
real en O(‖d³r/dt³‖·h³): por debajo del milímetro en órbita baja para
h < 1 s. El error dominante es el error global del integrador (del orden de
0,1 m en una transferencia LEO–GEO con las tolerancias por defecto), no esta
búsqueda. El resultado es una **estimación numérica con esa exactitud, no una cota**.

En el benchmark, el verificador de un escenario con restricciones hace dos
comprobaciones sobre el mejor plan de cada ejecución:

1. `verify_impulsive_plan`: el plan llega a la órbita objetivo;
2. las restricciones se vuelven a evaluar sobre el resumen numérico, con una
   tolerancia relativa de 1e-6 de cada escala. El error global del integrador
   en una transferencia LEO–GEO es de unos 1e-8 del radio; 1e-6 deja dos
   órdenes de margen.

Una ejecución solo cuenta como éxito si pasa ambas.

## Solución de referencia con límite de tiempo (`problems.references`)

Con un límite de tiempo menor que el de Hohmann no hay forma cerrada del
mínimo. La referencia se calcula por un método que no comparte nada con los
optimizadores ni con la cadena de evaluación de los problemas.

**Forma cerrada** (`two_impulse_closed_form`), para una salida desde la
órbita circular r1 con primer impulso Δv a un ángulo α y r2 > r1:

    v_t = v_c1 + Δv·cos α,   v_r = Δv·sin α,   h = r1·v_t
    ε = (v_t² + v_r²)/2 − μ/r1,   a = −μ/(2ε),   e² = 1 − h²/(μ·a)
    velocidad en r2:  v2² = 2(ε + μ/r2); transversal h/r2; radial √(v2² − (h/r2)²)
    segundo impulso:  Δv2 = hypot(v_c2 − h/r2, radial)
    e·sin E1 = r1·v_r/√(μ·a),   e·cos E1 = 1 − r1/a
    cos E2 = (1 − r2/a)/e,  E2 ∈ [0, π]
    τ = ((E2 − e·sin E2) − (E1 − e·sin E1)) mod 2π / √(μ/a³)
    radio mínimo: a(1 − e) si v_r < 0 (el vehículo desciende primero); r1 si no

**Conjunto factible para una dirección α.** La transferencia alcanza r2 si
Δv ≥ Δv_alcance(α), raíz positiva de

    A·Δv² + B·Δv + C = 0,   A = r2² − r1²·cos²α,
    B = 2·v_c1·cos α·(r2² − r1²),   C = −μ·(r2 − r1)²/r1

Por encima, el tiempo decrece al crecer Δv. Los candidatos con τ ≤ t_max son
Δv ≥ Δv_min(α): el de alcance si la llegada tangente ya cumple el límite, o
la raíz de τ(Δv, α) = t_max si no. La referencia es

    min_α J(Δv_min(α), α)

**Procedimiento** (`time_limited_reference`):

1. para 4001 direcciones en (−π, π), se halla el Δv factible más barato
   (llegada tangente, o bisección de τ = t_max);
2. la mejor dirección se refina por sección áurea entre sus dos vecinas;
3. se evalúa una malla exhaustiva de 401 × 401 sobre toda la caja de decisión
   como comprobación cruzada.

Los dos supuestos (el tiempo decrece y el coste crece con Δv) **se
comprueban, no se dan por buenos**: si un nodo factible de la malla es más
barato que el mínimo de la frontera en más de 1e-6 m/s, se lanza
`RuntimeError` en lugar de devolver una referencia errónea. El punto devuelto
es factible por construcción, así que es una cota superior del mínimo real.

**Por qué no una malla que se estrecha.** El primer diseño estrechaba una
malla 2-D alrededor del mejor nodo factible. Se descartó el 2026-11-10: a lo
largo de una restricción activa, el mejor nodo puede quedar a una distancia
del óptimo que escala con √h, mientras que la ventana se estrecha con h. Ese
método dio 4 103,020 m/s para un límite de 0,8·t_H; el método actual da
4 100,057 m/s, 2,96 m/s menos.

## Benchmark con restricciones

El benchmark de la Fase 3 admite un tipo de escenario más:

| Tipo | Parámetros | Referencia |
|---|---|---|
| `constrained_transfer` | Los de `coplanar_transfer`, más `spacecraft`, `limits` y, opcional, `commanded_thrust_n` | `hohmann` (las restricciones no se activan en el óptimo), `time_limited` o `infeasible` |

- Con `reference: {"type": "infeasible"}` el escenario no tiene solución
  factible por construcción. Una ejecución acierta si informa de que no
  encontró ninguna.
- Cada entrada de `algorithms` puede llevar `constraint_handling`, para
  comparar estrategias con el mismo optimizador, problema y semillas (RF-4.10).
- El éxito exige una solución **factible**, dentro del umbral y verificada.

## Known Limitations

- **Modelo impulsivo.** La restricción de empuje solo compara un nivel de
  empuje declarado con sus límites. No hay dinámica de encendido finito ni
  pérdidas por gravedad.
- **Propelente ideal.** Tsiolkovsky con un único impulso específico, sin
  reservas, residuos ni márgenes.
- **Altitud sobre una esfera.** Sin achatamiento, relieve ni atmósfera.
- **Clase de transferencias.** Las de la Fase 3: coplanarias entre órbitas circulares.
- **Referencia con límite de tiempo.** Solo para transferencias de subida de
  dos impulsos. No hay referencia independiente para tres impulsos ni para
  otras restricciones activas.
- **Penalización estática.** Es inexacta por construcción. Su resultado
  depende del coeficiente; no hay un valor recomendado.
- **Verificación numérica.** Es una estimación con la exactitud del
  integrador, no una cota. Una solución exactamente sobre un límite puede
  diferir de él en el error del integrador; por eso su tolerancia (1e-6) es
  mayor que la analítica (1e-9).

## Referencias

- Deb, K. (2000). *An Efficient Constraint Handling Method for Genetic Algorithms*. Computer Methods in Applied Mechanics and Engineering 186.
- Coello Coello, C. A. (2002). *Theoretical and Numerical Constraint-Handling Techniques Used with Evolutionary Algorithms: A Survey of the State of the Art*. Computer Methods in Applied Mechanics and Engineering 191.
- Nocedal, J. y Wright, S. *Numerical Optimization*, cap. 17 (métodos de penalización).
- Vallado, D. A. *Fundamentals of Astrodynamics and Applications*, cap. 6 (maniobras) y § 2.2 (ecuación de Kepler).

Las referencias se citan de memoria como fuente de los métodos; no se
consultaron en línea durante esta fase. La evidencia de implementación procede
del código y los tests del repositorio.
