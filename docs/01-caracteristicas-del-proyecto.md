# Autonomous Orbital Transfer Optimizer — Características del Proyecto

> Documento 1 de 3. Complementarios: `02-requerimientos.md` y `03-fundamentos-cientificos-y-algoritmicos.md`.
>
> **Estado del documento (2026-10-16):** alcance y arquitectura revisados sobre `fd9d785` (`v0.1.0`). Fase 1 implementada y validada para los casos documentados. Ver [validación](phase1-validation.md) y [cierre](phase1-closeout.md). Las fases posteriores describen objetivos futuros.

---

## 1. Identidad del proyecto

**Nombre:** DAEDALUS — Orbital Trajectory Optimization Engine (nombre inicial: Autonomous Orbital Transfer Optimizer).

**Definición en una frase:** motor de optimización algorítmica que, dada una órbita inicial, una órbita objetivo y un conjunto de restricciones de misión, calcula automáticamente la secuencia de maniobras que minimiza el coste de la transferencia.

**Función objetivo primaria:** minimizar el Δv total (equivalentemente, el propelente consumido vía la ecuación de Tsiolkovsky), considerando además tiempo de transferencia, energía y restricciones operativas de la nave.

**Naturaleza del proyecto:** ingeniería científica. El núcleo es matemático y algorítmico: dinámica orbital, simulación numérica, optimización global y local, control óptimo y benchmarking estadístico.

---

## 2. Alcance

### 2.1 Alcance objetivo del proyecto

La entrega actual cubre dos cuerpos, elementos orbitales, invariantes y costes Δv. Las extensiones, optimizadores y restricciones de esta tabla son trabajo futuro.

| Dimensión | Contenido |
|---|---|
| Dinámica orbital | Propagación numérica del problema de dos cuerpos, con extensiones perturbativas activables |
| Optimización | Metaheurísticas poblacionales, refinamiento local y control óptimo sobre un problema común |
| Multiobjetivo | Estudio del compromiso Δv / tiempo / energía / riesgo, con frentes de Pareto |
| Restricciones | Masa, empuje máximo, Δv máximo, tiempo de misión, altitud mínima, propelente disponible |
| Benchmark | Comparación reproducible de algoritmos sobre escenarios idénticos, con múltiples semillas |
| Producto final | Mission Optimization Engine: entra un requisito de misión, sale una estrategia de transferencia |

### 2.2 Lo que el proyecto NO ES

- **No es un proyecto de IA o de modelos de lenguaje.** No hay aprendizaje automático en el núcleo. Si en algún momento se incorpora, será como componente auxiliar claramente aislado y justificado por benchmark.
- **No es un visualizador de órbitas.** La visualización existe para inspeccionar y comunicar resultados, nunca como objetivo.
- **No es un sustituto de herramientas de vuelo operacional** (GMAT, STK, Orekit). No se pretende precisión de navegación operativa ni modelos de efemérides de alta fidelidad.
- **No es una interfaz gráfica.** Se prioriza código científico reproducible sobre cualquier front-end.

### 2.3 Fronteras explícitas (fuera de alcance en la versión actual)

- Problema de N cuerpos completo y dinámica de tres cuerpos restringida (CR3BP).
- Efemérides planetarias reales (SPICE/JPL DE).
- Transferencias interplanetarias con asistencias gravitatorias.
- Modelos de propulsión detallados (cámaras, ciclos, transitorios de encendido).
- Navegación, determinación de órbita y filtrado estocástico.
- Control en lazo cerrado y guiado en tiempo real.

---

## 3. Arquitectura

### 3.1 Arquitectura propuesta en la especificación original

```
orbital-optimizer/
├── src/
│   ├── dynamics/
│   ├── orbital/
│   ├── optimization/
│   ├── missions/
│   ├── simulation/
│   └── visualization/
├── tests/
├── experiments/
├── notebooks/
├── configs/
├── docs/
└── README.md
```

### 3.2 Arquitectura vigente

```text
src/orbital_optimizer/
  dynamics/two_body.py
  orbital/mechanics.py
  orbital/elements.py
  simulation/propagator.py
  visualization/plots.py
src/physics/                 adaptadores de compatibilidad
experiments/validate_two_body.py
configs/earth_orbit.json
  orbital/kepler.py          solver elíptico (ampliación local)
  simulation/kepler_propagator.py  propagación analítica (ampliación local)
tests/                      32 tests en cinco módulos, ampliación validada localmente
docs/                       documentación técnica
```

### 3.3 Resultado de la auditoría inicial

Los bloqueos históricos E1–E6 (dinámica y elementos ausentes, imports rotos,
README, tests y experimento vacíos) están resueltos en `v0.1.0`.
Se adoptó layout `src`, con física, geometría orbital y simulación separadas.
El empaquetado publica `orbital_optimizer*` y `physics*`; las reservas vacías
`optimization/`, `missions/` y `agent/` permanecen fuera del paquete instalado.
No son funcionalidades disponibles. El inventario y la evidencia del cierre
se encuentran en [phase1-closeout.md](phase1-closeout.md).

---

## 4. Propagación numérica implementada

`src/orbital_optimizer/simulation/propagator.py` implementa la propagación numérica. Características observadas:

**Dos integradores seleccionables**

- **RK45**: Dormand-Prince 5(4) adaptativo, delegado en `scipy.integrate.solve_ivp`. Control de error local por tolerancias.
- **RK4**: Runge-Kutta clásico de cuarto orden, implementado a mano, con paso fijo igual a `max_step`, acortado para alcanzar exactamente el instante final y cada instante de muestreo solicitado.

**Tolerancias por defecto**

| Parámetro | Valor por defecto | Unidad |
|---|---|---|
| Tolerancia relativa | 1e-10 | adimensional |
| Tolerancia absoluta, posición | 1e-3 | m |
| Tolerancia absoluta, velocidad | 1e-6 | m/s |
| Paso máximo | 60 | s |

**Diseño defensivo.** El código valida exhaustivamente las entradas y distingue dos familias de error: `ValueError` para entradas inválidas y `RuntimeError` para fallos de integración. Rechaza tolerancias relativas por debajo de 100 veces el épsilon de máquina en lugar de dejar que SciPy las recorte en silencio. Envuelve las operaciones aritméticas en `np.errstate` con `over`, `invalid` y `divide` en modo `raise`, de modo que un desbordamiento no se propaga como `NaN`.

**Resultado inmutable.** `PropagationResult` es un `dataclass` congelado que copia los arrays de entrada y los marca de solo lectura. Verifica que los instantes sean finitos y estrictamente crecientes y que la matriz de estados tenga forma `(N, 6)` y sea finita. La propiedad `final_state` devuelve una copia independiente.

**Limitaciones documentadas en el propio código**

- El modelo es de masa puntual autónoma, sin superficie del cuerpo central. No debe usarse para continuar trayectorias que atraviesen una colisión central.
- En RK4, pedir instantes de muestreo altera la malla de integración, porque el paso se acorta para aterrizar exactamente en cada muestra. Para estudios de convergencia con paso prescrito hay que usar `sample_times=None`.
- `rtol` y `atol` se validan siempre pero solo los usa RK45. RK4 no tiene control de error.
- `nfev` cuenta evaluaciones de la derivada. No es el número de muestras de salida ni una medida de tiempo de pared.

---

## 5. Stack tecnológico y convenciones

### 5.1 Entorno

| Elemento | Valor |
|---|---|
| Lenguaje | Python |
| Versión mínima declarada | 3.11 |
| Versión usada en la validación | 3.12.1 |
| Dependencias de ejecución | NumPy (>=2.2,<3), SciPy (>=1.16,<2) |
| Versiones congeladas | NumPy 2.3.3, SciPy 1.16.3 |
| Backend de empaquetado | setuptools >= 68 |

El fichero `requirements-reproducible.txt` es una instantánea de versiones, no un lock multiplataforma con hashes. Esta distinción está documentada en el propio fichero y debe mantenerse honesta.

### 5.2 Convención de unidades

**Regla absoluta: todo el interior del sistema trabaja en unidades SI y radianes.** Las conversiones a kilómetros, grados, días o unidades astronómicas ocurren exclusivamente en las fronteras de entrada y salida.

| Magnitud | Unidad interna | Símbolo |
|---|---|---|
| Longitud | metro | m |
| Tiempo | segundo | s |
| Masa | kilogramo | kg |
| Velocidad | metro por segundo | m/s |
| Aceleración | metro por segundo al cuadrado | m/s² |
| Ángulo | radián | rad |
| Parámetro gravitacional μ | metro cúbico por segundo al cuadrado | m³/s² |
| Energía específica | julio por kilogramo | J/kg = m²/s² |
| Momento angular específico | metro cuadrado por segundo | m²/s |
| Empuje | newton | N |
| Impulso específico | segundo | s |

La razón de esta rigidez es evitar la clase de error que ha destruido misiones reales: la pérdida del Mars Climate Orbiter en 1999 se debió a una discrepancia entre unidades imperiales y métricas en el software de navegación.

### 5.3 Convención de representación del estado

El estado de la nave es un vector de seis componentes en un marco inercial:

```
y = [x, y, z, vx, vy, vz]
```

Las tres primeras componentes en metros, las tres últimas en metros por segundo. Este orden es el que asume `propagator.py` y no debe alterarse.

### 5.4 Convenciones de código

- **Tipado estático** mediante anotaciones, incluyendo `numpy.typing.NDArray` y `ArrayLike` para arrays.
- **Docstrings** que declaran unidades explícitamente en cada función pública.
- **Argumentos solo por nombre** (`*`) para las opciones de configuración, forzando llamadas legibles.
- **Estructuras inmutables** para los resultados, mediante `dataclass(frozen=True)` y arrays de solo lectura.
- **Separación estricta** entre física, simulación y optimización. Un módulo de física no importa nunca un optimizador.
- **Sin abstracciones innecesarias.** No se introduce una jerarquía de clases donde una función basta.

---

## 6. Modelo de desarrollo por fases

El proyecto se construye de forma incremental. Cada fase debe estar validada antes de comenzar la siguiente. Esta es una decisión metodológica deliberada: un error en la capa de dinámica contamina silenciosamente cualquier resultado de optimización construido encima.

### Fase 1 — Orbital Dynamics Engine

Propagador de dos cuerpos fiable. Conversión entre elementos orbitales y vectores de estado. Cálculo de energía específica, momento angular y Δv. Tests contra soluciones analíticas conocidas.

**Estado actual: implementada y validada para los casos de Fase 1.** El baseline pasó 20 tests; la ampliación local Kepler pasa 32. El experimento circular cumple sus cinco límites. Ver [validación](phase1-validation.md) y [cierre](phase1-closeout.md).

### Fase 2 — Baseline analítico

Transferencia de Hohmann entre órbitas circulares coplanarias. Para cada misión: Δv total, Δv de cada impulso, tiempo de transferencia y masa de propelente. Sirve como referencia contra la que medir a los optimizadores.

**Estado: no iniciada.**

### Fase 3 — Optimization Engine

Representación de una trayectoria como vector de variables de decisión. Definición explícita de variables, función objetivo, restricciones y método de evaluación. Implementación de Differential Evolution. Arquitectura preparada para añadir algoritmos genéticos, enjambre de partículas, optimización local y control óptimo, todos comparables sobre el mismo problema.

**Estado: no iniciada.**

### Fase 4 — Constraints

Restricciones realistas añadidas progresivamente: masa máxima, empuje máximo, Δv máximo, tiempo máximo de transferencia, altitud mínima y propelente disponible. Toda violación debe penalizarse o descartarse de forma explícita y registrada.

**Estado: no iniciada.**

### Fase 5 — Perturbations

Extensiones al modelo, cada una activable y desactivable de forma independiente para permitir comparación entre modelos: encendidos de duración finita, masa variable, arrastre atmosférico y perturbación J2 terrestre.

**Estado: no iniciada.**

### Fase 6 — Mission Optimization Engine

Capa final que recibe requisitos de misión y genera automáticamente una estrategia de transferencia óptima, apoyándose en todo lo anterior.

**Estado: no iniciada.**

---

## 7. Sistema de experimentos

La Fase 1 dispone de un experimento determinista que registra configuración, versiones, estados, errores e invariantes. Esta tabla corresponde al sistema futuro de comparación de optimizadores:

| Métrica | Descripción |
|---|---|
| Mejor objetivo | Valor mínimo de J alcanzado |
| Δv | Desglosado por impulso y total, en m/s |
| Tiempo de transferencia | En segundos |
| Tiempo de cómputo | Tiempo de pared, en segundos |
| Evaluaciones de función | Número de llamadas a la función objetivo |
| Convergencia | Historia del mejor valor frente a evaluaciones |
| Violaciones de restricción | Magnitud y tipo de cada violación |
| Tasa de éxito | Fracción de ejecuciones que alcanzan una solución factible |

Las futuras ejecuciones estocásticas se repetirán con múltiples semillas; el experimento actual no utiliza aleatoriedad. Cada experimento queda definido por un fichero de configuración versionado.

**Gráficas previstas:** trayectoria orbital, curvas de convergencia, Δv frente a tiempo de cómputo, comparación entre optimizadores y frente de Pareto cuando exista optimización multiobjetivo.

---

## 8. Principios de ingeniería vinculantes

Estos principios gobiernan todas las decisiones técnicas del proyecto y no son negociables sin justificación explícita.

1. **No inventar resultados.** Ningún número aparece en un informe si no procede de una ejecución real y reproducible.
2. **No afirmar superioridad sin benchmark.** Ningún algoritmo se declara mejor que otro antes de haberlo medido con múltiples semillas y escenarios.
3. **Validar antes de construir encima.** Cada componente matemático se valida contra solución analítica o invariante conservada antes de que otro dependa de él.
4. **Separar física, simulación y optimización.** Las tres capas tienen fronteras nítidas y dependencias unidireccionales.
5. **Evitar abstracciones innecesarias.** La complejidad estructural debe estar justificada por una necesidad real y presente.
6. **Documentar las ecuaciones.** Cada implementación numérica cita la ecuación que realiza y su fuente.
7. **Tests unitarios para toda la matemática.** Ninguna función de física se da por buena sin test.
8. **Reproducibilidad por semilla y configuración.** Cualquier resultado debe poder regenerarse exactamente.
9. **Unidades explícitas y consistentes.** SI en el interior, conversión solo en las fronteras.
10. **Explicar las alternativas de diseño.** Cuando una decisión tenga varias opciones razonables, se documentan ventajas e inconvenientes antes de implementar.

---

## 9. Riesgos identificados

| Riesgo | Descripción | Mitigación |
|---|---|---|
| R1 | Deriva hacia un proyecto de IA a través del paquete `agent/` | Mantener el motor determinista y autosuficiente. Cualquier componente de IA debe ser opcional y estar aislado |
| R2 | Error silencioso en la capa de dinámica que contamina todos los resultados | Validación contra invariantes conservados y soluciones analíticas antes de la Fase 2 |
| R3 | Deriva energética por integración no simpléctica en propagaciones largas | Medir la deriva de energía explícitamente. Documentar el régimen de validez de cada integrador |
| R4 | Comparación injusta entre optimizadores por presupuestos de evaluación distintos | Fijar presupuesto idéntico de evaluaciones de función para todos los algoritmos comparados |
| R5 | Conclusiones estadísticamente infundadas por número insuficiente de semillas | Fijar un mínimo de semillas por configuración y aplicar tests no paramétricos |
| R6 | Sobreajuste del optimizador a un único escenario | Conjunto de escenarios de benchmark diverso y fijado de antemano |
| R7 | Irreproducibilidad por dependencias no fijadas | Instantánea de versiones y, a futuro, lock con hashes |

---

## 10. Estado global resumido

| Fase | Componente | Estado revisado el 2026-10-16 |
|---|---|---|
| 1 | Dos cuerpos, invariantes, elementos y RK4/RK45 | Implementados y probados |
| 1 | Tests del baseline | 20 superados localmente |
| 1 | Ampliación Kepler local | 32 tests superados; modelo de error de redondeo documentado |
| 1 | Experimento circular y gráficos | Ejecutado; cinco criterios cumplidos |
| 1 | API, contratos y cierre documental | Documentados en `docs/` y enlazados al vault |
| 2 | Hohmann e impulsos aplicados | Pendientes de ingeniería |
| 3–6 | Optimización, restricciones, perturbaciones y misiones | No implementadas |

La evidencia es local y se limita a los casos y al entorno registrados.
No se ha comprobado el estado remoto de GitHub Actions en esta revisión.
