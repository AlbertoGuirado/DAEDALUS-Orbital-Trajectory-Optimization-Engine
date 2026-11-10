# Autonomous Orbital Transfer Optimizer — Requerimientos

> Documento 2 de 3. Complementarios: `01-caracteristicas-del-proyecto.md` y `03-fundamentos-cientificos-y-algoritmicos.md`.
>
> **Convención de identificadores:** `RF` requisito funcional, `RNF` requisito no funcional, `RV` requisito de validación, `RR` requisito de reproducibilidad, `RD` requisito de datos y configuración.
>
> **Convención de prioridad:** `OBLIGATORIO` el proyecto no es correcto sin él. `RECOMENDADO` mejora sustancialmente la calidad. `OPCIONAL` deseable, no bloqueante.

---

## 1. Requisitos transversales

Aplican a todas las fases y a todo el código del proyecto.

| ID | Prioridad | Requisito |
|---|---|---|
| RNF-01 | OBLIGATORIO | Todo cálculo interno usa unidades SI y radianes. Las conversiones a otras unidades ocurren únicamente en las fronteras de entrada y salida del sistema |
| RNF-02 | OBLIGATORIO | Toda función pública declara sus unidades explícitamente en el docstring |
| RNF-03 | OBLIGATORIO | Todo el código está anotado con tipos, incluyendo formas de array cuando sea relevante |
| RNF-04 | OBLIGATORIO | Toda implementación numérica cita en su docstring la ecuación que realiza |
| RNF-05 | OBLIGATORIO | Las capas de física, simulación y optimización mantienen dependencias unidireccionales. Un módulo de física nunca importa un optimizador |
| RNF-06 | OBLIGATORIO | Las entradas inválidas producen `ValueError`. Los fallos de cómputo producen `RuntimeError`. Nunca se devuelve silenciosamente un resultado no finito |
| RNF-07 | OBLIGATORIO | El estado de la nave se representa siempre como `[x, y, z, vx, vy, vz]`, en ese orden |
| RNF-08 | RECOMENDADO | Las estructuras de resultado son inmutables y sus arrays de solo lectura |
| RNF-09 | RECOMENDADO | Ninguna abstracción se introduce sin una necesidad presente y demostrable |
| RNF-10 | OBLIGATORIO | Ningún resultado numérico aparece en documentación o informe si no procede de una ejecución real reproducible |
| RNF-11 | OBLIGATORIO | Ningún algoritmo se declara superior a otro sin benchmark con múltiples semillas |
| RNF-12 | RECOMENDADO | Python 3.11 o superior. Dependencias limitadas a NumPy y SciPy en el núcleo |

---

## 2. Fase 1 — Orbital Dynamics Engine

Estado y evidencia: [matriz de cierre](phase1-closeout.md). Los requisitos mantienen su prioridad original; los recomendados y opcionales no implican capacidades públicas implementadas.

### 2.1 Requisitos funcionales

| ID | Prioridad | Requisito |
|---|---|---|
| RF-1.01 | OBLIGATORIO | Implementar la función derivada del problema de dos cuerpos con la firma `f(t, y, mu) -> dy/dt`, compatible con `scipy.integrate.solve_ivp` |
| RF-1.02 | OBLIGATORIO | Definir la constante `EARTH_MU` con el valor estándar del parámetro gravitacional terrestre en m³/s², citando su fuente |
| RF-1.03 | OBLIGATORIO | Implementar `validate_state` que verifique dimensión 6, finitud, tipo real y radio no nulo |
| RF-1.04 | OBLIGATORIO | Implementar `validate_mu` que verifique escalar real, finito y estrictamente positivo |
| RF-1.05 | OBLIGATORIO | Propagar un estado durante un intervalo de tiempo, devolviendo la trayectoria muestreada |
| RF-1.06 | OBLIGATORIO | Ofrecer al menos dos integradores: RK4 de paso fijo y RK45 adaptativo |
| RF-1.07 | OBLIGATORIO | Permitir especificar instantes de muestreo arbitrarios dentro del intervalo de propagación |
| RF-1.08 | OBLIGATORIO | Convertir de elementos orbitales clásicos a vector de estado (`elements_to_state`) |
| RF-1.09 | OBLIGATORIO | Convertir de vector de estado a elementos orbitales clásicos (`state_to_elements`) |
| RF-1.10 | OBLIGATORIO | Calcular la energía específica de un estado |
| RF-1.11 | OBLIGATORIO | Calcular el vector de momento angular específico de un estado |
| RF-1.12 | OBLIGATORIO | Calcular el periodo orbital a partir del semieje mayor |
| RF-1.13 | OBLIGATORIO | Calcular el Δv de un cambio de velocidad y la suma de una secuencia de impulsos |
| RF-1.14 | OBLIGATORIO | Tratar correctamente los casos degenerados de los elementos orbitales: órbita circular, órbita ecuatorial y circular-ecuatorial |
| RF-1.15 | RECOMENDADO | Calcular el vector de excentricidad y el vector de nodo como funciones expuestas |
| RF-1.16 | RECOMENDADO | Resolver la ecuación de Kepler para obtener la anomalía excéntrica a partir de la anomalía media |
| RF-1.17 | OPCIONAL | Ofrecer propagación analítica kepleriana como alternativa exacta al integrador numérico |

### 2.2 Restricciones de dominio de la Fase 1

| ID | Prioridad | Requisito |
|---|---|---|
| RF-1.18 | OBLIGATORIO | El modelo es de masa puntual sin superficie del cuerpo central. Debe documentarse que no es válido para trayectorias que atraviesen una colisión central |
| RF-1.19 | OBLIGATORIO | El dominio soportado en esta fase es el de órbitas elípticas cerradas, con excentricidad estrictamente menor que 1. Cualquier extensión a parábola o hipérbola debe declararse explícitamente |
| RF-1.20 | OBLIGATORIO | Documentar que en RK4 la especificación de instantes de muestreo altera la malla de integración, y que los estudios de convergencia deben hacerse sin ellos |

### 2.3 Requisitos de validación de la Fase 1

Cada uno de estos requisitos define un test que debe existir y pasar.

| ID | Prioridad | Test de validación |
|---|---|---|
| RV-1.01 | OBLIGATORIO | **Conservación de energía.** La energía específica se conserva a lo largo de una propagación de varias órbitas, con deriva relativa por debajo de una tolerancia declarada |
| RV-1.02 | OBLIGATORIO | **Conservación del momento angular.** El vector de momento angular específico se conserva en módulo y dirección |
| RV-1.03 | OBLIGATORIO | **Conservación del vector de excentricidad.** El vector de Laplace-Runge-Lenz normalizado se conserva en módulo y dirección |
| RV-1.04 | OBLIGATORIO | **Cierre del periodo.** Propagar exactamente un periodo orbital devuelve el estado inicial, dentro de tolerancia |
| RV-1.05 | OBLIGATORIO | **Órbita circular.** Un estado circular mantiene radio y velocidad constantes durante la propagación |
| RV-1.06 | OBLIGATORIO | **Ida y vuelta de elementos.** `state_to_elements` seguido de `elements_to_state` recupera el estado original, y la composición inversa recupera los elementos |
| RV-1.07 | OBLIGATORIO | **Relación vis-viva.** La velocidad calculada a un radio dado coincide con la predicha por la ecuación de vis-viva |
| RV-1.08 | OBLIGATORIO | **Consistencia de apsides.** El radio en periapsis y apoapsis de la trayectoria propagada coincide con los valores analíticos derivados de los elementos |
| RV-1.09 | OBLIGATORIO | **Orden de convergencia de RK4.** Al reducir el paso a la mitad, el error global se reduce en un factor próximo a 16, confirmando el cuarto orden |
| RV-1.10 | OBLIGATORIO | **Coincidencia entre integradores.** RK4 con paso suficientemente pequeño y RK45 con tolerancias estrictas producen la misma trayectoria dentro de tolerancia |
| RV-1.11 | RECOMENDADO | **Validación contra solución analítica kepleriana.** La propagación numérica coincide con la propagación kepleriana exacta |
| RV-1.12 | OBLIGATORIO | **Casos degenerados.** Los tests cubren órbita circular inclinada, elíptica ecuatorial y circular ecuatorial |
| RV-1.13 | OBLIGATORIO | **Rechazo de entradas inválidas.** Estado de dimensión incorrecta, componentes no finitas, μ nulo o negativo, duración no positiva y tolerancias no positivas producen `ValueError` |
| RV-1.14 | RECOMENDADO | **Inmutabilidad del resultado.** Los arrays devueltos no son escribibles y `final_state` devuelve una copia independiente |

### 2.4 Criterio de aceptación de la Fase 1

La Fase 1 se considera cerrada cuando se cumplen simultáneamente estas cuatro condiciones:

1. El paquete `physics` se importa sin error y todos los símbolos declarados en su `__init__.py` existen.
2. La totalidad de los tests `RV-1.*` marcados como obligatorios pasa.
3. Existe un experimento ejecutable que propaga una órbita terrestre e informa numéricamente de la deriva de los invariantes conservados.
4. Los valores de tolerancia empleados en cada test están justificados y documentados, no ajustados a posteriori para que el test pase.

---

## 3. Fase 2 — Baseline analítico

| ID | Prioridad | Requisito |
|---|---|---|
| RF-2.01 | OBLIGATORIO | Implementar la transferencia de Hohmann entre dos órbitas circulares coplanarias |
| RF-2.02 | OBLIGATORIO | Devolver el Δv de cada uno de los dos impulsos por separado, además del total |
| RF-2.03 | OBLIGATORIO | Devolver el tiempo de transferencia, igual a medio periodo de la elipse de transferencia |
| RF-2.04 | OBLIGATORIO | Calcular la masa de propelente consumida mediante la ecuación de Tsiolkovsky, dada la masa inicial y el impulso específico |
| RF-2.05 | RECOMENDADO | Implementar la transferencia bielíptica y determinar el régimen de relación de radios en que supera a Hohmann |
| RF-2.06 | RECOMENDADO | Implementar el cambio de plano puro y la maniobra combinada de cambio de plano y de altitud |
| RF-2.07 | OPCIONAL | Implementar el solucionador del problema de Lambert para transferencias de tiempo fijo |
| RV-2.01 | OBLIGATORIO | El Δv calculado por la fórmula de Hohmann coincide con el obtenido aplicando los impulsos al propagador numérico y verificando la órbita final |
| RV-2.02 | OBLIGATORIO | El tiempo de transferencia coincide con la mitad del periodo de la elipse de transferencia calculado independientemente |
| RV-2.03 | RECOMENDADO | Para el caso de radios iguales, el Δv de Hohmann tiende a cero |

**Estado (2026-11-02):** RF-2.01 a RF-2.04 y RV-2.01 a RV-2.03 implementados y verificados en `tests/test_transfers.py` y `tests/test_impulses.py`; integrados en `develop` (commit `16c5ac0`). RF-2.05, RF-2.06 y RF-2.07 **no implementados**. API y resultados: [phase2-api.md](phase2-api.md).

**Criterio de aceptación:** el baseline produce, para cualquier par de radios circulares, una solución cuya validez se verifica propagando numéricamente la trayectoria resultante y comprobando que la órbita final coincide con la objetivo dentro de tolerancia.

---

## 4. Fase 3 — Optimization Engine

### 4.1 Definición del problema de optimización

| ID | Prioridad | Requisito |
|---|---|---|
| RF-3.01 | OBLIGATORIO | Definir explícitamente el vector de variables de decisión, con su significado físico, unidades y límites inferior y superior |
| RF-3.02 | OBLIGATORIO | Definir explícitamente la función objetivo, comenzando por J igual al Δv total |
| RF-3.03 | OBLIGATORIO | Definir explícitamente el conjunto de restricciones, distinguiendo igualdades y desigualdades |
| RF-3.04 | OBLIGATORIO | Definir explícitamente el método de evaluación de una solución candidata, es decir, cómo se traduce un vector de decisión en una trayectoria propagada y en un valor de objetivo |
| RF-3.05 | OBLIGATORIO | Las variables de decisión incluyen como mínimo: instantes de los impulsos, magnitud de los impulsos y dirección de los impulsos |
| RF-3.06 | RECOMENDADO | Soportar duración de los impulsos como variable, de cara a encendidos de duración finita |
| RF-3.07 | RECOMENDADO | Soportar el número de maniobras como parámetro del problema, permitiendo comparar estrategias de dos, tres y más impulsos |
| RF-3.08 | OBLIGATORIO | Normalizar internamente las variables de decisión a un dominio común, típicamente el hipercubo unidad, para evitar que magnitudes de escalas muy distintas sesguen la búsqueda |

### 4.2 Algoritmos

| ID | Prioridad | Requisito |
|---|---|---|
| RF-3.09 | OBLIGATORIO | Implementar Differential Evolution como primer optimizador |
| RF-3.10 | OBLIGATORIO | Definir una interfaz común de optimizador que permita ejecutar cualquier algoritmo sobre exactamente el mismo problema, sin modificar el problema |
| RF-3.11 | RECOMENDADO | Implementar un algoritmo genético |
| RF-3.12 | RECOMENDADO | Implementar optimización por enjambre de partículas |
| RF-3.13 | RECOMENDADO | Implementar refinamiento local, por ejemplo Nelder-Mead o programación cuadrática secuencial |
| RF-3.14 | RECOMENDADO | Implementar un esquema híbrido que use una metaheurística para la exploración global y un método local para el refinamiento |
| RF-3.15 | OPCIONAL | Implementar control óptimo, por método directo de colocación o por método indirecto basado en el principio del máximo |
| RF-3.16 | OBLIGATORIO | Todo optimizador acepta una semilla y produce resultados idénticos para la misma semilla y configuración |
| RF-3.17 | OBLIGATORIO | Todo optimizador acepta un presupuesto máximo de evaluaciones de la función objetivo y lo respeta |
| RF-3.18 | OBLIGATORIO | Todo optimizador registra la historia de convergencia, es decir, el mejor valor encontrado en función del número de evaluaciones |

### 4.3 Función objetivo multiobjetivo

| ID | Prioridad | Requisito |
|---|---|---|
| RF-3.19 | OBLIGATORIO | Preparar la función objetivo para la forma ponderada con términos de Δv, tiempo, energía y riesgo |
| RF-3.20 | OBLIGATORIO | Todo peso y todo parámetro de la función objetivo está documentado: su significado, su unidad y la razón de su valor. No se admiten constantes arbitrarias sin justificación |
| RF-3.21 | OBLIGATORIO | Los términos de la suma ponderada se adimensionalizan antes de sumarse, mediante escalas de referencia declaradas |
| RF-3.22 | RECOMENDADO | Implementar optimización multiobjetivo genuina basada en dominancia de Pareto, sin colapsar a un escalar |
| RF-3.23 | RECOMENDADO | Calcular métricas de calidad del frente de Pareto, como el hipervolumen |

### 4.4 Validación de la Fase 3

| ID | Prioridad | Test de validación |
|---|---|---|
| RV-3.01 | OBLIGATORIO | Sobre un problema de transferencia entre órbitas circulares coplanarias con dos impulsos, el optimizador converge a una solución cuyo Δv no es significativamente mejor que el de Hohmann, que es el óptimo analítico en ese caso |
| RV-3.02 | OBLIGATORIO | Los optimizadores se validan primero sobre funciones de prueba estándar con óptimo conocido, antes de aplicarse al problema orbital |
| RV-3.03 | OBLIGATORIO | La reproducibilidad por semilla se verifica ejecutando dos veces la misma configuración y comparando resultados bit a bit o dentro de tolerancia numérica |
| RV-3.04 | OBLIGATORIO | Toda solución devuelta por un optimizador se verifica propagándola con el simulador y comprobando que alcanza efectivamente la órbita objetivo |

---

## 5. Fase 4 — Restricciones

| ID | Prioridad | Requisito |
|---|---|---|
| RF-4.01 | OBLIGATORIO | Restricción de masa máxima de la nave |
| RF-4.02 | OBLIGATORIO | Restricción de empuje máximo |
| RF-4.03 | OBLIGATORIO | Restricción de Δv máximo disponible |
| RF-4.04 | OBLIGATORIO | Restricción de tiempo máximo de misión |
| RF-4.05 | OBLIGATORIO | Restricción de altitud mínima sobre la superficie del cuerpo central, evaluada a lo largo de toda la trayectoria y no solo en los extremos |
| RF-4.06 | OBLIGATORIO | Restricción de propelente disponible |
| RF-4.07 | OBLIGATORIO | Toda solución que viole una restricción se penaliza o se descarta de forma explícita, nunca implícita |
| RF-4.08 | OBLIGATORIO | La magnitud de cada violación se registra, no solo el hecho binario de que exista |
| RF-4.09 | OBLIGATORIO | El método de manejo de restricciones está documentado, y si es de penalización, el coeficiente de penalización está justificado |
| RF-4.10 | RECOMENDADO | Comparar al menos dos estrategias de manejo de restricciones sobre el mismo problema, por ejemplo penalización frente a reglas de factibilidad |
| RV-4.01 | OBLIGATORIO | Existe un test que construye deliberadamente una solución infactible y verifica que el sistema la rechaza o penaliza según lo especificado |
| RV-4.02 | OBLIGATORIO | La restricción de altitud mínima se verifica sobre una trayectoria cuyo periapsis cae por debajo del umbral aun cuando sus extremos son válidos |

---

## 6. Fase 5 — Perturbaciones

| ID | Prioridad | Requisito |
|---|---|---|
| RF-5.01 | OBLIGATORIO | Cada modelo de perturbación puede activarse y desactivarse de forma independiente, para permitir comparación directa entre modelos |
| RF-5.02 | OBLIGATORIO | Implementar encendidos de duración finita, sustituyendo el modelo de impulso instantáneo |
| RF-5.03 | OBLIGATORIO | Implementar masa variable de la nave, acoplada al consumo de propelente durante los encendidos |
| RF-5.04 | OBLIGATORIO | Implementar la perturbación J2 del achatamiento terrestre |
| RF-5.05 | RECOMENDADO | Implementar arrastre atmosférico, con un modelo de densidad documentado y su rango de validez |
| RF-5.06 | OBLIGATORIO | La configuración de perturbaciones activas se registra junto con todo resultado experimental |
| RV-5.01 | OBLIGATORIO | Con todas las perturbaciones desactivadas, los resultados coinciden exactamente con los del modelo de dos cuerpos de la Fase 1 |
| RV-5.02 | OBLIGATORIO | La perturbación J2 reproduce la tasa de regresión nodal analítica conocida para una órbita dada |
| RV-5.03 | OBLIGATORIO | En el límite de empuje muy alto y duración muy corta, el encendido finito converge al resultado del impulso instantáneo |
| RV-5.04 | OBLIGATORIO | Con masa variable, la relación entre masa final e inicial coincide con la predicha por la ecuación de Tsiolkovsky |

---

## 7. Fase 6 — Benchmark y Mission Optimization Engine

La Fase 6 se divide en dos subfases (decisión del 2026-11-02). Los antiguos
identificadores `RF-6.0n` pasan a ser `RF-6A.0n` sin cambio de contenido ni de
prioridad: `RF-6.01` → `RF-6A.01`, …, `RF-6.08` → `RF-6A.08`.

### 7.1 Fase 6A — Benchmark y validación experimental de algoritmos

| ID | Prioridad | Requisito |
|---|---|---|
| RF-6A.01 | OBLIGATORIO | Cada experimento registra: mejor objetivo, Δv, tiempo de transferencia, tiempo de cómputo, número de evaluaciones de función, historia de convergencia, violaciones de restricción y tasa de éxito |
| RF-6A.02 | OBLIGATORIO | Cada configuración se ejecuta con múltiples semillas independientes |
| RF-6A.03 | OBLIGATORIO | Todos los algoritmos comparados reciben idéntico presupuesto de evaluaciones de función. El tiempo de pared se registra pero no se usa como criterio de parada en la comparación |
| RF-6A.04 | OBLIGATORIO | Los escenarios de benchmark se definen y fijan antes de ejecutar las comparaciones, para evitar selección a posteriori |
| RF-6A.05 | RECOMENDADO | Las comparaciones entre algoritmos se acompañan de un test estadístico no paramétrico apropiado, dado que las distribuciones de resultados no son normales |
| RF-6A.06 | OBLIGATORIO | Generar gráficas de trayectoria orbital, curvas de convergencia, Δv frente a tiempo de cómputo y comparación entre optimizadores |
| RF-6A.07 | RECOMENDADO | Generar el frente de Pareto cuando exista optimización multiobjetivo |
| RF-6A.08 | OBLIGATORIO | Toda gráfica lleva ejes etiquetados con unidades |

### 7.2 Fase 6B — Mission Optimization Engine

> **BORRADOR, pendiente de aprobación.** Requisitos propuestos para delimitar
> el alcance; no son vinculantes hasta su revisión en la puerta de la Fase 6B.

| ID | Prioridad | Requisito propuesto |
|---|---|---|
| RF-6B.01 | OBLIGATORIO | Recibir una definición de misión conforme a RD-02 a RD-05 y validarla al cargarla |
| RF-6B.02 | OBLIGATORIO | Traducir la misión a un problema de optimización usando exclusivamente las interfaces aprobadas de dinámica (Fases 1, 2 y 5), optimización (Fase 3) y restricciones (Fase 4) |
| RF-6B.03 | OBLIGATORIO | Devolver una estrategia de transferencia: secuencia de maniobras con instantes, vectores Δv, coste, tiempo, propelente y violaciones de restricción |
| RF-6B.04 | OBLIGATORIO | Verificar la estrategia devuelta propagándola con el simulador antes de informarla |
| RF-6B.05 | OBLIGATORIO | Ser determinista dada la definición de misión, la configuración y la semilla; registrar ambas con el resultado |
| RF-6B.06 | OBLIGATORIO | No presentar un óptimo local como global; informar del método, el criterio de parada y la comparación con el baseline analítico cuando exista |
| RF-6B.07 | RECOMENDADO | Comparar estrategias alternativas (número de impulsos, optimizadores) mediante la infraestructura de la Fase 6A |

---

## 8. Requisitos de reproducibilidad

| ID | Prioridad | Requisito |
|---|---|---|
| RR-01 | OBLIGATORIO | Todo componente estocástico acepta una semilla explícita. No se usa el estado global del generador de números aleatorios |
| RR-02 | OBLIGATORIO | Cada experimento se define mediante un fichero de configuración versionado, no mediante parámetros incrustados en el código |
| RR-03 | OBLIGATORIO | Cada resultado registra la configuración completa que lo produjo, incluyendo semilla, versión del código y versiones de las dependencias |
| RR-04 | OBLIGATORIO | Las versiones de las dependencias están fijadas en un fichero de instantánea |
| RR-05 | RECOMENDADO | Evolucionar el fichero de instantánea a un lock multiplataforma con hashes |
| RR-06 | OBLIGATORIO | Los experimentos se ejecutan mediante un punto de entrada estándar y documentado |
| RR-07 | OBLIGATORIO | Los artefactos de salida se escriben a un directorio de salida separado, excluido del control de versiones |
| RR-08 | RECOMENDADO | Repetir un experimento con la misma configuración produce resultados idénticos dentro de la tolerancia numérica de la plataforma |

---

## 9. Requisitos de datos y configuración

| ID | Prioridad | Requisito |
|---|---|---|
| RD-01 | OBLIGATORIO | Las constantes físicas se definen en un único lugar, con su valor, su unidad y su fuente bibliográfica |
| RD-02 | OBLIGATORIO | La definición de una misión incluye: órbita inicial, órbita objetivo, propiedades de la nave y restricciones aplicables |
| RD-03 | OBLIGATORIO | Las propiedades de la nave incluyen como mínimo: masa seca, masa de propelente, empuje máximo e impulso específico |
| RD-04 | OBLIGATORIO | Las órbitas pueden especificarse tanto por elementos orbitales como por vector de estado |
| RD-05 | OBLIGATORIO | Toda configuración se valida al cargarse, rechazando valores físicamente imposibles como excentricidad negativa o masa no positiva |
| RD-06 | RECOMENDADO | El formato de configuración es legible por humanos y diferenciable en control de versiones |

---

## 10. Requisitos de calidad de código y pruebas

| ID | Prioridad | Requisito |
|---|---|---|
| RNF-13 | OBLIGATORIO | Toda función matemática tiene al menos un test unitario |
| RNF-14 | OBLIGATORIO | Los tests de física validan contra solución analítica o invariante conservado, no contra valores previamente producidos por el propio código |
| RNF-15 | OBLIGATORIO | Las tolerancias de los tests están justificadas por análisis de error, no ajustadas hasta que el test pase |
| RNF-16 | OBLIGATORIO | Los tests cubren los casos degenerados y las entradas inválidas, no solo el camino feliz |
| RNF-17 | RECOMENDADO | La suite de tests se ejecuta en un tiempo que permita su uso frecuente durante el desarrollo |
| RNF-18 | RECOMENDADO | El código pasa un verificador estático de tipos y un analizador de estilo |
| RNF-19 | OBLIGATORIO | Existe un README que describe el proyecto, su instalación y su ejecución. `pyproject.toml` ya lo referencia como fichero de descripción larga |

---

## 11. Seguimiento de la deuda técnica inicial

Revisión de cierre del 2026-10-16 sobre `fd9d785`. Se conservan los
identificadores de la auditoría inicial para mantener trazabilidad.

| ID | Estado | Evidencia o decisión |
|---|---|---|
| DT-01 | Resuelto | `src/orbital_optimizer/dynamics/two_body.py` e imports compatibles |
| DT-02 | Resuelto | `src/orbital_optimizer/orbital/elements.py` y `src/physics/elements.py` |
| DT-03 | Resuelto | `README.md` existe y documenta instalación y ejecución |
| DT-04 | Resuelto | Baseline: 20 tests en cuatro módulos. Ampliación Kepler: 32 tests en cinco módulos, superados localmente |
| DT-05 | Resuelto | Experimento circular ejecutado con aceptación positiva |
| DT-06 | Resuelto | Layout `src/` documentado en `phase1.md` |
| DT-07 | Resuelto para Fase 1 | Solo se empaquetan capacidades implementadas y adaptadores; reservas futuras excluidas deliberadamente |
| DT-08 | Diferido | `agent/` sigue como reserva vacía; no implica un planificador implementado |
| DT-09 | Resuelto para Fase 1 | `configs/earth_orbit.json` existe; notebooks no son necesarios para ejecutar la validación |

### 11.1 Revisión de integración del 2026-11-02

Integración de `feature/phase-1-kepler` y `feature/phase-2-hohmann` en `develop`.
Estado de DT-01 a DT-09 sin cambios. Nueva deuda identificada en la auditoría
previa a la Fase 3:

| ID | Prioridad | Descripción | Fase afectada |
|---|---|---|---|
| DT-10 | ALTA | `propagate()` valida el dominio con `state_to_elements()` y lanza `ValueError` ante estados hiperbólicos, parabólicos o radiales. Un candidato del optimizador puede producirlos: la función objetivo necesita una política explícita para evaluaciones no válidas | 3 |
| DT-11 | ALTA | `propagate_with_impulses()` devuelve solo el estado final. La restricción de altitud mínima a lo largo de la trayectoria (RF-4.05, RV-4.02) exige la trayectoria o el periapsis de cada arco | 4 |
| DT-12 | MEDIA | No existe convención de época ni de escala de tiempo: el tiempo es transcurrido desde cero. Necesaria para J2/arrastre y escenarios de misión | 5, 6B |
| DT-13 | MEDIA | Los instantes de impulso deben ser estrictamente crecientes; las variables de decisión deben reparar o penalizar impulsos coincidentes | 3 |
| DT-14 | MEDIA | Coste por evaluación de RK45 con `rtol=1e-10` elevado para optimizadores poblacionales; valorar `propagate_kepler` para arcos balísticos de dos cuerpos | 3 |
| DT-15 | BAJA | La rama `feature/streamlit-ui` importa `altair` y `pandas`, presentes solo como dependencias transitivas de Streamlit; no declaradas en el extra `ui` | UI |
| DT-16 | BAJA | `build/` y `dist/` locales contienen un paquete de la Fase 1 sin la Fase 2 (ignorados por Git); regenerar o eliminar antes de cualquier publicación | Release |
| DT-08 | Diferido | `agent/` sigue vacío; su destino se decidirá en la puerta de la Fase 6B, como capa determinista, no de IA | 6B |

La [matriz de cierre](phase1-closeout.md) concreta la cobertura real de
RF-1 y RV-1, incluyendo el alcance del acuerdo entre integradores y los
requisitos recomendados/opcionales que no forman parte de la API.

---

## 12. Matriz de trazabilidad por fase

| Fase | Requisitos funcionales | Requisitos de validación | Criterio de cierre |
|---|---|---|---|
| 1 | RF-1.01 a RF-1.20 | RV-1.01 a RV-1.14 | Paquete importable, tests obligatorios en verde, experimento de validación ejecutable con deriva de invariantes informada |
| 2 | RF-2.01 a RF-2.07 | RV-2.01 a RV-2.03 | Baseline verificado por propagación numérica de la trayectoria resultante |
| 3 | RF-3.01 a RF-3.23 | RV-3.01 a RV-3.04 | Differential Evolution reproduce el óptimo de Hohmann en el caso donde es el óptimo conocido |
| 4 | RF-4.01 a RF-4.10 | RV-4.01 a RV-4.02 | Toda violación de restricción detectada, cuantificada y tratada explícitamente |
| 5 | RF-5.01 a RF-5.06 | RV-5.01 a RV-5.04 | Con perturbaciones desactivadas se recupera exactamente el modelo de dos cuerpos |
| 6A | RF-6A.01 a RF-6A.08 | — | Benchmark completo con múltiples semillas y escenarios fijados de antemano |
| 6B | RF-6B.01 a RF-6B.07 (borrador) | Por definir en su puerta de fase | Misión → estrategia verificada por propagación, reproducible y comparada con el baseline |

---

## 13. Definición de "hecho"

Un componente se considera terminado cuando cumple las siete condiciones siguientes:

1. Está implementado con tipado completo y docstrings que declaran unidades.
2. Cita la ecuación que implementa y su fuente.
3. Tiene tests unitarios que validan contra solución analítica o invariante conservado.
4. Sus casos degenerados y sus entradas inválidas están cubiertos por test.
5. Sus tolerancias están justificadas por análisis, no ajustadas a posteriori.
6. No rompe ningún test existente.
7. Su decisión de diseño está documentada cuando existían alternativas razonables.

## Estado de implementación tras Phase 2 — 2026-11-02

- **Implemented:** RF-1 y su ampliación Kepler según el contrato de Phase 1;
  RF-2.01–04 y RV-2.01–03 con los tests de transferencias e impulsos.
- **Planned:** RF-2.05–07 y requisitos de fases posteriores. La condición
  OBLIGATORIO define el criterio de su fase, no una capacidad ya disponible.
- **Known Limitations:** DT-10, DT-11 y DT-13 se observan en las APIs actuales.
  DT-14 es una preocupación de rendimiento pendiente de benchmark, no una
  medición demostrada por esta revisión. DT-15 y DT-16 no se revalidaron aquí.

La [auditoría de integración](phase1-phase2-audit.md) registra 47 tests OK,
reproducción de casos y límites. No se han modificado requisitos ni tolerancias.

## Estado de implementación tras Phase 3 — 2026-11-10

Rama `feature/phase-3-optimization`. Los requisitos no cambian; esta sección
solo registra su estado. Trazabilidad requisito a requisito y evidencia:
[phase3-validation.md](phase3-validation.md).

- **Implemented:** RF-3.01 a RF-3.05, RF-3.07 a RF-3.14 y RF-3.16 a RF-3.21;
  RV-3.01 a RV-3.04. También RF-2.05 (transferencia bielíptica) y RD-01
  (`constants.py` como único lugar de las constantes).
- **No implementado:** RF-3.06 (recomendado), RF-3.15 (opcional), RF-3.22 y
  RF-3.23 (recomendados); RF-2.06 y RF-2.07.
- **Matiz de RF-3.03:** la igualdad "llegar a la órbita objetivo" se elimina
  por construcción (el último impulso se deriva) y los predicados de dominio
  actúan como desigualdades. Las restricciones de misión corresponden a la Fase 4.
- **Matiz de RF-3.19:** los términos de energía y riesgo son indicadores
  declarados y documentados, no modelos validados; su peso por defecto es cero.

### Deuda técnica: revisión del 2026-11-10

| ID | Estado | Evidencia o decisión |
|---|---|---|
| DT-10 | Resuelto | Estados `INVALID_DOMAIN` y `NUMERICAL_FAILURE`; predicados analíticos previos a la física |
| DT-11 | Resuelto | `propagate_impulsive_arcs` conserva los arcos; `radius_extrema_on_arc` da el mínimo exacto de cada arco |
| DT-13 | Resuelto | Los impulsos coincidentes se rechazan de forma explícita como `INVALID_DOMAIN` |
| DT-14 | Medido | El objetivo usa propagación analítica; coste por evaluación en [phase3-validation.md](phase3-validation.md) |
| DT-12, DT-15, DT-16, DT-08 | Sin cambios | Corresponden a fases posteriores, a la interfaz o a la publicación |
