# Mapa de carpetas y ficheros

Guía breve del árbol de `orbital-optimizer`, revisada el **17 de octubre de
2026**. Cada entrada explica su papel actual. Incluye los cambios locales
de Kepler; no todos pertenecen al tag `v0.1.0`.

Para orientarte: `src/` calcula, `configs/` configura, `experiments/` ejecuta,
`tests/` comprueba, `outputs/` guarda resultados y `docs/` explica el proyecto.
Los archivos generados se describen por su función, sin enumerar cada copia
de paquetes, cachés o resultados de ejecuciones anteriores.

## Raíz del repositorio

**`README.md`.** Presenta DAEDALUS y explica cómo instalarlo, ejecutar los
tests y lanzar el experimento de referencia.

**`CONTRIBUTING.md`.** Define cómo trabajar en el repositorio: ramas,
verificación local, organización científica y criterios de colaboración.

**`CHANGELOG.md`.** Resume los cambios publicados por versión. Su entrada
0.1.0 describe el baseline original; no es un inventario de cambios locales.

**`pyproject.toml`.** Declara el paquete Python, su versión, dependencias y
reglas de construcción. Publica los paquetes de `src/` que corresponden al núcleo.

**`requirements-reproducible.txt`.** Fija las versiones de NumPy y SciPy
utilizadas en la validación de referencia para facilitar su reproducción.

**`.gitignore`.** Indica qué archivos locales, resultados regenerables y
datos sensibles deben quedar fuera del seguimiento de Git.

**`.gitattributes`.** Establece cómo trata Git los archivos de texto,
los finales de línea y las imágenes binarias.

## Código científico: `src/`

**`src/`.** Contiene el código que se instala como paquete Python. Es la
fuente de implementación del motor científico.

**`src/orbital_optimizer/`.** Es el paquete principal: reúne dinámica,
geometría orbital, propagación y gráficos científicos.

**`src/orbital_optimizer/__init__.py`.** Expone la API principal para poder
importar funciones como `propagate` o `solve_kepler` desde un único lugar.

**`src/orbital_optimizer/dynamics/`.** Agrupa el modelo de aceleración y las
validaciones básicas, independientemente de los métodos de propagación.

**`src/orbital_optimizer/dynamics/__init__.py`.** Identifica el subpaquete
de dinámica y documenta su responsabilidad; no añade cálculos.

**`src/orbital_optimizer/dynamics/two_body.py`.** Define la gravedad de dos
cuerpos, la constante terrestre y la validación del estado y de μ.

**`src/orbital_optimizer/orbital/`.** Agrupa las representaciones y relaciones
matemáticas de las órbitas, sin encargarse de integrar sus ecuaciones de movimiento.

**`src/orbital_optimizer/orbital/__init__.py`.** Identifica el subpaquete
de geometría orbital y deja constancia de sus convenciones de unidades.

**`src/orbital_optimizer/orbital/elements.py`.** Define los elementos
orbitales clásicos y su conversión a posición/velocidad y viceversa,
incluidas las convenciones para órbitas circulares y ecuatoriales.

**`src/orbital_optimizer/orbital/mechanics.py`.** Calcula energía, momento
angular, excentricidad, nodo, periodo y costes Δv. Calcula el coste de los
impulsos, pero no los aplica a una trayectoria.

**`src/orbital_optimizer/orbital/kepler.py`.** Resuelve la ecuación elíptica
de Kepler mediante Newton con salvaguarda por bisección, relacionando
anomalía media y excéntrica.

**`src/orbital_optimizer/simulation/`.** Contiene los métodos que hacen
evolucionar un estado orbital durante un intervalo de tiempo.

**`src/orbital_optimizer/simulation/__init__.py`.** Identifica el subpaquete
de simulación; las implementaciones se encuentran en sus módulos.

**`src/orbital_optimizer/simulation/propagator.py`.** Propaga órbitas con
RK4 o RK45 y define `PropagationResult`, que contiene tiempos, estados y
datos del método utilizado.

**`src/orbital_optimizer/simulation/kepler_propagator.py`.** Propaga
analíticamente el modelo de dos cuerpos mediante coeficientes de Lagrange
f/g, con resolución numérica de la raíz de Kepler y aritmética float64.

**`src/orbital_optimizer/visualization/`.** Contiene los gráficos científicos
que muestran las trayectorias calculadas y sus errores de validación.

**`src/orbital_optimizer/visualization/__init__.py`.** Identifica el
subpaquete gráfico, manteniendo Matplotlib como dependencia opcional.

**`src/orbital_optimizer/visualization/plots.py`.** Exporta las figuras
`trajectory.png` y `consistency.png` a partir de los resultados calculados;
no implementa otra física ni una interfaz interactiva de misión.

## Compatibilidad: `src/physics/`

**`src/physics/`.** Conserva imports antiguos para usar el núcleo actual
desde el nombre de paquete `physics`, sin duplicar sus algoritmos.

**`src/physics/__init__.py`.** Reexporta los símbolos de compatibilidad
desde `orbital_optimizer`; no expone necesariamente toda su API nueva.

**`src/physics/dynamics.py`.** Reexporta la dinámica, las validaciones y
algunas funciones de invariantes y costes desde sus módulos actuales.

**`src/physics/elements.py`.** Reexporta los elementos orbitales, sus
conversiones y el cálculo del periodo.

**`src/physics/propagator.py`.** Reexporta `propagate` y `PropagationResult`
para mantener compatibles los imports del propagador numérico.

## Configuración y experimentos

**`configs/`.** Guarda los parámetros de escenarios reproducibles separados
del código que los ejecuta.

**`configs/earth_orbit.json`.** Define la órbita circular terrestre del
experimento, el integrador, las muestras, las tolerancias y los cinco límites
de aceptación.

**`experiments/`.** Contiene los puntos de ejecución que combinan
configuración, motor científico y exportación de resultados.

**`experiments/__init__.py`.** Identifica el paquete de experimentos para
poder ejecutar sus módulos con `python -m experiments...`.

**`experiments/validate_two_body.py`.** Ejecuta el caso circular, lo compara
con una referencia analítica y guarda métricas y trayectoria. Devuelve un
código de salida que indica aceptación, incumplimiento o error.

## Pruebas: `tests/`

**`tests/`.** Reúne las comprobaciones automáticas del motor, sus contratos
de entrada/salida y las referencias independientes de validación.

**`tests/__init__.py`.** Identifica el paquete de tests y permite imports
entre sus módulos, como el uso de referencias analíticas compartidas.

**`tests/test_dynamics.py`.** Comprueba gravedad, invariantes circulares,
periodo, costes Δv, validaciones de entrada y ciertos desbordamientos numéricos.

**`tests/test_elements.py`.** Comprueba las conversiones orbitales,
los apsides, vis-viva, las geometrías degeneradas y el rechazo de estados
fuera del dominio soportado.

**`tests/test_propagator.py`.** Contrasta RK4/RK45 con referencias
analíticas y entre sí, y comprueba conservación, convergencia, muestreo,
apsides y fallos del integrador.

**`tests/test_kepler.py`.** Comprueba el solver de Kepler y la propagación
analítica mediante bisección independiente, geometrías variadas, invariantes
y cierre de múltiples vueltas.

**`tests/test_experiment.py`.** Comprueba que el experimento acepta el caso
de referencia y rechaza límites incumplidos o configuraciones inválidas.

## Documentación: `docs/`

**`docs/`.** Contiene la documentación técnica del proyecto: alcance,
requisitos, ecuaciones, interfaces y evidencia de validación.

**`docs/README.md`.** Es el índice de lectura de la documentación y enlaza
la base de conocimiento científica complementaria.

**`docs/mapa-del-proyecto.md`.** Es esta guía: explica brevemente para qué
sirve cada carpeta y fichero del árbol de trabajo relevante.

**`docs/01-caracteristicas-del-proyecto.md`.** Describe la identidad,
arquitectura, alcance actual y objetivos futuros de DAEDALUS.

**`docs/02-requerimientos.md`.** Define requisitos identificados, sus
prioridades y criterios de aceptación, además del seguimiento de la deuda inicial.

**`docs/03-fundamentos-cientificos-y-algoritmicos.md`.** Desarrolla las bases
teóricas de mecánica orbital y optimización para las distintas fases; parte
de ese contenido corresponde a capacidades futuras.

**`docs/phase1.md`.** Documenta las ecuaciones, unidades, convenciones,
métodos y criterios de validación de la primera fase.

**`docs/phase1-api.md`.** Explica cómo usar las interfaces actuales y cuál
es el contrato de sus estados, parámetros, resultados y archivos exportados.

**`docs/phase1-validation.md`.** Registra las ejecuciones observadas,
entornos, métricas y resultados de tests, incluido el historial del criterio
de error de la ampliación Kepler.

**`docs/phase1-closeout.md`.** Relaciona los requisitos con implementación
y pruebas, y recoge la auditoría documental y el alcance del cierre local.

**`docs/informe-visual.html`.** Es un documento HTML local que presenta
un mapa interactivo del sistema; es material de presentación y no ejecuta
el motor de dinámica orbital.

## Automatización y control de versiones

**`.github/`.** Guarda la configuración de colaboración y automatización
utilizada por GitHub.

**`.github/pull_request_template.md`.** Propone una plantilla para explicar
el problema, la validación y los límites de cada cambio enviado a revisión.

**`.github/workflows/`.** Contiene los flujos de ejecución automática de
GitHub Actions.

**`.github/workflows/tests.yml`.** Configura tests en varias versiones de
Python y sistemas operativos, y un experimento que publica artefactos de
validación. La configuración no demuestra por sí sola que una ejecución haya pasado.

**`.git/`.** Almacena historial, ramas y metadatos del repositorio. Su
contenido lo gestiona Git y no forma parte del motor científico.

## Reservas de fases futuras

**`optimization/`.** Reserva módulos para futuros optimizadores. Sus
cuatro archivos siguientes están vacíos; todavía no proporcionan algoritmos.

**`optimization/differential_evolution.py`.** Reserva el lugar histórico
previsto para evolución diferencial; su implementación está pendiente.

**`optimization/genetic_algorithm.py`.** Reserva un posible algoritmo
genético para etapas posteriores.

**`optimization/hybrid_optimizer.py`.** Reserva un posible método que
combine búsqueda global y refinamiento local.

**`optimization/optimal_control.py`.** Reserva el desarrollo futuro de
métodos de control óptimo.

**`missions/`.** Reserva la definición y coordinación de misiones futuras;
los dos archivos siguientes están vacíos.

**`missions/mission_definition.py`.** Reserva la futura representación de
los datos y requisitos de una misión.

**`missions/orbit_transfer.py`.** Reserva la futura coordinación de
transferencias entre órbitas.

**`agent/`.** Contiene una reserva histórica para planificación; no
implementa actualmente un agente ni un planificador.

**`agent/mission_planner.py`.** Es un archivo vacío reservado para una
posible capa futura de planificación de misiones.

**`benchmarks/`.** Es una carpeta actualmente vacía, disponible para
comparaciones de rendimiento futuras; no contiene benchmarks ejecutables.

**`physics/` en la raíz.** Es una carpeta local vacía. La compatibilidad
activa está en `src/physics/`.

## Material local de visualización y conocimiento

**`visualization/` en la raíz.** Guarda ideas e imágenes de referencia para
una posible interfaz. Es distinta del módulo científico
`src/orbital_optimizer/visualization/`.

**`visualization/IDEA.md`.** Describe una propuesta de entorno de diseño
y análisis de misiones; sus pantallas y funciones son ideas, no una aplicación
implementada en esta carpeta.

**`visualization/C.png`.** Es una imagen local de referencia visual;
el motor y el experimento no la utilizan.

**`visualization/S.png`.** Es otra imagen local de referencia visual,
separada de los gráficos producidos por la simulación.

**`visualization/image.png`.** Es material gráfico local de apoyo al diseño;
no forma parte de los cálculos científicos.

**`visualization/world map1.jpg`.** Es un recurso gráfico local de mapa
para exploración visual; no interviene en el modelo orbital.

**`visualization/world map2.jpg`.** Es un segundo recurso gráfico local
de mapa, también independiente de la propagación actual.

**`vault/`.** Contiene únicamente configuración local de Obsidian. La base
de conocimiento con notas científicas está fuera del repositorio, en la
carpeta hermana `../DAEDALUS-Vault/`.

**`vault/.obsidian/`.** Guarda preferencias del espacio local de Obsidian;
no contiene algoritmos ni evidencia científica.

**`vault/.obsidian/app.json`.** Guarda opciones generales de Obsidian para
ese espacio local.

**`vault/.obsidian/appearance.json`.** Guarda preferencias de apariencia
del espacio de Obsidian.

**`vault/.obsidian/core-plugins.json`.** Guarda la selección de plugins
integrados de Obsidian para ese espacio.

**`vault/.obsidian/workspace.json`.** Conserva la disposición local de
paneles y archivos abiertos en Obsidian.

## Entorno y archivos generados

**`.venv/`.** Contiene el intérprete y las dependencias del entorno Python
local; no es código desarrollado para DAEDALUS.

**`outputs/`.** Guarda resultados regenerables de experimentos y auxiliares
locales de revisión. Está ignorada por Git y cada ejecución puede usar su
propia subcarpeta.

**`outputs/<ejecución>/report.json`.** Registra configuración efectiva,
entorno, métricas y aceptación del experimento correspondiente.

**`outputs/<ejecución>/trajectory.csv`.** Guarda los tiempos y estados de
la trayectoria muestreada, en unidades SI.

**`outputs/<ejecución>/trajectory.png`.** Muestra la trayectoria numérica
frente a la referencia analítica circular.

**`outputs/<ejecución>/consistency.png`.** Muestra derivas de invariantes
y error de posición para inspeccionar la precisión de esa ejecución.

**`build/`.** Contiene copias intermedias generadas al construir el paquete.
Pueden corresponder a una construcción anterior; la fuente editable sigue en `src/`.

**`dist/`.** Guarda distribuciones construidas, como el wheel
`autonomous_orbital_transfer_optimizer-0.1.0-py3-none-any.whl`; su existencia
no demuestra que incluya los últimos cambios locales.

**`src/*.egg-info/`.** Contiene metadatos generados por la instalación o
construcción del paquete: nombre, versión, dependencias y archivos asociados.

**`__pycache__/` y `*.pyc`.** Son cachés de bytecode que Python genera
dentro de varias carpetas para acelerar las importaciones.

## Carpetas hermanas del espacio de trabajo

**`../DAEDALUS-Vault/`.** Es la base científica de Obsidian: reúne teoría,
métodos numéricos, optimización, validación y decisiones del proyecto.
Su navegación comienza en [00 Inicio](../../DAEDALUS-Vault/00%20Inicio.md).

**`../.agents/`.** Guarda las instrucciones y recursos locales de Codex;
el rol de documentación está en `.agents/skills/daedalus-engineer/rol.md`.

**`../.claude/`.** Guarda las instrucciones locales de Claude, incluido
su rol de ingeniería en `.claude/rol.md`.
