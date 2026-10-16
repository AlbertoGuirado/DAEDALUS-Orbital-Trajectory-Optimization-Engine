# Autonomous Orbital Transfer Optimizer — Fundamentos Científicos y Algorítmicos

> Documento 3 de 3. Complementarios: `01-caracteristicas-del-proyecto.md` y `02-requerimientos.md`.
>
> **Naturaleza del documento:** referencia teórica. Recoge las ecuaciones, los algoritmos y los criterios de validación que sustentan el proyecto. Toda ecuación aquí recogida es teoría establecida de mecánica orbital y optimización numérica. Los resultados experimentales del proyecto se registran por separado en [phase1-validation.md](phase1-validation.md).

> **Alcance revisado (2026-10-16):** este documento contiene teoría y diseño de todas las fases. La implementación actual se limita a la Fase 1; consultar [API](phase1-api.md) y [cierre](phase1-closeout.md) para capacidades verificadas.

---

## Índice

1. [Marco de referencia, notación y constantes](#1-marco-de-referencia-notación-y-constantes)
2. [El problema de los dos cuerpos](#2-el-problema-de-los-dos-cuerpos)
3. [Invariantes del movimiento](#3-invariantes-del-movimiento)
4. [Elementos orbitales clásicos](#4-elementos-orbitales-clásicos)
5. [La ecuación de Kepler y la propagación analítica](#5-la-ecuación-de-kepler-y-la-propagación-analítica)
6. [Integración numérica](#6-integración-numérica)
7. [Análisis de error y estrategia de validación](#7-análisis-de-error-y-estrategia-de-validación)
8. [Maniobras impulsivas y transferencias analíticas](#8-maniobras-impulsivas-y-transferencias-analíticas)
9. [Propulsión y consumo de propelente](#9-propulsión-y-consumo-de-propelente)
10. [Formulación del problema de optimización](#10-formulación-del-problema-de-optimización)
11. [Algoritmos de optimización global](#11-algoritmos-de-optimización-global)
12. [Optimización local y esquemas híbridos](#12-optimización-local-y-esquemas-híbridos)
13. [Control óptimo](#13-control-óptimo)
14. [Optimización multiobjetivo](#14-optimización-multiobjetivo)
15. [Manejo de restricciones](#15-manejo-de-restricciones)
16. [Perturbaciones orbitales](#16-perturbaciones-orbitales)
17. [Metodología de benchmark](#17-metodología-de-benchmark)
18. [Referencias](#18-referencias)

---

## 1. Marco de referencia, notación y constantes

### 1.1 Marco de referencia

Todo el cálculo se realiza en un **marco inercial geocéntrico ecuatorial**, con origen en el centro de masas del cuerpo central. En el caso terrestre esto corresponde aproximadamente al marco ECI (Earth-Centered Inertial):

- Eje **x**: hacia el equinoccio vernal.
- Eje **z**: hacia el polo norte celeste, perpendicular al plano ecuatorial.
- Eje **y**: completa el triedro directo.

La aproximación inercial es la que legitima aplicar las leyes de Newton sin términos ficticios. Es válida para las escalas temporales de una transferencia orbital típica y no incluye precesión ni nutación del eje terrestre, efectos que quedan fuera del alcance del proyecto.

### 1.2 Notación

| Símbolo | Significado | Unidad |
|---|---|---|
| $\mathbf{r}$ | Vector de posición | m |
| $\mathbf{v}$ | Vector de velocidad | m/s |
| $r = \lVert\mathbf{r}\rVert$ | Distancia radial | m |
| $v = \lVert\mathbf{v}\rVert$ | Rapidez | m/s |
| $\mu$ | Parámetro gravitacional estándar, $GM$ | m³/s² |
| $\mathbf{h}$ | Momento angular específico | m²/s |
| $\varepsilon$ | Energía específica | J/kg |
| $\mathbf{e}$ | Vector de excentricidad | adimensional |
| $a$ | Semieje mayor | m |
| $e$ | Excentricidad | adimensional |
| $i$ | Inclinación | rad |
| $\Omega$ | Ascensión recta del nodo ascendente | rad |
| $\omega$ | Argumento del periapsis | rad |
| $\nu$ | Anomalía verdadera | rad |
| $E$ | Anomalía excéntrica | rad |
| $M$ | Anomalía media | rad |
| $p$ | Semilatus rectum | m |
| $n$ | Movimiento medio | rad/s |
| $T$ | Periodo orbital | s |

### 1.3 Constantes físicas

| Constante | Valor | Unidad | Fuente |
|---|---|---|---|
| $\mu_\oplus$ | 3.986004418 × 10¹⁴ | m³/s² | WGS-84 / EGM-96 |
| $R_\oplus$ | 6 378 137 | m | Radio ecuatorial WGS-84 |
| $J_2$ | 1.08262668 × 10⁻³ | adimensional | EGM-96 |
| $g_0$ | 9.806 65 | m/s² | Aceleración estándar, definición CGPM |

**Nota sobre $g_0$:** no es la gravedad local terrestre sino una constante de definición, empleada exclusivamente para convertir impulso específico en velocidad de escape efectiva. Su uso en la ecuación de Tsiolkovsky es puramente convencional.

---

## 2. El problema de los dos cuerpos

### 2.1 Derivación

Partimos de la ley de gravitación universal de Newton aplicada a dos masas puntuales $M$ y $m$. En el marco del centro de masas, y bajo la hipótesis $m \ll M$, el movimiento relativo del cuerpo menor satisface:

$$\ddot{\mathbf{r}} = -\frac{\mu\,\mathbf{r}}{\lVert\mathbf{r}\rVert^{3}}, \qquad \mu = G(M+m) \approx GM$$

Esta es la ecuación fundamental del proyecto. Es una ecuación diferencial ordinaria vectorial de segundo orden.

### 2.2 Reducción a primer orden

Para la integración numérica se reduce a un sistema de seis ecuaciones de primer orden, definiendo el vector de estado:

$$\mathbf{y} = [x,\; y,\; z,\; v_x,\; v_y,\; v_z]^{\mathsf{T}}$$

con derivada:

$$\dot{\mathbf{y}} = \mathbf{f}(t, \mathbf{y}) = \left[v_x,\; v_y,\; v_z,\; -\frac{\mu x}{r^{3}},\; -\frac{\mu y}{r^{3}},\; -\frac{\mu z}{r^{3}}\right]^{\mathsf{T}}, \qquad r = \sqrt{x^2+y^2+z^2}$$

Esta es la forma que consume `two_body_dynamics` y la que espera `scipy.integrate.solve_ivp`.

### 2.3 Propiedades relevantes del campo vectorial

**Es autónomo.** $\mathbf{f}$ no depende explícitamente de $t$. El argumento temporal se conserva en la firma únicamente por compatibilidad con la interfaz de SciPy.

**Es singular en el origen.** Cuando $r \to 0$ la aceleración diverge como $r^{-2}$. Por esto `validate_state` rechaza un radio nulo y el propagador comprueba que la norma de la posición no se anule tras cada paso.

**Es conservativo.** Deriva del potencial $U(r) = -\mu/r$. Esta propiedad es la base de todos los invariantes de la sección siguiente y, por tanto, de toda la estrategia de validación.

**Tiene solución analítica cerrada.** Las órbitas son secciones cónicas. Esto es excepcional y extraordinariamente valioso: permite validar el integrador numérico contra una verdad exacta, algo que rara vez está disponible en simulación científica.

### 2.4 Condicionamiento numérico

El problema está mal condicionado cerca del periapsis de órbitas muy excéntricas, donde la velocidad y la aceleración crecen bruscamente. Un integrador de paso fijo con paso adecuado para el apoapsis cometerá error grande en el periapsis. Esta es la razón principal para preferir un integrador adaptativo en órbitas excéntricas, y una de las motivaciones del límite de paso máximo configurable en el propagador.

---

## 3. Invariantes del movimiento

Los invariantes son el instrumento central de validación del proyecto. Una cantidad que la física conserva exactamente y que el integrador numérico solo conserva aproximadamente proporciona una medida directa y honesta del error acumulado.

### 3.1 Energía específica

$$\varepsilon = \frac{v^{2}}{2} - \frac{\mu}{r}$$

Es la suma de energía cinética y potencial por unidad de masa. Se relaciona con el semieje mayor por:

$$\varepsilon = -\frac{\mu}{2a}$$

De donde se despeja la **ecuación de vis-viva**, una de las relaciones más útiles de la mecánica orbital:

$$v^{2} = \mu\left(\frac{2}{r} - \frac{1}{a}\right)$$

El signo de $\varepsilon$ clasifica la órbita: negativo para elipse (órbita cerrada), nulo para parábola, positivo para hipérbola.

### 3.2 Momento angular específico

$$\mathbf{h} = \mathbf{r} \times \mathbf{v}$$

Se conserva como vector completo, no solo en módulo. Su conservación en dirección implica que el movimiento es plano, y el plano orbital es el perpendicular a $\mathbf{h}$. Su módulo se relaciona con la geometría de la órbita mediante:

$$h = \sqrt{\mu p}, \qquad p = a(1-e^{2})$$

### 3.3 Vector de excentricidad

También llamado vector de Laplace-Runge-Lenz normalizado:

$$\mathbf{e} = \frac{\mathbf{v} \times \mathbf{h}}{\mu} - \frac{\mathbf{r}}{r} = \frac{1}{\mu}\left[\left(v^{2}-\frac{\mu}{r}\right)\mathbf{r} - (\mathbf{r}\cdot\mathbf{v})\,\mathbf{v}\right]$$

Apunta del foco hacia el periapsis y su módulo es la excentricidad. Su conservación es específica del potencial $1/r$ y se rompe en cuanto se introduce cualquier perturbación. Esto lo convierte en un detector especialmente sensible: **es el primer invariante que se degrada cuando el modelo o el integrador tienen un problema**.

### 3.4 Por qué los tres, y no solo uno

Los tres invariantes detectan familias de error distintas:

| Invariante | Qué detecta especialmente bien |
|---|---|
| Energía | Deriva secular del integrador, errores en la magnitud de la aceleración |
| Momento angular (dirección) | Errores que sacan la trayectoria de su plano |
| Momento angular (módulo) | Errores en la componente transversal de la velocidad |
| Vector de excentricidad | Rotación espuria de la línea de apsides, error en la orientación dentro del plano |

Un integrador puede conservar energía razonablemente bien y aun así hacer precesar la línea de apsides de forma no física. Por eso la validación debe comprobar los tres.

### 3.5 La ecuación orbital

Integrando el problema se obtiene la ecuación de la cónica en coordenadas polares:

$$r(\nu) = \frac{p}{1 + e\cos\nu} = \frac{a(1-e^{2})}{1+e\cos\nu}$$

De donde los radios de apsides:

$$r_p = a(1-e), \qquad r_a = a(1+e)$$

---

## 4. Elementos orbitales clásicos

### 4.1 Motivación

El vector de estado describe la posición instantánea pero oculta la geometría. Los seis elementos orbitales clásicos separan los parámetros que definen la órbita, constantes en el problema de dos cuerpos, del único parámetro que varía con el tiempo.

| Elemento | Qué describe | Constante en dos cuerpos |
|---|---|---|
| $a$ | Tamaño de la órbita | Sí |
| $e$ | Forma de la órbita | Sí |
| $i$ | Inclinación del plano orbital | Sí |
| $\Omega$ | Orientación del plano en el espacio | Sí |
| $\omega$ | Orientación de la órbita dentro del plano | Sí |
| $\nu$ | Posición de la nave en la órbita | **No** |

Esta separación es lo que hace a los elementos orbitales la representación natural para especificar una misión: la órbita inicial y la objetivo se declaran por elementos, no por vectores de estado.

### 4.2 De estado a elementos

Dado $(\mathbf{r}, \mathbf{v})$ y $\mu$:

**1. Momento angular y vector de nodo**

$$\mathbf{h} = \mathbf{r}\times\mathbf{v}, \qquad \mathbf{n} = \hat{\mathbf{z}} \times \mathbf{h} = [-h_y,\; h_x,\; 0]^{\mathsf{T}}$$

**2. Vector de excentricidad y excentricidad**

$$\mathbf{e} = \frac{1}{\mu}\left[\left(v^{2}-\frac{\mu}{r}\right)\mathbf{r} - (\mathbf{r}\cdot\mathbf{v})\,\mathbf{v}\right], \qquad e = \lVert\mathbf{e}\rVert$$

**3. Semieje mayor**

$$\varepsilon = \frac{v^2}{2}-\frac{\mu}{r}, \qquad a = -\frac{\mu}{2\varepsilon}$$

**4. Inclinación**

$$i = \arccos\left(\frac{h_z}{h}\right) \in [0, \pi]$$

**5. Ascensión recta del nodo ascendente**

$$\Omega = \arccos\left(\frac{n_x}{n}\right), \qquad \text{si } n_y < 0 \Rightarrow \Omega \leftarrow 2\pi - \Omega$$

**6. Argumento del periapsis**

$$\omega = \arccos\left(\frac{\mathbf{n}\cdot\mathbf{e}}{n\,e}\right), \qquad \text{si } e_z < 0 \Rightarrow \omega \leftarrow 2\pi - \omega$$

**7. Anomalía verdadera**

$$\nu = \arccos\left(\frac{\mathbf{e}\cdot\mathbf{r}}{e\,r}\right), \qquad \text{si } \mathbf{r}\cdot\mathbf{v} < 0 \Rightarrow \nu \leftarrow 2\pi - \nu$$

Las comprobaciones de signo son necesarias porque el arcocoseno solo devuelve el rango $[0,\pi]$ y hay que resolver la ambigüedad de cuadrante.

### 4.3 Casos degenerados

Este es el punto donde más implementaciones fallan silenciosamente. Tres situaciones rompen las fórmulas anteriores por división por cero:

**Órbita circular ($e \to 0$).** El periapsis no está definido, luego $\omega$ y $\nu$ pierden sentido individualmente. Se sustituyen por el **argumento de latitud**:

$$u = \arccos\left(\frac{\mathbf{n}\cdot\mathbf{r}}{n\,r}\right), \qquad \text{si } r_z < 0 \Rightarrow u \leftarrow 2\pi - u$$

**Órbita ecuatorial ($i \to 0$ o $i \to \pi$).** El nodo ascendente no está definido, luego $\Omega$ y $\omega$ pierden sentido individualmente. Se sustituyen por la **longitud del periapsis**:

$$\tilde{\omega} = \arccos\left(\frac{e_x}{e}\right), \qquad \text{si } e_y < 0 \Rightarrow \tilde{\omega} \leftarrow 2\pi - \tilde{\omega}$$

**Órbita circular y ecuatorial.** Ni periapsis ni nodo están definidos. Solo tiene sentido la **longitud verdadera**:

$$\ell = \arccos\left(\frac{x}{r}\right), \qquad \text{si } y < 0 \Rightarrow \ell \leftarrow 2\pi - \ell$$

Estos casos no son exóticos. Una órbita geoestacionaria es circular y ecuatorial. Una órbita de aparcamiento típica es casi circular. **Ignorar los casos degenerados produce `NaN` en escenarios perfectamente ordinarios**, razón por la cual el requisito RV-1.12 los exige cubiertos por test.

**Alternativa de diseño:** los **elementos equinocciales modificados** eliminan estas singularidades por construcción, a cambio de perder interpretación geométrica directa. La recomendación del proyecto es implementar primero los elementos clásicos con tratamiento explícito de degeneraciones, por ser la representación que el usuario entiende, y considerar los equinocciales más adelante si la optimización sufre por discontinuidades en el espacio de búsqueda.

### 4.4 De elementos a estado

Se construye primero el estado en el **marco perifocal**, cuyo eje $\hat{p}$ apunta al periapsis y cuyo eje $\hat{w}$ es paralelo a $\mathbf{h}$:

$$\mathbf{r}_{pqw} = \frac{p}{1+e\cos\nu}\begin{bmatrix}\cos\nu\\ \sin\nu\\ 0\end{bmatrix}, \qquad \mathbf{v}_{pqw} = \sqrt{\frac{\mu}{p}}\begin{bmatrix}-\sin\nu\\ e+\cos\nu\\ 0\end{bmatrix}$$

Y se rota al marco inercial mediante tres rotaciones elementales sucesivas:

$$\mathbf{R} = \mathbf{R}_z(-\Omega)\,\mathbf{R}_x(-i)\,\mathbf{R}_z(-\omega)$$

$$\mathbf{r} = \mathbf{R}\,\mathbf{r}_{pqw}, \qquad \mathbf{v} = \mathbf{R}\,\mathbf{v}_{pqw}$$

La composición de esta transformación con su inversa debe recuperar la entrada. Ese es el test de ida y vuelta del requisito RV-1.06.

---

## 5. La ecuación de Kepler y la propagación analítica

### 5.1 Anomalías

Tres ángulos describen la posición en la órbita, relacionados entre sí de forma no trivial:

- **Anomalía verdadera** $\nu$: ángulo geométrico real desde el periapsis.
- **Anomalía excéntrica** $E$: construcción auxiliar sobre el círculo circunscrito.
- **Anomalía media** $M$: ángulo ficticio que **crece linealmente con el tiempo**.

$$M = n\,(t - t_p), \qquad n = \sqrt{\frac{\mu}{a^{3}}}$$

### 5.2 La ecuación de Kepler

$$M = E - e\sin E$$

Es **trascendente en $E$**. No tiene solución cerrada, y ese hecho es históricamente el origen de buena parte del análisis numérico. Se resuelve por iteración de Newton-Raphson:

$$E_{k+1} = E_k - \frac{E_k - e\sin E_k - M}{1 - e\cos E_k}$$

Con semilla $E_0 = M$ converge rápidamente para excentricidad moderada. Para $e$ próximo a 1 la convergencia se degrada y conviene una semilla mejor, por ejemplo $E_0 = \pi$, o un método más robusto.

### 5.3 Relación entre anomalía verdadera y excéntrica

$$\tan\frac{\nu}{2} = \sqrt{\frac{1+e}{1-e}}\,\tan\frac{E}{2}$$

En implementación conviene usar la forma con arcotangente de dos argumentos para preservar el cuadrante sin ambigüedad.

### 5.4 Periodo orbital

Tercera ley de Kepler:

$$T = 2\pi\sqrt{\frac{a^{3}}{\mu}}$$

Esta ecuación sostiene uno de los tests de validación más potentes del proyecto: propagar numéricamente durante exactamente un periodo debe devolver el estado inicial. Es un test cerrado y exacto que no requiere ningún valor de referencia externo.

### 5.5 Propagación analítica frente a numérica

| Aspecto | Propagación kepleriana | Integración numérica |
|---|---|---|
| Exactitud | Exacta, salvo error de redondeo | Acumula error de truncamiento |
| Coste | Constante, independiente del intervalo | Proporcional al intervalo |
| Extensibilidad | Solo dos cuerpos puros | Admite cualquier perturbación |
| Uso en el proyecto | Referencia de validación | Motor de simulación |

**Decisión de diseño:** el motor extensible de simulación es numérico, porque las fases 4 y 5 requieren perturbaciones y encendidos finitos. La ampliación local de Fase 1 incorpora una API analítica de dos cuerpos como referencia adicional. Se contrasta con oráculos independientes; su resultado de validación está en [phase1-validation.md](phase1-validation.md). Una solución analítica evaluada en float64 conserva errores de redondeo y de resolución de raíces.

---

## 6. Integración numérica

### 6.1 Runge-Kutta clásico de cuarto orden (RK4)

Con paso $h$ y estado $\mathbf{y}_n$ en el instante $t_n$:

$$\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(t_n,\; \mathbf{y}_n)\\
\mathbf{k}_2 &= \mathbf{f}\!\left(t_n+\tfrac{h}{2},\; \mathbf{y}_n + \tfrac{h}{2}\mathbf{k}_1\right)\\
\mathbf{k}_3 &= \mathbf{f}\!\left(t_n+\tfrac{h}{2},\; \mathbf{y}_n + \tfrac{h}{2}\mathbf{k}_2\right)\\
\mathbf{k}_4 &= \mathbf{f}(t_n+h,\; \mathbf{y}_n + h\,\mathbf{k}_3)\\[2mm]
\mathbf{y}_{n+1} &= \mathbf{y}_n + \frac{h}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)
\end{aligned}$$

**Orden:** el error local de truncamiento es $\mathcal{O}(h^{5})$ y el error global acumulado $\mathcal{O}(h^{4})$. Cuesta cuatro evaluaciones de la derivada por paso.

**Consecuencia verificable:** al reducir el paso a la mitad, el error global debe reducirse en un factor próximo a $2^4 = 16$. Este es el test de orden de convergencia del requisito RV-1.09, y es una comprobación exigente que detecta errores sutiles en la implementación del esquema.

**Ventajas:** simple, determinista, con coste exactamente predecible. Su malla temporal es fija y conocida de antemano, lo cual facilita comparaciones controladas.

**Inconvenientes:** sin control de error. Un paso adecuado en apoapsis puede ser gravemente insuficiente en periapsis de una órbita excéntrica. No es simpléctico, de modo que la energía deriva de forma secular en propagaciones muy largas.

### 6.2 Runge-Kutta-Fehlberg adaptativo, Dormand-Prince 5(4) (RK45)

Es el método por defecto de `solve_ivp`. Emplea siete etapas para producir simultáneamente una solución de quinto orden y una estimación embebida de cuarto orden. La diferencia entre ambas estima el error local, que gobierna la adaptación del paso:

$$h_{\text{nuevo}} = h \cdot \sigma \left(\frac{\text{tol}}{\lVert \mathbf{err}\rVert}\right)^{1/5}$$

donde $\sigma$ es un factor de seguridad menor que uno. La tolerancia combina las componentes relativa y absoluta:

$$\text{tol}_i = \text{atol}_i + \text{rtol}\cdot\lvert y_i \rvert$$

La propiedad FSAL (*First Same As Last*) permite reutilizar la última evaluación del paso anterior, de modo que el coste efectivo es de seis evaluaciones por paso aceptado.

**Por qué las tolerancias absolutas son vectoriales en este proyecto.** Posición y velocidad difieren en varios órdenes de magnitud: una órbita baja terrestre tiene $r \sim 7\times10^{6}$ m y $v \sim 7.5\times10^{3}$ m/s. Una tolerancia absoluta escalar sería absurdamente laxa para la velocidad o absurdamente estricta para la posición. Por eso el propagador usa por defecto 10⁻³ m para posición y 10⁻⁶ m/s para velocidad, es decir, precisión de milímetro en posición y de micrómetro por segundo en velocidad.

**Ventajas:** control de error, eficiencia en órbitas excéntricas al concentrar pasos donde la dinámica es rápida, y robustez.

**Inconvenientes:** malla temporal no predecible de antemano, coste variable, y el control de error es **local**, no global. Una tolerancia local estricta no garantiza por sí sola un error global pequeño, matiz que el propio docstring del propagador ya recoge correctamente.

### 6.3 Comparación y criterio de selección

| Criterio | RK4 | RK45 |
|---|---|---|
| Control de error | No | Sí, local |
| Coste por paso | 4 evaluaciones | 6 efectivas |
| Malla temporal | Fija, predecible | Adaptativa |
| Órbitas casi circulares | Adecuado | Adecuado |
| Órbitas muy excéntricas | Ineficiente o impreciso | Claramente preferible |
| Estudios de convergencia | Ideal, paso prescrito | No aplicable |
| Reproducibilidad exacta | Total | Alta, dependiente de la implementación |

**Recomendación del proyecto:** RK45 como opción por defecto para simulación de misión, RK4 reservado para los estudios de orden de convergencia, donde el paso prescrito es precisamente lo que se quiere controlar.

### 6.4 Un matiz de implementación ya documentado en el código

En el RK4 del propagador, el paso se acorta para aterrizar exactamente sobre cada instante de muestreo solicitado. Esto significa que **pedir muestras altera la malla de integración** y, por tanto, altera el error. Un estudio de convergencia que solicite muestras intermedias medirá algo distinto de lo que pretende medir. El código lo advierte explícitamente y el requisito RF-1.20 lo eleva a norma.

### 6.5 Integradores simplécticos: alternativa no adoptada

Para propagaciones de muchas órbitas, los integradores simplécticos como Verlet de velocidad o los de Yoshida de orden superior preservan la estructura hamiltoniana del problema. Su ventaja decisiva es que **el error de energía permanece acotado y oscilante en lugar de derivar secularmente**.

**Por qué no se adoptan como opción por defecto.** Los simplécticos clásicos requieren paso fijo, lo que los hace ineficientes en órbitas excéntricas, y se complican al introducir fuerzas no conservativas como el empuje y el arrastre, que son precisamente el objeto del proyecto en las fases 4 y 5. La recomendación es no adoptarlos por ahora, medir la deriva de energía de RK45 en los regímenes de interés y reconsiderar únicamente si esa deriva resulta ser limitante. Esta es una decisión que debe tomarse con datos, no a priori.

---

## 7. Análisis de error y estrategia de validación

### 7.1 Taxonomía del error

| Tipo | Origen | Comportamiento |
|---|---|---|
| Truncamiento local | Truncar la serie de Taylor en cada paso | $\mathcal{O}(h^{p+1})$ |
| Truncamiento global | Acumulación de los errores locales | $\mathcal{O}(h^{p})$ |
| Redondeo | Aritmética de punto flotante finita | Crece al reducir $h$, aproximadamente como $\mathcal{O}(1/h)$ |
| Modelado | Física omitida del modelo | Independiente de $h$ |

La competencia entre truncamiento y redondeo implica que **existe un paso óptimo**. Reducir el paso por debajo de él empeora el resultado. Este fenómeno debe ser visible en el estudio de convergencia como un suelo de error, y su aparición es en sí misma una señal de que la implementación es correcta.

### 7.2 Por qué la validación se apoya en invariantes

En simulación científica rara vez se dispone de la solución exacta. El problema de dos cuerpos es la excepción afortunada, y además ofrece tres cantidades conservadas exactamente. Esto permite construir una batería de validación que no depende de ningún valor tabulado externo ni de resultados previos del propio código.

**Principio metodológico vinculante:** un test nunca debe comparar contra valores producidos anteriormente por el propio código. Eso solo verifica que el código no ha cambiado, no que sea correcto. Los tests deben comparar contra solución analítica, contra invariante conservado o contra relación física independiente.

### 7.3 Batería de validación de la Fase 1

| Test | Qué comprueba | Referencia usada |
|---|---|---|
| Deriva de energía | Calidad global del integrador | Invariante exacto |
| Deriva de momento angular | Planaridad y componente transversal | Invariante exacto |
| Deriva del vector de excentricidad | Orientación de la línea de apsides | Invariante exacto |
| Cierre del periodo | Precisión sobre una órbita completa | Tercera ley de Kepler |
| Órbita circular | Caso límite | Geometría exacta |
| Ida y vuelta de elementos | Consistencia de las conversiones | Composición identidad |
| Vis-viva | Relación velocidad-radio | Ecuación analítica |
| Radios de apsides | Geometría de la cónica | $r_p = a(1-e)$, $r_a = a(1+e)$ |
| Orden de convergencia | Corrección del esquema RK4 | Factor 16 al halvar el paso |
| Acuerdo entre integradores | Consistencia mutua | RK4 fino frente a RK45 estricto |

### 7.4 Sobre la elección de tolerancias

Cada tolerancia de test debe justificarse antes de ejecutarlo, a partir del orden del método, el número de pasos y la precisión de la aritmética. **Ajustar una tolerancia hasta que un test pase invalida el test**, porque convierte una comprobación de corrección en una descripción del comportamiento actual, incluido su posible error. El requisito RNF-15 recoge esta norma.

---

## 8. Maniobras impulsivas y transferencias analíticas

### 8.1 El modelo de impulso instantáneo

Se supone que el encendido del motor ocurre en un tiempo despreciable frente al periodo orbital, de modo que la posición no cambia y la velocidad experimenta un salto:

$$\mathbf{r}^{+} = \mathbf{r}^{-}, \qquad \mathbf{v}^{+} = \mathbf{v}^{-} + \Delta\mathbf{v}$$

El coste de una maniobra es $\lVert\Delta\mathbf{v}\rVert$, y el de una secuencia:

$$\Delta v_{\text{total}} = \sum_{k=1}^{N} \lVert\Delta\mathbf{v}_k\rVert$$

**La suma es de módulos, no de vectores.** Dos impulsos opuestos no se cancelan en coste: cuestan el doble. Este detalle, elemental pero fácil de implementar mal, es la diferencia entre una función objetivo correcta y una que premia maniobras absurdas.

**Validez del modelo.** Es excelente para propulsión química de alto empuje, donde un encendido dura decenas de segundos frente a periodos de miles. Es inadecuado para propulsión eléctrica, donde los encendidos duran días. La Fase 5 lo sustituye por encendidos de duración finita.

### 8.2 Transferencia de Hohmann

Es la transferencia de dos impulsos entre órbitas circulares coplanarias, mediante una elipse tangente a ambas. Walter Hohmann la publicó en 1925.

Con radios inicial $r_1$ y final $r_2$, la elipse de transferencia tiene semieje mayor:

$$a_t = \frac{r_1 + r_2}{2}$$

Aplicando vis-viva en cada tangencia:

$$\Delta v_1 = \sqrt{\frac{\mu}{r_1}}\left(\sqrt{\frac{2r_2}{r_1+r_2}} - 1\right)$$

$$\Delta v_2 = \sqrt{\frac{\mu}{r_2}}\left(1 - \sqrt{\frac{2r_1}{r_1+r_2}}\right)$$

$$\Delta v_{\text{total}} = \lvert\Delta v_1\rvert + \lvert\Delta v_2\rvert$$

El tiempo de transferencia es medio periodo de la elipse:

$$t_{\text{transferencia}} = \pi\sqrt{\frac{a_t^{3}}{\mu}}$$

**Su papel en el proyecto.** La transferencia de Hohmann es óptima en Δv entre órbitas circulares coplanarias con dos impulsos, dentro de un rango amplio de relaciones de radios. Esto la convierte en el baseline ideal: **un optimizador que no la alcance está mal implementado o mal configurado**, y eso proporciona una prueba de corrección del motor de optimización, no solo una comparación de rendimiento. Es el contenido del requisito RV-3.01.

### 8.3 Transferencia bielíptica

Usa tres impulsos y un radio intermedio $r_b > r_2$. Contra la intuición, **para relaciones de radios suficientemente grandes consume menos Δv que Hohmann**, a costa de un tiempo de transferencia mucho mayor.

Los umbrales clásicos, presentes en la literatura estándar:

| Relación $r_2/r_1$ | Resultado |
|---|---|
| Menor que ≈ 11.94 | Hohmann es siempre mejor |
| Entre ≈ 11.94 y ≈ 15.58 | Depende del radio intermedio elegido |
| Mayor que ≈ 15.58 | La bielíptica es mejor con radio intermedio suficientemente grande |

**Por qué importa al proyecto.** Es el ejemplo canónico de que la solución de dos impulsos no siempre es óptima, lo cual justifica que el número de maniobras sea un parámetro del problema de optimización y no una constante. También es un excelente caso de prueba del compromiso Δv frente a tiempo, central en la optimización multiobjetivo.

### 8.4 Cambio de plano

Un cambio puro de inclinación $\Delta i$ a rapidez $v$ cuesta:

$$\Delta v = 2v\sin\frac{\Delta i}{2}$$

Es una maniobra **cara**. Un cambio de plano de 30 grados en órbita baja terrestre cuesta aproximadamente la mitad de lo que costaría llegar a la órbita geoestacionaria. De ahí la regla operativa de realizar los cambios de plano donde la velocidad es mínima, es decir, en el apoapsis.

La maniobra combinada de cambio de plano y de rapidez cuesta:

$$\Delta v = \sqrt{v_1^{2} + v_2^{2} - 2v_1v_2\cos\Delta i}$$

que por la desigualdad triangular es **menor que la suma de las dos maniobras separadas**. Esta es una de las primeras no linealidades interesantes que un optimizador puede descubrir por sí mismo, y un buen caso de prueba de que la búsqueda global está funcionando.

### 8.5 El problema de Lambert

Dados dos vectores de posición $\mathbf{r}_1$ y $\mathbf{r}_2$ y un tiempo de vuelo $\Delta t$, determinar la órbita que los conecta. Es un problema de contorno de dos puntos, y su solución da directamente las velocidades requeridas en ambos extremos, de donde se obtienen los dos impulsos.

Su relevancia para el proyecto es estructural: **permite reducir la dimensión del espacio de búsqueda**. En lugar de optimizar sobre magnitudes y direcciones de impulso, se optimiza únicamente sobre instantes de salida y llegada, y Lambert determina el resto. Esto reduce típicamente de seis a dos las variables de decisión de una transferencia de dos impulsos, lo cual transforma la dificultad del problema de optimización.

**Coste de la alternativa:** exige un solucionador de Lambert robusto, incluidos los casos de transferencia múltiple revolución y las configuraciones casi colineales. Es trabajo no trivial, y por eso figura como requisito opcional en la Fase 2.

---

## 9. Propulsión y consumo de propelente

### 9.1 Ecuación de Tsiolkovsky

Derivada de la conservación del momento para un cuerpo de masa variable que expulsa propelente a velocidad de escape efectiva constante:

$$\Delta v = v_e \ln\frac{m_0}{m_f} = I_{sp}\,g_0 \ln\frac{m_0}{m_f}$$

De donde la masa final y el propelente consumido:

$$m_f = m_0\,\exp\!\left(-\frac{\Delta v}{I_{sp}g_0}\right), \qquad m_{\text{prop}} = m_0\left[1 - \exp\!\left(-\frac{\Delta v}{I_{sp}g_0}\right)\right]$$

**La consecuencia central del proyecto:** la relación entre Δv y masa es **exponencial**. Un ahorro porcentual modesto en Δv se traduce en un ahorro desproporcionado de masa de propelente, y por tanto en capacidad de carga útil. Esto es exactamente lo que justifica la existencia de un motor de optimización: en este dominio, optimizar rinde muchísimo más que en uno con relación lineal.

**Nota sobre la monotonía.** Como la exponencial es monótona, minimizar el Δv total equivale a minimizar el propelente **si y solo si el impulso específico es el mismo para todas las maniobras**. Con varios sistemas de propulsión de impulso específico distinto la equivalencia se rompe, y la función objetivo debe expresarse directamente en masa. Este matiz debe quedar registrado en la documentación de la función objetivo.

### 9.2 Encendidos de duración finita

Un encendido real tiene duración $\Delta t$ y produce una aceleración que se añade a la gravitatoria:

$$\ddot{\mathbf{r}} = -\frac{\mu\mathbf{r}}{r^{3}} + \frac{F}{m(t)}\,\hat{\mathbf{u}}(t), \qquad \dot{m} = -\frac{F}{I_{sp}g_0}$$

El sistema pasa de seis a siete dimensiones al incorporar la masa como variable de estado, y deja de ser autónomo, puesto que el empuje depende explícitamente del tiempo a través del programa de encendidos.

**Pérdida por gravedad.** Un encendido finito es siempre menos eficiente que el impulso equivalente, porque parte del empuje se consume contrarrestando la gravedad durante un tiempo no nulo en lugar de traducirse íntegramente en cambio de velocidad útil. La diferencia se denomina pérdida por gravedad y crece con la duración del encendido.

**Test de consistencia entre modelos.** En el límite de empuje muy alto y duración muy corta, el resultado del encendido finito debe converger al del impulso instantáneo. Es el requisito RV-5.03, y es la forma correcta de verificar que el nuevo modelo generaliza al anterior en lugar de contradecirlo.

---

## 10. Formulación del problema de optimización

### 10.1 Forma general

$$\min_{\mathbf{x}\in\mathcal{X}} \; J(\mathbf{x}) \quad \text{sujeto a} \quad g_j(\mathbf{x}) \le 0,\; h_k(\mathbf{x}) = 0, \quad \mathbf{x}_L \le \mathbf{x} \le \mathbf{x}_U$$

Los cuatro elementos que el requisito RF-3.01 a RF-3.04 exige definir explícitamente son: variables de decisión, función objetivo, restricciones y método de evaluación.

### 10.2 Variables de decisión

Para una transferencia de $N$ impulsos, la parametrización directa:

$$\mathbf{x} = \left[t_1, \Delta\mathbf{v}_1, \; t_2, \Delta\mathbf{v}_2, \;\ldots,\; t_N, \Delta\mathbf{v}_N\right]$$

con dimensión $4N$ si cada impulso se describe por magnitud y dos ángulos más su instante, o $4N$ igualmente si se usan tres componentes cartesianas más el instante.

**Dos alternativas de parametrización de la dirección:**

| Parametrización | Ventaja | Inconveniente |
|---|---|---|
| Componentes cartesianas $(\Delta v_x, \Delta v_y, \Delta v_z)$ | Sin singularidades, límites simples | La magnitud no es una variable directa, dificulta imponer un Δv máximo por impulso |
| Esférica: magnitud más dos ángulos | La magnitud es variable directa, límites físicos naturales | Singularidad en los polos, topología periódica de los ángulos |

**Recomendación:** empezar por la parametrización esférica, porque las restricciones de la Fase 4 se expresan naturalmente sobre la magnitud, y tratar la periodicidad angular de forma explícita en el operador de mutación o mediante envoltura del dominio.

**Normalización.** Los instantes se miden en miles de segundos y las magnitudes de impulso en cientos de metros por segundo. Sin normalización, un operador de mutación que use un mismo factor de escala trata ambas variables de forma incoherente. El requisito RF-3.08 exige mapear todas las variables al hipercubo unidad internamente.

### 10.3 Función objetivo

**Fase inicial, objetivo único:**

$$J(\mathbf{x}) = \Delta v_{\text{total}}(\mathbf{x}) = \sum_{k=1}^{N}\lVert\Delta\mathbf{v}_k\rVert$$

**Fase posterior, suma ponderada:**

$$J(\mathbf{x}) = w_1\,\frac{\Delta v}{\Delta v_{\text{ref}}} + w_2\,\frac{t}{t_{\text{ref}}} + w_3\,\frac{E}{E_{\text{ref}}} + w_4\,\frac{R}{R_{\text{ref}}}$$

**Cada término debe adimensionalizarse antes de sumarse.** Sumar metros por segundo con segundos carece de sentido dimensional, y los pesos acabarían absorbiendo silenciosamente una conversión de unidades no documentada. Las escalas de referencia deben declararse, y la elección natural es el valor del baseline de Hohmann para cada magnitud. Esto es el requisito RF-3.21.

**Toda constante debe justificarse.** El requisito RF-3.20 prohíbe introducir pesos o parámetros arbitrarios sin documentar su significado, su unidad y la razón de su valor.

### 10.4 Método de evaluación

La cadena que traduce un vector de decisión en un valor de objetivo:

1. Denormalizar $\mathbf{x}$ desde el hipercubo unidad al dominio físico.
2. Construir el programa de maniobras.
3. Propagar desde el estado inicial, aplicando cada impulso en su instante.
4. Comparar el estado final con la órbita objetivo, obteniendo el residuo.
5. Evaluar las restricciones a lo largo de toda la trayectoria, no solo en los extremos.
6. Calcular el objetivo y aplicar penalizaciones.

**El paso 5 es el más fácil de implementar mal.** La restricción de altitud mínima debe comprobarse en toda la trayectoria: una transferencia puede empezar y terminar en órbitas válidas y aun así atravesar la atmósfera en su periapsis intermedio. El requisito RV-4.02 lo exige verificado por test.

**El coste computacional está dominado por la propagación.** Cada evaluación de la función objetivo requiere integrar numéricamente una trayectoria completa. Esto sitúa el problema en la categoría de **optimización con función objetivo cara**, lo cual condiciona toda la elección algorítmica: el presupuesto de evaluaciones es el recurso escaso, no el tiempo de CPU del optimizador en sí.

### 10.5 Características del paisaje de optimización

| Característica | Implicación algorítmica |
|---|---|
| No convexo, con múltiples óptimos locales | Un método puramente local queda atrapado. Se requiere búsqueda global |
| Objetivo caro de evaluar | El presupuesto de evaluaciones es la métrica de comparación justa |
| Gradientes no disponibles en forma analítica | Favorece métodos sin derivadas, o exige diferenciación numérica costosa |
| Posibles discontinuidades por restricciones | Complica los métodos basados en gradiente |
| Variables de escalas muy dispares | Exige normalización |
| Topología periódica en los ángulos | Requiere operadores que respeten la periodicidad |

Estas características son precisamente las que justifican el uso de metaheurísticas poblacionales, y en particular la elección de Differential Evolution como primer algoritmo.

---

## 11. Algoritmos de optimización global

### 11.1 Differential Evolution

Propuesto por Storn y Price en 1997. Es el primer algoritmo del proyecto por su simplicidad, su robustez y su escaso número de parámetros.

**Idea central:** la mutación no usa una distribución de probabilidad fija, sino **diferencias entre miembros de la propia población**. La escala de la perturbación se adapta así automáticamente a la dispersión actual de la población, que es grande al principio de la búsqueda y pequeña al converger. Esta autoadaptación implícita es lo que hace al método notablemente robusto.

**Mutación, variante DE/rand/1:**

$$\mathbf{v}_i = \mathbf{x}_{r_1} + F\,(\mathbf{x}_{r_2} - \mathbf{x}_{r_3})$$

con $r_1, r_2, r_3$ índices distintos entre sí y distintos de $i$, y $F \in [0,2]$ el factor de escala diferencial.

**Cruce binomial:**

$$u_{i,j} = \begin{cases} v_{i,j} & \text{si } \text{rand}_j \le CR \;\text{ o }\; j = j_{\text{rand}}\\[1mm] x_{i,j} & \text{en otro caso}\end{cases}$$

La condición $j = j_{\text{rand}}$ garantiza que al menos una componente proceda del mutante, evitando que el descendiente sea idéntico al padre.

**Selección voraz:**

$$\mathbf{x}_i^{(g+1)} = \begin{cases}\mathbf{u}_i & \text{si } J(\mathbf{u}_i) \le J(\mathbf{x}_i)\\ \mathbf{x}_i & \text{en otro caso}\end{cases}$$

La selección es elitista por construcción: el mejor individuo nunca se pierde.

**Parámetros y su significado:**

| Parámetro | Rango habitual | Efecto |
|---|---|---|
| $NP$, tamaño de población | 5 a 10 veces la dimensión | Mayor población explora más, consume más evaluaciones por generación |
| $F$, factor diferencial | 0.4 a 0.9 | Mayor valor favorece exploración, menor favorece convergencia |
| $CR$, probabilidad de cruce | 0.1 a 0.9 | Alto valor favorece problemas no separables |

**Variantes de mutación:** DE/best/1 converge más rápido a costa de mayor riesgo de estancamiento. DE/rand/2 explora más. Existen variantes autoadaptativas que ajustan $F$ y $CR$ durante la ejecución. El proyecto debe implementar DE/rand/1 primero por ser la referencia canónica, y comparar variantes por benchmark antes de afirmar cuál conviene.

### 11.2 Algoritmos genéticos

Operan sobre una población mediante selección, cruce y mutación inspirados en la evolución biológica.

| Componente | Opciones habituales |
|---|---|
| Selección | Torneo, ruleta proporcional, por ranking |
| Cruce | Aritmético, SBX (simulado binario), de un punto o multipunto |
| Mutación | Gaussiana, polinómica, uniforme |
| Elitismo | Preservación explícita de los mejores individuos |

**Frente a Differential Evolution:** los algoritmos genéticos ofrecen mayor flexibilidad de operadores, lo cual los hace preferibles cuando la representación no es un vector real, por ejemplo si el número de maniobras es una variable entera del problema. A cambio, tienen más parámetros que ajustar y su rendimiento es más sensible a ese ajuste.

### 11.3 Optimización por enjambre de partículas

Cada partícula mantiene posición y velocidad y se actualiza según:

$$\mathbf{v}_i^{(k+1)} = w\,\mathbf{v}_i^{(k)} + c_1 r_1\left(\mathbf{p}_i - \mathbf{x}_i^{(k)}\right) + c_2 r_2\left(\mathbf{g} - \mathbf{x}_i^{(k)}\right)$$

$$\mathbf{x}_i^{(k+1)} = \mathbf{x}_i^{(k)} + \mathbf{v}_i^{(k+1)}$$

donde $\mathbf{p}_i$ es la mejor posición histórica de la partícula, $\mathbf{g}$ la mejor global, $w$ la inercia y $c_1, c_2$ los coeficientes cognitivo y social.

**Carácter:** converge típicamente rápido, con riesgo asociado de convergencia prematura. El peso de inercia decreciente es el mecanismo habitual para equilibrar exploración inicial y explotación final.

### 11.4 Sobre la comparación entre metaheurísticas

El teorema de **no free lunch** de Wolpert y Macready establece que, promediado sobre todos los problemas posibles, ningún algoritmo de optimización supera a otro. La consecuencia práctica es inequívoca: **la pregunta "cuál es el mejor algoritmo" no tiene respuesta en abstracto**. Solo tiene sentido preguntar cuál funciona mejor en esta clase concreta de problemas, con este presupuesto de evaluaciones y estos criterios.

Esto es exactamente lo que el benchmark del proyecto debe responder, y la razón por la que el principio 2 del proyecto prohíbe afirmar superioridad antes de medirla.

---

## 12. Optimización local y esquemas híbridos

### 12.1 Métodos locales

| Método | Requiere gradiente | Maneja restricciones | Observación |
|---|---|---|---|
| Nelder-Mead | No | No directamente | Robusto ante ruido, converge lentamente |
| Powell | No | No directamente | Direcciones conjugadas sin derivadas |
| L-BFGS-B | Sí | Solo límites | Eficiente con gradiente disponible |
| SLSQP | Sí | Sí, igualdades y desigualdades | Programación cuadrática secuencial, apropiado al problema con restricciones |
| Trust-region | Sí | Sí | Robusto, más costoso por iteración |

Cuando el gradiente no está disponible analíticamente puede obtenerse por diferencias finitas, a coste de $n+1$ evaluaciones adicionales por gradiente, o por diferenciación automática si la cadena de propagación se implementa de forma diferenciable. Esta segunda vía es prometedora pero exige reescribir el propagador, y no se contempla en el alcance actual.

### 12.2 Esquemas híbridos

La combinación habitual, conocida como algoritmo memético, usa la metaheurística para localizar la cuenca de atracción del óptimo global y un método local para refinar dentro de ella.

**Tres estrategias de acoplamiento:**

| Estrategia | Descripción | Compromiso |
|---|---|---|
| Refinamiento final | Búsqueda global completa, luego un refinamiento local del mejor | Simple, no aprovecha el refinamiento durante la búsqueda |
| Lamarckiana | El refinamiento local sustituye al individuo en la población | Convergencia rápida, riesgo de pérdida de diversidad |
| Baldwiniana | El refinamiento mejora la aptitud evaluada pero no el genotipo | Preserva diversidad, aprovecha menos el refinamiento |

El proyecto debe implementar primero el refinamiento final, por ser el más simple y el que menos interfiere con la comparación limpia entre metaheurísticas.

---

## 13. Control óptimo

El control óptimo aborda el problema en su forma continua: hallar la ley de empuje $\mathbf{u}(t)$ que minimiza un funcional de coste sujeto a la dinámica.

$$\min_{\mathbf{u}(\cdot)} \; J = \phi(\mathbf{x}(t_f), t_f) + \int_{t_0}^{t_f} L(\mathbf{x},\mathbf{u},t)\,dt \quad \text{sujeto a} \quad \dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u},t)$$

### 13.1 Métodos indirectos: el principio del máximo de Pontryagin

Se define el hamiltoniano:

$$H = L + \boldsymbol{\lambda}^{\mathsf{T}}\mathbf{f}$$

Las condiciones necesarias de optimalidad son:

$$\dot{\mathbf{x}} = \frac{\partial H}{\partial\boldsymbol{\lambda}}, \qquad \dot{\boldsymbol{\lambda}} = -\frac{\partial H}{\partial\mathbf{x}}, \qquad \mathbf{u}^{*} = \arg\min_{\mathbf{u}\in\mathcal{U}} H$$

**Teoría del vector primer.** Para el problema de mínimo combustible, Lawden definió en 1963 el vector primer $\mathbf{p} = -\boldsymbol{\lambda}_v$, el costado asociado a la velocidad. Sus condiciones necesarias para una trayectoria impulsiva óptima son notablemente concretas:

1. $\lVert\mathbf{p}\rVert \le 1$ en todo instante del intervalo.
2. $\lVert\mathbf{p}\rVert = 1$ exactamente en los instantes de impulso.
3. La derivada de $\lVert\mathbf{p}\rVert$ se anula en los instantes de impulso.
4. La dirección del impulso coincide con la de $\mathbf{p}$.

**Por qué esto importa al proyecto, incluso si no se implementa el método completo.** Estas condiciones permiten **verificar la optimalidad de una solución producida por un optimizador heurístico**, y además indican si añadir un impulso adicional mejoraría el resultado: si el módulo del vector primer excede la unidad en algún punto, la solución no es óptima y existe margen de mejora insertando una maniobra. Es un diagnóstico riguroso que ninguna metaheurística proporciona por sí sola.

**Ventajas de los métodos indirectos:** alta precisión y satisfacción de condiciones necesarias de optimalidad demostrables.

**Inconvenientes:** requieren derivar analíticamente las ecuaciones de costado, presentan un radio de convergencia muy estrecho y exigen una estimación inicial de los costados, magnitudes que carecen de interpretación física intuitiva. Son notoriamente difíciles de hacer converger.

### 13.2 Métodos directos: transcripción

Se discretiza el problema continuo en un problema de programación no lineal de dimensión finita.

| Método | Descripción |
|---|---|
| Disparo simple | Se parametriza el control y se integra toda la trayectoria de una vez |
| Disparo múltiple | Se divide el intervalo en tramos con restricciones de continuidad entre ellos |
| Colocación | Se discretizan estado y control, y la dinámica se impone como restricción algebraica en puntos de colocación |
| Pseudoespectral | Colocación con polinomios globales en nodos de cuadratura de Gauss |

**Ventajas:** radio de convergencia mucho mayor, no requieren derivar costados, manejan restricciones de forma natural.

**Inconvenientes:** menor precisión que los indirectos, y el problema resultante puede ser de dimensión muy alta.

**Recomendación del proyecto.** El control óptimo figura como requisito opcional de la Fase 3 y debe abordarse, si se aborda, por la vía directa. La razón es pragmática: converge sin necesidad de estimaciones iniciales de costado, y encaja con la arquitectura ya prevista al reducirse a un problema de programación no lineal que puede resolverse con las mismas herramientas que el resto del proyecto. El uso más valioso a corto plazo de la teoría indirecta no es resolver, sino **verificar**, mediante las condiciones del vector primer.

---

## 14. Optimización multiobjetivo

### 14.1 Dominancia de Pareto

Una solución $\mathbf{x}_1$ **domina** a $\mathbf{x}_2$ si es al menos igual de buena en todos los objetivos y estrictamente mejor en al menos uno:

$$\mathbf{x}_1 \prec \mathbf{x}_2 \iff \forall i:\; J_i(\mathbf{x}_1) \le J_i(\mathbf{x}_2) \;\;\wedge\;\; \exists j:\; J_j(\mathbf{x}_1) < J_j(\mathbf{x}_2)$$

El **conjunto de Pareto** es el de las soluciones no dominadas, y su imagen en el espacio de objetivos es el **frente de Pareto**.

### 14.2 Por qué importa aquí

El compromiso entre Δv y tiempo de transferencia es real y no reducible. La transferencia bielíptica lo ilustra de forma nítida: ahorra propelente a costa de un tiempo de vuelo mucho mayor. **No existe una respuesta única correcta**: la elección depende de prioridades de misión que el algoritmo no conoce.

Presentar un frente de Pareto en lugar de una única solución traslada esa decisión a quien corresponde, que es el diseñador de la misión, y es científicamente más honesto que ocultarla dentro de unos pesos elegidos arbitrariamente.

### 14.3 Suma ponderada frente a métodos de Pareto

| Aspecto | Suma ponderada | Métodos de Pareto |
|---|---|---|
| Implementación | Trivial, reutiliza el optimizador escalar | Requiere ordenación por dominancia y mecanismos de diversidad |
| Salida | Una solución por vector de pesos | Frente completo en una ejecución |
| Regiones no convexas del frente | **No las alcanza nunca** | Las alcanza |
| Interpretabilidad de los parámetros | Los pesos son difíciles de fijar a priori | No requiere pesos |

La limitación de la suma ponderada frente a regiones no convexas es un resultado teórico establecido, no una cuestión de ajuste. Es la razón por la que el requisito RF-3.22 recomienda implementar optimización multiobjetivo genuina y no conformarse con escalarizar.

### 14.4 NSGA-II

El algoritmo de referencia en optimización multiobjetivo evolutiva. Sus dos mecanismos centrales:

**Ordenación rápida por no dominancia.** Clasifica la población en frentes sucesivos: el primero contiene las soluciones no dominadas, el segundo las no dominadas tras eliminar el primero, y así sucesivamente.

**Distancia de apiñamiento.** Dentro de un mismo frente, prefiere las soluciones situadas en regiones menos pobladas del espacio de objetivos. Esto mantiene el frente bien distribuido en lugar de agrupado en una zona.

### 14.5 Métricas de calidad del frente

| Métrica | Qué mide |
|---|---|
| Hipervolumen | Volumen dominado respecto a un punto de referencia. Captura convergencia y diversidad simultáneamente |
| Distancia generacional | Distancia media al frente verdadero, cuando se conoce |
| Espaciado | Uniformidad de la distribución de soluciones |

El hipervolumen es la métrica preferida por ser la única de uso común que es estrictamente monótona respecto a la dominancia de Pareto. Requiere fijar y documentar el punto de referencia, cuya elección afecta al valor obtenido.

---

## 15. Manejo de restricciones

### 15.1 El problema

Las metaheurísticas operan de forma natural sobre problemas sin restricciones. Incorporarlas requiere una decisión de diseño explícita, y esa decisión afecta materialmente al resultado.

### 15.2 Estrategias

**Penalización estática.**

$$J_{\text{pen}}(\mathbf{x}) = J(\mathbf{x}) + \sum_j \rho_j \max(0, g_j(\mathbf{x}))^{2}$$

Simple, pero el coeficiente $\rho$ es delicado: si es demasiado pequeño, el optimizador acepta soluciones infactibles; si es demasiado grande, el paisaje se vuelve abrupto y la búsqueda se degrada.

**Penalización dinámica.** El coeficiente crece con las generaciones, permitiendo explorar regiones infactibles al principio y forzando factibilidad al final. Introduce un parámetro más de ajuste.

**Penalización de muerte.** Toda solución infactible se descarta. Es simple y no requiere coeficientes, pero desperdicia información y falla cuando la región factible es pequeña o difícil de alcanzar, porque la población inicial puede no contener ninguna solución válida.

**Reglas de factibilidad de Deb.** No usan coeficiente alguno. Al comparar dos soluciones:

1. Si ambas son factibles, gana la de mejor objetivo.
2. Si una es factible y la otra no, gana la factible.
3. Si ambas son infactibles, gana la de menor violación total.

**Es la estrategia recomendada como punto de partida del proyecto**, porque elimina el parámetro arbitrario que el principio de documentación obligaría a justificar y porque su comportamiento es predecible.

**Reparación.** Se proyecta la solución infactible sobre la región factible. Es eficaz cuando la proyección es barata y bien definida, lo cual no es el caso general aquí.

### 15.3 Restricciones del proyecto

| Restricción | Expresión | Naturaleza |
|---|---|---|
| Δv disponible | $\Delta v_{\text{total}} \le \Delta v_{\max}$ | Escalar, evaluable al final |
| Propelente | $m_{\text{prop}} \le m_{\text{prop,disponible}}$ | Escalar, equivalente a la anterior con impulso específico único |
| Tiempo de misión | $t_f - t_0 \le t_{\max}$ | Escalar, evaluable al final |
| Empuje máximo | $F \le F_{\max}$ | Relevante solo con encendidos finitos |
| Altitud mínima | $\min_t r(t) \ge R_\oplus + h_{\min}$ | **De trayectoria: debe evaluarse en todo el intervalo** |
| Masa máxima | $m_0 \le m_{\max}$ | Escalar, de configuración |

La restricción de altitud mínima es cualitativamente distinta de las demás por ser una restricción de trayectoria y no de punto final. Evaluarla solo en los extremos es un error grave y silencioso: la solución parece válida y es físicamente imposible.

---

## 16. Perturbaciones orbitales

### 16.1 Formulación general

El modelo perturbado añade aceleraciones al término kepleriano:

$$\ddot{\mathbf{r}} = -\frac{\mu\mathbf{r}}{r^{3}} + \mathbf{a}_{J_2} + \mathbf{a}_{\text{arrastre}} + \mathbf{a}_{\text{empuje}} + \ldots$$

Esta formulación directa se conoce como **método de Cowell**. Su alternativa, el **método de Encke**, integra únicamente la desviación respecto a una órbita kepleriana de referencia, lo que reduce el error cuando las perturbaciones son pequeñas, a costa de mayor complejidad y de la necesidad de rectificar la referencia periódicamente.

**Recomendación:** adoptar Cowell por su simplicidad y por ser directamente compatible con la arquitectura ya existente del propagador. Encke solo se justificaría si se demuestra una limitación de precisión concreta.

### 16.2 Perturbación J2

El achatamiento terrestre es, con diferencia, la perturbación dominante en órbita baja. Su aceleración:

$$\mathbf{a}_{J_2} = -\frac{3J_2\mu R_\oplus^{2}}{2r^{5}}\begin{bmatrix}x\left(1-\dfrac{5z^{2}}{r^{2}}\right)\\[3mm] y\left(1-\dfrac{5z^{2}}{r^{2}}\right)\\[3mm] z\left(3-\dfrac{5z^{2}}{r^{2}}\right)\end{bmatrix}$$

**Efectos seculares.** J2 produce dos derivas sistemáticas que no promedian a cero:

$$\frac{d\Omega}{dt} = -\frac{3}{2}\,n\,J_2\left(\frac{R_\oplus}{p}\right)^{2}\cos i \qquad \text{(regresión nodal)}$$

$$\frac{d\omega}{dt} = \frac{3}{4}\,n\,J_2\left(\frac{R_\oplus}{p}\right)^{2}\left(5\cos^{2} i - 1\right) \qquad \text{(rotación de apsides)}$$

**Dos consecuencias de ingeniería notables.** La regresión nodal, lejos de ser solo una molestia, permite diseñar **órbitas heliosíncronas**: eligiendo la inclinación adecuada, el plano orbital precesa exactamente un giro por año y la órbita mantiene una geometría constante respecto al Sol. Y la rotación de apsides se anula cuando $5\cos^2 i = 1$, es decir, en la **inclinación crítica** de aproximadamente 63.4 grados, propiedad que explotan las órbitas Molniya para mantener el apoapsis sobre latitudes altas.

**Validación:** las fórmulas seculares proporcionan una referencia analítica independiente contra la que verificar la implementación numérica de J2. Es el requisito RV-5.02.

### 16.3 Arrastre atmosférico

$$\mathbf{a}_{\text{arrastre}} = -\frac{1}{2}\rho\,\frac{C_D A}{m}\,\lVert\mathbf{v}_{\text{rel}}\rVert\,\mathbf{v}_{\text{rel}}$$

donde $\mathbf{v}_{\text{rel}}$ es la velocidad relativa a la atmósfera, que corrota con la Tierra, y $C_D A/m$ es el coeficiente balístico inverso.

**Es no conservativa.** Disipa energía de forma monótona, reduce el semieje mayor y circulariza la órbita. A diferencia de J2, no admite una validación por invariante conservado, sino por comportamiento cualitativo esperado: decaimiento monótono de la energía.

**La dificultad real es la densidad.** $\rho(h)$ varía en órdenes de magnitud con la actividad solar y es la mayor fuente de incertidumbre en predicción orbital en órbita baja. Un modelo exponencial simple basta para estudios comparativos; cualquier pretensión de precisión operativa exigiría un modelo atmosférico completo, fuera del alcance del proyecto. **Esta limitación debe declararse explícitamente junto a cualquier resultado que involucre arrastre.**

### 16.4 Régimen de relevancia

Órdenes de magnitud típicos de la aceleración perturbativa relativa a la kepleriana:

| Perturbación | Régimen de mayor relevancia |
|---|---|
| J2 | Dominante en órbita baja, decreciente con la altitud como $r^{-4}$ frente a $r^{-2}$ |
| Arrastre | Significativo por debajo de unos 600 km, despreciable por encima de unos 1000 km |
| Tercer cuerpo (Luna, Sol) | Creciente con la altitud, relevante en órbita geoestacionaria y superiores |
| Presión de radiación solar | Relevante para relaciones área-masa altas |
| Armónicos de orden superior | Correcciones sobre J2, orden de magnitud menor |

La Fase 5 prevé implementar J2 y arrastre; actualmente ninguno está implementado. Los demás modelos quedan fuera del alcance previsto y esa exclusión deberá constar en cualquier resultado futuro.

---

## 17. Metodología de benchmark

### 17.1 Principio rector

Comparar algoritmos estocásticos requiere disciplina estadística. Una única ejecución de cada algoritmo no demuestra nada: la varianza entre ejecuciones de un mismo algoritmo puede superar ampliamente la diferencia entre algoritmos distintos.

### 17.2 Diseño experimental

| Elemento | Requisito metodológico |
|---|---|
| Semillas | Múltiples ejecuciones independientes por configuración, con semillas registradas |
| Presupuesto | **Idéntico número de evaluaciones de función para todos los algoritmos comparados** |
| Escenarios | Fijados antes de ejecutar, para evitar selección a posteriori |
| Criterio de parada | Presupuesto agotado, no convergencia subjetiva |

**Por qué el presupuesto se mide en evaluaciones y no en tiempo de pared.** El coste está dominado por la propagación numérica, que es común a todos los algoritmos. Comparar por tiempo de pared mide, además del algoritmo, la calidad de su implementación y el estado de la máquina. Comparar por evaluaciones aísla la eficiencia algorítmica, que es lo que se quiere medir. El tiempo de pared se registra igualmente, pero como dato complementario y nunca como criterio de parada.

### 17.3 Métricas

| Métrica | Qué informa |
|---|---|
| Mejor objetivo alcanzado | Calidad de la mejor solución |
| Mediana y rango intercuartílico | Comportamiento típico y dispersión, más robustos que media y desviación |
| Tasa de éxito | Fracción de ejecuciones que alcanzan una solución factible bajo un umbral declarado |
| Evaluaciones hasta el objetivo | Velocidad de convergencia a calidad fijada |
| Curva de convergencia | Mejor valor frente a evaluaciones, promediada entre semillas |
| Violaciones de restricción | Magnitud y tipo |

**Se prefieren mediana y rango intercuartílico a media y desviación típica** porque las distribuciones de resultados de metaheurísticas son marcadamente asimétricas y la media es sensible a valores atípicos.

### 17.4 Tests estadísticos

Las distribuciones de resultados no son normales, por lo que los tests paramétricos no son apropiados.

| Test | Uso |
|---|---|
| Wilcoxon de rangos con signo | Comparar dos algoritmos sobre los mismos problemas, datos emparejados |
| Mann-Whitney U | Comparar dos muestras independientes |
| Friedman | Comparar varios algoritmos sobre varios problemas |
| Post-hoc de Holm o Nemenyi | Comparaciones múltiples tras un Friedman significativo |

**Corrección por comparaciones múltiples.** Comparar cinco algoritmos entre sí supone diez comparaciones por pares. Sin corrección, la probabilidad de obtener al menos un falso positivo al nivel habitual del 5 por ciento supera el 40 por ciento. La corrección no es opcional.

**Tamaño del efecto.** La significación estadística no implica relevancia práctica. Con suficientes ejecuciones, cualquier diferencia arbitrariamente pequeña resulta significativa. Debe reportarse también la magnitud de la diferencia, y valorarse si es operativamente relevante.

### 17.5 Presentación de resultados

| Gráfica | Contenido |
|---|---|
| Trayectoria orbital | Órbita inicial, de transferencia y objetivo, con los impulsos marcados |
| Curva de convergencia | Mejor objetivo frente a evaluaciones, con banda intercuartílica entre semillas |
| Δv frente a tiempo de cómputo | Compromiso entre calidad y coste computacional |
| Comparación de optimizadores | Diagramas de caja del objetivo final por algoritmo |
| Frente de Pareto | Objetivos enfrentados, con el baseline de Hohmann como punto de referencia |

Toda gráfica lleva ejes etiquetados con unidades. Toda curva agregada sobre varias semillas indica explícitamente el estadístico representado y la dispersión.

### 17.6 Amenazas a la validez

| Amenaza | Mitigación |
|---|---|
| Ajuste de hiperparámetros desigual entre algoritmos | Mismo esfuerzo de ajuste para todos, o valores por defecto de la literatura en todos |
| Selección de escenarios favorable | Escenarios fijados antes de ejecutar |
| Semillas insuficientes | Número de semillas fijado a priori, no ampliado hasta obtener el resultado deseado |
| Conclusión generalizada en exceso | Declarar explícitamente la clase de problemas sobre la que se concluye |
| Comparación contra implementaciones ajenas mal configuradas | Usar implementaciones de referencia con configuración documentada |

---

## 18. Referencias

**Mecánica orbital**

- Vallado, D. A. *Fundamentals of Astrodynamics and Applications*. Referencia estándar para conversiones de elementos, casos degenerados y modelos de perturbación.
- Curtis, H. D. *Orbital Mechanics for Engineering Students*. Tratamiento pedagógico de las transferencias y los elementos orbitales.
- Battin, R. H. *An Introduction to the Mathematics and Methods of Astrodynamics*. Tratamiento matemático profundo, incluido el problema de Lambert.
- Bate, Mueller y White. *Fundamentals of Astrodynamics*. Texto clásico y conciso.

**Control óptimo y optimización de trayectorias**

- Lawden, D. F. *Optimal Trajectories for Space Navigation*, 1963. Origen de la teoría del vector primer.
- Betts, J. T. *Practical Methods for Optimal Control Using Nonlinear Programming*. Referencia de métodos directos.
- Conway, B. A. (ed.) *Spacecraft Trajectory Optimization*. Panorámica del estado del arte.

**Optimización numérica y metaheurísticas**

- Storn, R. y Price, K. "Differential Evolution: A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces", 1997.
- Deb, K. et al. "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II", 2002.
- Deb, K. "An Efficient Constraint Handling Method for Genetic Algorithms", 2000. Origen de las reglas de factibilidad.
- Nocedal, J. y Wright, S. *Numerical Optimization*. Referencia de métodos locales.
- Wolpert, D. y Macready, W. "No Free Lunch Theorems for Optimization", 1997.

**Análisis numérico**

- Hairer, Nørsett y Wanner. *Solving Ordinary Differential Equations I: Nonstiff Problems*. Referencia de métodos Runge-Kutta y control de paso.
- Hairer, Lubich y Wanner. *Geometric Numerical Integration*. Referencia de métodos simplécticos y conservación de invariantes.

**Metodología experimental**

- Derrac, J. et al. "A Practical Tutorial on the Use of Nonparametric Statistical Tests as a Methodology for Comparing Evolutionary and Swarm Intelligence Algorithms", 2011.

---

## Apéndice: tabla de referencia rápida de ecuaciones

| Concepto | Ecuación |
|---|---|
| Dinámica de dos cuerpos | $\ddot{\mathbf{r}} = -\mu\mathbf{r}/r^{3}$ |
| Energía específica | $\varepsilon = v^{2}/2 - \mu/r = -\mu/(2a)$ |
| Vis-viva | $v^{2} = \mu(2/r - 1/a)$ |
| Momento angular | $\mathbf{h} = \mathbf{r}\times\mathbf{v}$, $h = \sqrt{\mu p}$ |
| Vector de excentricidad | $\mathbf{e} = (\mathbf{v}\times\mathbf{h})/\mu - \mathbf{r}/r$ |
| Semilatus rectum | $p = a(1-e^{2})$ |
| Ecuación orbital | $r = p/(1+e\cos\nu)$ |
| Radios de apsides | $r_p = a(1-e)$, $r_a = a(1+e)$ |
| Periodo | $T = 2\pi\sqrt{a^{3}/\mu}$ |
| Movimiento medio | $n = \sqrt{\mu/a^{3}}$ |
| Ecuación de Kepler | $M = E - e\sin E$ |
| Anomalías | $\tan(\nu/2) = \sqrt{(1+e)/(1-e)}\,\tan(E/2)$ |
| RK4 | $\mathbf{y}_{n+1} = \mathbf{y}_n + \tfrac{h}{6}(\mathbf{k}_1+2\mathbf{k}_2+2\mathbf{k}_3+\mathbf{k}_4)$ |
| Hohmann, primer impulso | $\Delta v_1 = \sqrt{\mu/r_1}\left(\sqrt{2r_2/(r_1+r_2)}-1\right)$ |
| Hohmann, segundo impulso | $\Delta v_2 = \sqrt{\mu/r_2}\left(1-\sqrt{2r_1/(r_1+r_2)}\right)$ |
| Hohmann, tiempo | $t = \pi\sqrt{a_t^{3}/\mu}$, $a_t = (r_1+r_2)/2$ |
| Cambio de plano | $\Delta v = 2v\sin(\Delta i/2)$ |
| Maniobra combinada | $\Delta v = \sqrt{v_1^{2}+v_2^{2}-2v_1v_2\cos\Delta i}$ |
| Tsiolkovsky | $\Delta v = I_{sp}g_0\ln(m_0/m_f)$ |
| Regresión nodal J2 | $\dot{\Omega} = -\tfrac{3}{2}nJ_2(R/p)^{2}\cos i$ |
| Rotación de apsides J2 | $\dot{\omega} = \tfrac{3}{4}nJ_2(R/p)^{2}(5\cos^{2}i-1)$ |
| Arrastre | $\mathbf{a} = -\tfrac{1}{2}\rho(C_DA/m)\lVert\mathbf{v}_{\text{rel}}\rVert\mathbf{v}_{\text{rel}}$ |
| Mutación DE | $\mathbf{v}_i = \mathbf{x}_{r_1} + F(\mathbf{x}_{r_2}-\mathbf{x}_{r_3})$ |
| Actualización PSO | $\mathbf{v}^{k+1} = w\mathbf{v}^{k} + c_1r_1(\mathbf{p}-\mathbf{x}) + c_2r_2(\mathbf{g}-\mathbf{x})$ |
