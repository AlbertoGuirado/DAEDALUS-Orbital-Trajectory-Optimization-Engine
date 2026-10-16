# Fase 1: diseño y contrato de validación

La implementación vigente está en `src/orbital_optimizer`. Los documentos
01–03 describen el proyecto general, sus requisitos y fundamentos; sus referencias
a `physics` se resuelven mediante una capa de imports compatible, sin duplicar
ecuaciones. Los directorios vacíos originales de fases futuras son reservas.

El contrato original corresponde a `v0.1.0`. La ampliación local
`feature/phase-1-kepler` añade solver y propagación analítica; su estado de
validación se registra por separado en [phase1-validation.md](phase1-validation.md).

## Arquitectura y decisiones

```text
src/orbital_optimizer/
  dynamics/two_body.py        constante, validación del estado y campo gravitatorio
  orbital/mechanics.py       invariantes, periodo y costes de impulsos
  orbital/elements.py        elementos clásicos y conversiones
  orbital/kepler.py          ecuación de Kepler (ampliación local)
  simulation/propagator.py    RK4, RK45 y resultado de propagación
  simulation/kepler_propagator.py  coeficientes f/g (ampliación local)
  visualization/plots.py     gráficos científicos opcionales
src/physics/                 compatibilidad con imports anteriores
tests/                       tests unitarios y referencias analíticas
experiments/validate_two_body.py
configs/earth_orbit.json
docs/phase1.md
```

Un espacio de nombres propio evita paquetes genéricos como `dynamics` en el
entorno de Python. El layout `src` exige instalar el paquete, pero comprueba
que los imports funcionan como paquete instalado. Física no importa simulación
ni optimización. No hay clases de optimizadores o misiones sin uso presente.

RK45 usa Dormand–Prince 5(4) de SciPy: adapta el paso pero sus tolerancias son
locales, no una cota del error global. RK4 permite estudiar orden con pasos
prescritos, a cambio de carecer de control automático. Ninguno es simpléctico.
En RK4 las muestras pedidas acortan pasos; estudiar convergencia sin muestreo
explícito. En RK45 las muestras se interpolan sobre la malla adaptativa.

Elementos clásicos son interpretables pero singulares en órbitas circulares
y ecuatoriales. Se adopta una convención explícita en esos casos en lugar de
introducir elementos equinocciales en esta fase.

## Ecuaciones y unidades

Marco inercial dextrógiro centrado en la masa central, sin época o efemérides.
Estado `y=[x,y,z,vx,vy,vz]`: posición en m, velocidad en m/s, tiempo en s.
Todos los ángulos están en radianes. No se confunde radio con altitud.

\[
\dot{\mathbf r}=\mathbf v,\quad
\dot{\mathbf v}=-\mu\mathbf r/\|\mathbf r\|^3,\qquad
\mu_\oplus=3.986004418\times10^{14}\;\mathrm{m^3/s^2}.
\]

\[
\epsilon=\|\mathbf v\|^2/2-\mu/r,\quad
\mathbf h=\mathbf r\times\mathbf v,\quad
\mathbf e=(\mathbf v\times\mathbf h)/\mu-\mathbf r/r.
\]

Energía específica en J/kg, momento angular específico en m²/s, excentricidad
adimensional. El nodo es `n=z_hat × h`. `a=-mu/(2*epsilon)` y
`T=2*pi*sqrt(a³/mu)`. Vis-viva: `v²=mu*(2/r-1/a)`.

En el plano perifocal, con `p=a*(1-e²)`:

\[
\mathbf r_P=\frac p{1+e\cos\nu}(\cos\nu,\sin\nu,0),\qquad
\mathbf v_P=\sqrt{\mu/p}(-\sin\nu,e+\cos\nu,0).
\]

Rotación al marco inercial: `Rz(RAAN) Rx(i) Rz(argumento_periapsis)`.
Los ángulos inversos usan `atan2` y productos orientados respecto a `h`.
Circular: argumento de periapsis cero y anomalía como argumento de latitud.
Ecuatorial: RAAN cero y periapsis medido desde +x en el sentido de movimiento.
Circular ecuatorial: ambos cero y anomalía como longitud en ese sentido;
esto cubre también `i=pi` retrógrado.

Umbral de degeneración: `32*eps_float64`, adimensional, para `e`, `|n|/|h|`
y `|h|/(|r||v|)`. Solo absorbe redondeo; cerca de singularidades se compara
el estado reconstruido, no ángulos individualmente mal condicionados.
Se rechazan órbitas radiales, no ligadas o indistinguibles numéricamente de
parabólicas. La fase valida elipses no radiales (`a>0`, `0<=e<1`). No hay
superficie terrestre, colisiones, propulsión ni perturbaciones.

Coste impulsivo: `delta_v=|v_after-v_before|`, `total_delta_v=sum_k |dv_k|`,
en m/s. La suma de normas no es la norma del vector suma. Aquí solo se calcula
el coste: el propagador es de vuelo libre y no aplica burns.

## Criterios fijados antes de las pruebas

Los límites siguientes son criterios de aceptación para escalas terrestres,
no garantías universales del solver. Se separa posición de velocidad para
no combinar magnitudes con unidades distintas.

| Prueba | Criterio y fundamento |
| --- | --- |
| Circular analítica durante un periodo | Error máximo <1 m, <0.001 m/s: resolución submétrica para este ejemplo |
| Elipse frente a Kepler resuelto por bisección independiente | Mismos límites, incluyendo tiempos intermedios y apsides |
| 10 periodos inclinados, e=0.2 | Error relativo de energía y vector h <1e-7; error absoluto de vector e <1e-7 (evita dividir por e=0) |
| Cierre orbital | <1 m y <0.001 m/s tras un periodo |
| Conversión de elementos | Estado reconstruido con rtol=1e-12, atol=1e-7 m y 1e-10 m/s, margen sobre redondeo |
| RK4 con T/100,T/200,T/400 | Orden observado entre 3.7 y 4.5: entorno asintótico de orden cuatro, admitiendo términos superiores |
| RK4 fino frente a RK45 | Diferencia <1 m y <0.001 m/s |
| Validación de entradas | ValueError para forma, tipo, finitud, dominio y parámetros inválidos |
| Fallos del integrador | RuntimeError, nunca éxito parcial |

RK45 inicial: `rtol=1e-10`, `atol=[1e-3]*3+[1e-6]*3` (m y m/s), paso máximo
60 s. Los límites de aceptación son mayores que las escalas locales para
permitir acumulación; la comparación analítica mide el error global real.
El límite de 60 s es un techo explícito para ejemplos terrestres y se configura.

Experimento determinista: no hay seeds porque no interviene aleatoriedad.
Se guardan configuración efectiva, versiones, estados CSV, métricas JSON y
gráficos opcionales. Los tiempos de cálculo son mediciones individuales.
No se reportan Δv óptimo, convergencia de optimizadores, tasa de éxito de
misiones o Pareto: todavía no existen optimizadores en esta fase.

Referencias: [SciPy solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html),
[NASA/JPL NAIF oscelt](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/oscelt_c.html),
[IERS TN36](https://iers-conventions.obspm.fr/content/tn36.pdf).
NAIF usa km y anomalía media en su interfaz; aquí se usan SI y anomalía verdadera.

## Guía de consulta

[API y contratos de datos](phase1-api.md), [validación ejecutada](phase1-validation.md)
y [auditoría de cierre y trazabilidad](phase1-closeout.md).

## Formulación de la ampliación analítica local

`solve_kepler(M,e)` resuelve `E-e*sin(E)=M`, con anomalías en radianes,
`0<=e<1` y M escalar finito. La derivada `1-e*cos(E)>0` da una raíz única;
el acondicionamiento `dE/dM=1/(1-e*cos(E))` empeora cerca del periapsis
cuando e se aproxima a uno. El solver reduce vueltas internamente y las
restaura al devolver E. Newton se sustituye por bisección si la propuesta
sale del intervalo de raíz. El límite interno es 100 iteraciones; termina
por residuo cero, intervalo o corrección a escala de épsilon float64,
y lanza `RuntimeError` si no converge.

Para propagar desde un estado cartesiano se define
`a=-mu/(2*epsilon)`, `n=sqrt(mu/a^3)`, `c=1-r0/a` y
`s=dot(r0_vector,v0_vector)/sqrt(mu*a)`. El incremento x de anomalía excéntrica
satisface `n*t=x-c*sin(x)+s*(1-cos(x))`; su derivada es `r/a>0`.
La posición y velocidad se obtienen con coeficientes de Lagrange:


\[
f=1-\frac{a}{r_0}(1-\cos x),\qquad
g=t-\frac{x-\sin x}{n},\qquad
r=a(1-c\cos x+s\sin x),
\]
\[
\dot f=-\frac{\sqrt{\mu a}\sin x}{rr_0},\qquad
\dot g=1-\frac{a}{r}(1-\cos x),\qquad
\mathbf r=f\mathbf r_0+g\mathbf v_0,\quad
\mathbf v=\dot f\mathbf r_0+\dot g\mathbf v_0.
\]

Aquí r0 y r sin negrita son normas, en metros; f y gdot son
adimensionales, g está en segundos y fdot en s⁻¹. El código elimina vueltas
completas de tiempo/anomalía antes de evaluar f/g. Evita construir ángulos
para evolucionar el estado, pero valida inicialmente el dominio mediante
`state_to_elements`. El resultado conserva errores de raíz y redondeo;
su denominación analítica no garantiza precisión universal.

Los tests contrastan raíces contra bisección independiente y trayectorias
contra referencias circular y elíptica independientes. El modelo de error analítico de `tests/test_kepler.py` suma un término
de evaluación (`1e-12*a`, `1e-12*v_periapsis`) y otro de desfase temporal
por redondeo de energía y semieje mayor. La fórmula y su revisión tras
el primer fallo están en [phase1-validation.md](phase1-validation.md).
Los resultados observados, incluido cualquier incumplimiento, se registran
en el informe de validación sin cambiar esos criterios desde documentación.
