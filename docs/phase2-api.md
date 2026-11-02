# Fase 2: baseline analítico (API y resultados)

Alcance: RF-2.01 a RF-2.04 y RV-2.01 a RV-2.03 de
[02-requerimientos.md](02-requerimientos.md). Quedan fuera los requisitos
recomendados u opcionales (bielíptica, cambio de plano, Lambert).
Teoría: [03-fundamentos](03-fundamentos-cientificos-y-algoritmicos.md), secciones 8.2 y 9.1.

Todas las funciones se exportan desde `orbital_optimizer`. Unidades SI:
metros, m/s, segundos, kg, m³/s².

## Transferencia de Hohmann

```python
hohmann_transfer(initial_radius, final_radius, mu=EARTH_MU) -> HohmannTransfer
```

Radios circulares en metros, positivos y finitos. Devuelve un
`HohmannTransfer` congelado con:

| Campo | Significado |
| --- | --- |
| `first_delta_v` | Módulo del primer impulso, m/s |
| `second_delta_v` | Módulo del segundo impulso, m/s |
| `total_delta_v` | Suma de ambos, m/s |
| `transfer_time` | Medio periodo de la elipse de transferencia, s |
| `transfer_semi_major_axis` | `(r1+r2)/2`, m |

Los impulsos se dan como módulos. Al intercambiar los radios se conservan
el coste total y el tiempo, pero se intercambian las magnitudes del primer
y segundo impulso. Las velocidades salen de la
ecuación de vis-viva en las órbitas circulares y en los ápsides de la elipse.

## Ecuación de Tsiolkovsky

```python
propellant_mass(initial_mass, delta_v, specific_impulse) -> float
```

Devuelve la masa de propelente consumida, `m0*(1-exp(-dv/(Isp*g0)))`, en kg.
`m0` en kg, `dv >= 0` en m/s, `Isp > 0` en s. `STANDARD_GRAVITY` es
`g0 = 9.80665 m/s²`, la constante de definición del impulso específico, no la
gravedad local. Entradas no válidas producen `ValueError`.

## Impulsos aplicados al propagador

```python
propagate_with_impulses(initial_state, impulses, duration, mu=EARTH_MU, *,
                        method="RK45", max_step=60.0, rtol=1e-10, atol=None)
```

`impulses` es una secuencia de `(time_s, delta_v)` con `delta_v` de forma
`(3,)` en m/s. Los instantes deben ser estrictamente crecientes y estar en
`[0, duration]`. Cada impulso suma `delta_v` a la velocidad y no cambia la
posición. Entre impulsos se llama a `propagate` con las mismas opciones.
Devuelve el estado `(6,)` en `duration`, incluido un impulso situado en ese
instante. No devuelve la trayectoria intermedia.

## Validación

| Requisito | Prueba | Resultado |
| --- | --- | --- |
| RF-2.01, RF-2.02 | Impulsos de Hohmann LEO (6 678 km) a GEO (42 164 km) frente a la fórmula cerrada de 03 § 8.2 | `assertAlmostEqual`, `places=9` |
| RF-2.03, RV-2.02 | Tiempo frente a `pi*sqrt(a³/mu)` calculado de forma independiente | `assertAlmostEqual`, `places=6` |
| RF-2.04 | Propelente frente a la forma exponencial de Tsiolkovsky | `assert_allclose`, `rtol=1e-12`, `atol` predeterminado |
| RV-2.01 | Dos impulsos aplicados al propagador; órbita final frente a la circular GEO | Error de posición 0,07 m; de velocidad 1,4e-5 m/s; excentricidad < 1e-6 |
| RV-2.03 | Radios iguales | Δv = 0; decrece al acercar `r2` a `r1` |

Sin el segundo impulso, la misma propagación termina con `e = 0,73`, de modo
que la comprobación de RV-2.01 distingue una transferencia completa de una
incompleta.

Referencia de resultados (LEO→GEO, `mu` = `EARTH_MU`): Δv₁ = 2 425,8 m/s,
Δv₂ = 1 466,8 m/s, total 3 892,6 m/s, tiempo 18 990 s. Con 2 000 kg y
`Isp = 320 s`, el propelente es 1 421,5 kg.

Pruebas: `tests/test_transfers.py` y `tests/test_impulses.py`.

Las cifras anteriores se reprodujeron el **2026-11-02** sobre `develop`
(`27b4116`), Python 3.12.1, NumPy 2.3.3 y SciPy 1.16.3. Procedimiento,
resultados sin redondeo de presentación y criterios reales de tests:
[auditoría de integración](phase1-phase2-audit.md). Las mediciones no son
cotas universales del integrador ni resultados de optimización.

## Formulación, signos y unidades

Para radios circulares concéntricos y coplanarios \(r_1,r_2>0\), en m,
y parámetro gravitatorio constante \(\mu>0\), en m³/s², se usa la elipse
tangente de semieje \(a_t=(r_1+r_2)/2\). De la energía específica
\(\epsilon=v^2/2-\mu/r=-\mu/(2a_t)\) se obtiene vis-viva:

\[
v_{c,j}=\sqrt{\mu/r_j},\qquad
v_{t,j}=\sqrt{\mu(2/r_j-1/a_t)},\qquad j\in\{1,2\}.
\]

Los incrementos tangenciales con signo son
\(\delta v_1=v_{t,1}-v_{c,1}\) y
\(\delta v_2=v_{c,2}-v_{t,2}\), ambos en m/s. La API devuelve
\(\Delta v_1=|\delta v_1|\), \(\Delta v_2=|\delta v_2|\) y su suma.
La duración es \(t_H=\pi\sqrt{a_t^3/\mu}\), en s. La excentricidad
teórica del arco es \(e_t=|r_2-r_1|/(r_1+r_2)\); no es un campo del resultado.

En una subida, ambos incrementos son positivos respecto a la velocidad
local; en una bajada son negativos. Los vectores inerciales de los impulsos
no tienen por qué apuntar en la misma dirección. El llamante debe construirlos:
`hohmann_transfer` no genera una secuencia de maniobras orientada.
Si `r1 == r2`, devuelve Δv cero y **medio periodo circular**, no tiempo cero.
Es el resultado de la fórmula implementada, no una optimización del tiempo.

## Contrato de propagación impulsiva

**Implemented:** saltos \(\mathbf r^+=\mathbf r^-\),
\(\mathbf v^+=\mathbf v^-+\Delta\mathbf v\), con cada vector expresado
en los ejes inerciales del estado. El tiempo es transcurrido desde cero,
en s; `duration` debe ser positivo y finito. No hay transformación RTN/LVLH,
ordenación automática ni fusión de impulsos simultáneos.

El wrapper admite `RK4` y `RK45`, no `Kepler`. Cada tramo de duración
positiva llama a `propagate` y hereda su dominio elíptico no radial y sus
tolerancias: `rtol=1e-10`, `atol=[1e-3]*3+[1e-6]*3`, en m y m/s,
y `max_step=60 s`. RK4 valida esas tolerancias pero no las usa para adaptar
el paso. El resultado es un array `(6,)`, no un `PropagationResult`.

**Known Limitations:** un impulso puede dejar un estado fuera del dominio
del siguiente arco y provocar `ValueError`. Un impulso en el instante final
se aplica y devuelve sin una nueva validación de elementos orbitales ni
una comprobación explícita de finitud tras la suma. Por ello, no se promete
que todo estado final devuelto sea una elipse válida. Los tests actuales
verifican la llegada LEO→GEO, no todas las geometrías ni escalas extremas.
El código Hohmann tampoco comprueba explícitamente la finitud de todos los
Δv calculados: validar entradas finitas no garantiza resultados finitos
para cualquier escala float64. Estos límites se documentan sin cambiar código.

## Propelente y supuestos físicos

Con velocidad efectiva de escape constante \(v_e=I_{sp}g_0\), integrar
\(dv=-v_e\,dm/m\) entre \(m_0\) y \(m_f\) da
\(\Delta v=v_e\ln(m_0/m_f)\). Despejando:

\[
m_f=m_0e^{-\Delta v/(I_{sp}g_0)},\qquad m_p=m_0-m_f.
\]

`propellant_mass` evalúa `-m0*expm1(-dv/(Isp*g0))`. Es una evaluación
algebraica ideal con `m0 > 0` en kg e `Isp > 0` en s, ambos finitos,
y Δv escalar finito no negativo. No integra masa, no impone masa seca,
no comprueba reservas y no incluye pérdidas por encendido finito.
Usar el Δv total presupone un mismo Isp para los impulsos; no se suman
consumos calculados repetidamente con la misma masa inicial.

El modelo usa gravedad central de masa puntual, sin perturbaciones,
superficie, colisiones ni encuentro con un vehículo de fase orbital fijada.
Hohmann calcula un baseline; no resuelve una misión ni certifica optimalidad
global frente a otros números de impulsos, tiempos o restricciones.

**Planned:** bielíptica, cambio de plano, Lambert, optimización y restricciones
operativas. No se incluyen en el cierre de los requisitos obligatorios de Phase 2.

## Referencias

- [Orbital Mechanics & Astrodynamics — Hohmann Transfer](https://orbital-mechanics.space/orbital-maneuvers/hohmann-transfer.html): geometría de tangencia y energía del arco.
- [NASA Glenn — Ideal Rocket Equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/): balance de momento y relación de masas ideal.

Consultadas el 2026-11-02. Las derivaciones anteriores explican las ecuaciones;
la evidencia de implementación procede del código y los tests del repositorio.
