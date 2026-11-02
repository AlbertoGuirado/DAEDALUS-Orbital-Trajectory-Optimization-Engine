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

Los impulsos se dan como módulos: el resultado es el mismo para una
transferencia ascendente o descendente. Las velocidades salen de la
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
| RF-2.01, RF-2.02 | Impulsos de Hohmann LEO (6 678 km) a GEO (42 164 km) frente a la fórmula cerrada de 03 § 8.2 | Coinciden a 1e-9 m/s |
| RF-2.03, RV-2.02 | Tiempo frente a `pi*sqrt(a³/mu)` calculado de forma independiente | Coincide a 1e-6 s |
| RF-2.04 | Propelente frente a la forma exponencial de Tsiolkovsky | Coincide con rtol 1e-12 |
| RV-2.01 | Dos impulsos aplicados al propagador; órbita final frente a la circular GEO | Error de posición 0,07 m; de velocidad 1,4e-5 m/s; excentricidad < 1e-6 |
| RV-2.03 | Radios iguales | Δv = 0; decrece al acercar `r2` a `r1` |

Sin el segundo impulso, la misma propagación termina con `e = 0,73`, de modo
que la comprobación de RV-2.01 distingue una transferencia completa de una
incompleta.

Referencia de resultados (LEO→GEO, `mu` = `EARTH_MU`): Δv₁ = 2 425,8 m/s,
Δv₂ = 1 466,8 m/s, total 3 892,6 m/s, tiempo 18 990 s. Con 2 000 kg y
`Isp = 320 s`, el propelente es 1 421,5 kg.

Pruebas: `tests/test_transfers.py` y `tests/test_impulses.py`.
