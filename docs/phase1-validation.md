# Validación ejecutada — 14 de septiembre de 2026

> La reproducción del baseline y la validación de la ampliación local de
> Kepler se registran al final de este documento por separado.

20 tests de `unittest` pasaron. Cubren los criterios de la Fase 1 documentados
en [phase1.md](phase1.md): soluciones analíticas, invariantes durante diez
periodos, conversiones y degeneraciones, orden RK4, acuerdo entre integradores,
entradas inválidas, fallo del solver y aceptación/rechazo del experimento.

Se construyó correctamente el wheel y se comprobaron los imports de
`orbital_optimizer` y los símbolos públicos de compatibilidad de `physics`.

El ejemplo `configs/earth_orbit.json` se ejecutó con RK45: radio de 7 000 000 m,
inclinación de 28.5°, un periodo de 5828.516637686015 s y 1001 muestras.

| Medida de la ejecución | Valor |
| --- | ---: |
| Error máximo muestreado de posición | 0.015833654949214758 m |
| Error máximo muestreado de velocidad | 1.59444751686348e-5 m/s |
| Deriva relativa máxima de energía | 9.429883984432415e-10 |
| Deriva relativa máxima del vector de momento angular | 4.714942585888468e-10 |
| Deriva absoluta máxima del vector de excentricidad | 1.002666319192962e-9 |
| Evaluaciones de la dinámica | 788 |

Todos los criterios de la configuración se cumplieron. Estos datos corresponden
a este caso y configuración, no son cotas universales ni resultados de optimización.

Entorno: Windows, Python 3.12.1, NumPy 2.3.3, SciPy 1.16.3. Para esta ejecución
se creó `.venv` con `--system-site-packages` y se instaló el proyecto editable
con `--no-deps --no-build-isolation`, reutilizando las dependencias disponibles.
El README incluye la instalación normal en un entorno aislado. La instantánea
de dependencias fija versiones del núcleo; no es un lock multiplataforma con hashes.

Repetir desde la raíz del proyecto:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/earth_orbit.json --output outputs/earth_orbit --plots
```

Artefactos locales generados en `outputs/earth_orbit/` (regenerables e ignorados
por Git): `report.json`, `trajectory.csv`, `trajectory.png` y `consistency.png`.
En GitHub, las ejecuciones del workflow Tests los publican como artefacto
descargable `phase1-earth-orbit`.

## Reproducción de cierre — 16 de octubre de 2026

Sobre `fd9d785` (`v0.1.0`), rama `develop`, se ejecutaron de nuevo los 20 tests
con el Python de `.venv`: **20 OK**, código de salida 0. No se modificaron
código científico ni tests. Se ejecutó el experimento con `--plots` y salida
en `outputs/phase1-closeout-2026-10-16/`: código 0,
`validation_passed=true` y cinco comprobaciones cumplidas.

Las cinco métricas y las 788 evaluaciones coinciden con la tabla anterior.
Entorno efectivo: Python 3.12.1, NumPy 2.3.3, SciPy 1.16.3,
Windows 11 (`Windows-11-10.0.26200-SP0`). Se reutilizó el entorno existente;
no se realizó una instalación limpia ni se reconstruyó el wheel.

Se generaron `report.json`, `trajectory.csv`, `trajectory.png` y
`consistency.png`. El informe registra fecha UTC y configuración completa.
Los artefactos son locales, regenerables e ignorados por Git; no se ha
verificado su publicación en CI. En el baseline, la comparación entre integradores cubría el estado final.
La ampliación siguiente añade comparación en una malla común.

```powershell
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/earth_orbit.json --output outputs/phase1-closeout-2026-10-16 --plots
```

## Ampliación local de Kepler — historial y resultado final

Durante la revisión aparecieron cambios de ingeniería en
`feature/phase-1-kepler`, aún sin commit: `solve_kepler`, `propagate_kepler`,
`tests/test_kepler.py` y dos tests adicionales de propagación. No forman
parte del tag `v0.1.0` ni de la ejecución histórica de 20 tests.

Al ejecutar de nuevo `unittest discover -s tests -v` se observaron
**32 tests, con 1 fallo**:

| Dato | Resultado |
| --- | --- |
| Test | `test_kepler.KeplerPropagationTests.test_three_dimensional_and_degenerate_geometry` |
| Caso | `a=30e6 m`, `e=0.95`, `i=1.2 rad` |
| Error máximo de velocidad | `4.508439805646895e-08 m/s` |
| Límite predefinido | `2.2763579998321872e-08 m/s` |

Los dos tests nuevos de propagación numérica frente a Kepler y de paso por
apsides sí pasaron. Esa primera ejecución no validó la ampliación analítica. Este resultado no invalida
el registro previo del experimento circular. Codex no modifica el algoritmo
ni el límite de aceptación; el historial y la trazabilidad se mantienen en
[phase1-closeout.md](phase1-closeout.md).

### Cierre de la ampliación

Ingeniería revisó el modelo de error del test analítico: el primer límite
solo consideraba evaluación de f/g y raíz; omitía que el estado inicial
redondeado cambia la energía, el semieje mayor y, por tanto, el periodo.
No se modificó el propagador analítico entre ambas ejecuciones. El criterio
de los tests sí cambió y esta revisión no lo presenta como prefijado desde
la primera ejecución.

En `tests/test_kepler.py`, con u=eps_float64, se usa el modelo conservador
de primer orden:

`delta_a/a <= 4*u*(3+e)/(1-e)`;
`timing = 1.5*t_final*4*u*(3+e)/(1-e)`;
`position_limit = 1e-12*a + v_periapsis*timing`;
`velocity_limit = 1e-12*v_periapsis + (mu/r_periapsis**2)*timing`.

Las unidades de los términos son m y m/s, respectivamente. El segundo
término propaga el cambio de periodo como desplazamiento temporal; no es
una garantía universal de error para cualquier escala o e cercano a uno.
La justificación figura en el docstring del test y conserva el registro del
incumplimiento inicial de e=0.95.

La última ejecución de
`.\.venv\Scripts\python.exe -m unittest discover -s tests -v`
terminó con código **0: 32 tests OK**. Además de Kepler, incluye apsides
propagados, vis-viva durante la elipse inclinada y comparación directa
RK4/RK45 en **201 muestras comunes**, tanto circular como elíptica inclinada,
con máximos <1 m y <0.001 m/s. Estas pruebas pasan con los cambios de
ingeniería locales, aún sin commit; no se atribuyen a `v0.1.0`.

Se repitió además el experimento sobre esta ampliación, con código 0,
los mismos cinco valores de métricas y 788 evaluaciones:

```powershell
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/earth_orbit.json --output outputs/phase1-kepler-closeout-2026-10-16 --plots
```

La carpeta nueva conserva el informe y los gráficos de la última ejecución
sin sobrescribir los artefactos de reproducción del baseline.
