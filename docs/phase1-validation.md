# Validación ejecutada — 14 de septiembre de 2026

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
