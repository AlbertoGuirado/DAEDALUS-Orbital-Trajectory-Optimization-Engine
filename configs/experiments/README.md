# Configuraciones exploratorias

Estas configuraciones **no** son la referencia validada de la Fase 1 y no se
ejecutan en CI. La referencia es `configs/earth_orbit.json` (r = 7 000 km,
i = 28,5°, 1 periodo, `max_step` = 60 s).

| Fichero | Escenario | Notas |
| --- | --- | --- |
| `leo_iss_like_1period.json` | LEO circular similar a la ISS: r = 6 778 km (≈ 400 km sobre el radio ecuatorial WGS84 de 6 378,137 km), i = 0,9005899 rad (51,6°), 1 periodo, RK45 | Mismas tolerancias de integración y umbrales de aceptación que la referencia. Órbita circular de dos cuerpos: no modela drag, J2 ni la órbita real de la ISS |

Ejecución:

```powershell
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/experiments/leo_iss_like_1period.json --output outputs/leo_iss_like
```
