---
tags:
  - proyecto
  - validacion
  - resultados
fase: 1
fecha: 2026-09-14
---

# Resultados de validación de la Fase 1

Ejecución registrada en `docs/phase1-validation.md`. **20 tests de `unittest` pasaron.**

## Experimento de referencia
`configs/earth_orbit.json` con [[Dormand-Prince RK45]]: órbita circular de radio 7 000 000 m, inclinación 28.5°, un periodo de 5828.516637686015 s, 1001 muestras.

| Medida | Valor |
|---|---:|
| Error máximo muestreado de posición | 0.0158 m |
| Error máximo muestreado de velocidad | 1.59 × 10⁻⁵ m/s |
| Deriva relativa máxima de [[Energía específica]] | 9.43 × 10⁻¹⁰ |
| Deriva relativa máxima de [[Momento angular específico]] | 4.71 × 10⁻¹⁰ |
| Deriva absoluta máxima de [[Vector de excentricidad]] | 1.00 × 10⁻⁹ |
| Evaluaciones de la dinámica | 788 |

Entorno: Windows, Python 3.12.1, NumPy 2.3.3, SciPy 1.16.3.

> [!caution] Alcance de estos números
> Corresponden a **este caso y esta configuración**. No son cotas universales del solver ni resultados de optimización. Los máximos son sobre las muestras, no sobre el tiempo continuo. Ver [[Niveles de evidencia]].

## Reproducir
```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe -m experiments.validate_two_body --config configs/earth_orbit.json --output outputs/earth_orbit --plots
```

Criterios de aceptación: [[Batería de validación de la Fase 1]].

## Reproducción del 2026-10-16

Sobre `fd9d785` (`v0.1.0`), con el mismo Python 3.12.1, NumPy 2.3.3 y
SciPy 1.16.3, volvieron a pasar los **20 tests** y las cinco comprobaciones
del experimento. Las métricas coinciden con la tabla anterior.
Se generaron JSON, CSV y dos PNG en
`orbital-optimizer/outputs/phase1-closeout-2026-10-16/` (ignorados por Git).
La ejecución reutilizó el entorno local; no reconstruyó el wheel ni comprobó
CI remota. Ver [[Cierre documental de la Fase 1]] para trazabilidad y alcance.

La ampliación local Kepler posterior pasó **32 tests** en su ejecución
final, tras revisar justificadamente el modelo de error del test que falló
inicialmente en e=0.95. El historial se registra en el documento técnico;
los 20 tests y las métricas de esta nota pertenecen al baseline anterior.

## Contexto posterior — 2026-11-02

Las menciones anteriores a ampliaciones locales o Phase 2 futura son históricas.
Kepler y Phase 2 están integrados en `develop` (`27b4116`); la suite actual
pasó 47 tests. No se reemplazan las cifras, fechas ni el historial del criterio
de Kepler de esta nota. Ver [[Validación de la Fase 2 e integración]].
