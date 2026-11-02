---
tags:
  - proyecto
  - validacion
estado: implementado
fase: [1, 2]
---

# Validación de la Fase 2 e integración

Corte **2026-11-02**, código `27b4116` en `develop`. **Implemented:**
dos cuerpos, RK4/RK45, Kepler, Hohmann, propelente ideal e impulsos.
**Planned:** optimización, perturbaciones, restricciones de misión,
bielíptica, cambio de plano y Lambert.

## Evidencia observada

47 tests OK (32 de Phase 1 y 15 de Phase 2), comando
`.\.venv\Scripts\python.exe -B -m unittest discover -s tests -v`.
Python 3.12.1, NumPy 2.3.3, SciPy 1.16.3; Windows 11, instalación existente.
El experimento circular ejecutado en memoria pasó sus cinco checks con
788 evaluaciones; sus métricas coinciden con el registro de Phase 1.

Para radios 6 678 000 m y 42 164 000 m, μ terrestre, RK45 por defecto:

| Medida | Resultado observado |
| --- | ---: |
| Δv total | 3892.6077435913103 m/s |
| Tiempo de transferencia | 18990.051838481286 s |
| Consumo, m0=2000 kg e Isp=320 s | 1421.4766255239515 kg |
| Norma del error final de posición | 0.07199410476581594 m |
| Norma del error final de velocidad | 1.4275188021570504e-5 m/s |
| Excentricidad final | 3.7036434861574957e-9 |
| Excentricidad sin circularización | 0.7265468247853418 |

El procedimiento completo está en `docs/phase1-phase2-audit.md` del
repositorio `orbital-optimizer`. Las normas son mediciones de este caso,
no criterios añadidos a los tests ni cotas universales. Se preservaron los
criterios y el historial del fallo y revisión del modelo de error de Kepler.

## Known Limitations

- Las fórmulas de [[Transferencia de Hohmann]] dan módulos; el llamante
  orienta los impulsos. Invertir radios intercambia sus magnitudes.
- El wrapper de [[Maniobra impulsiva]] devuelve solo el estado final;
  cada arco exige una elipse no radial. No revalida orbitalmente el salto final.
- [[Ecuación de Tsiolkovsky]] calcula masa aparte; no la integra ni verifica
  masa seca o reservas. No existen perturbaciones ni encendidos finitos.
- La inversión de Hohmann tiene tests de coste/tiempo, pero no una prueba
  completa de propagación descendente. No se afirma esa cobertura.
- No se reconstruyeron distribuciones ni se verificaron CI remota o UI.
  El docstring de Hohmann conserva una imprecisión sobre la simetría del
  resultado; queda registrada para ingeniería, sin modificar código.

Relacionado: [[Roadmap por fases]], [[Arquitectura del código]],
[[Cierre documental de la Fase 1]], [[Niveles de evidencia]].
