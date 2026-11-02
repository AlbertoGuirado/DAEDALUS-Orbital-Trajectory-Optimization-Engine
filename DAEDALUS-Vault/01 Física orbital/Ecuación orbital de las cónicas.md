---
tags:
  - fisica/geometria
aliases:
  - Ecuación de la cónica
  - Radios de apsides
  - Periapsis
  - Apoapsis
estado: implementado
fase: 1
---

# Ecuación orbital de las cónicas

Integrando el [[Problema de dos cuerpos]] se obtiene la trayectoria en polares respecto al foco:

$$r(\nu) = \frac{p}{1 + e\cos\nu} = \frac{a(1-e^{2})}{1+e\cos\nu}$$

con $\nu$ la anomalía verdadera (ver [[Anomalías orbitales]]).

## Radios de apsides
$$r_p = a(1-e) \quad\text{(periapsis, } \nu = 0\text{)}, \qquad r_a = a(1+e) \quad\text{(apoapsis, } \nu = \pi\text{)}$$

- En el periapsis la velocidad es **máxima**; en el apoapsis, **mínima** (de ahí que los [[Cambio de plano|cambios de plano]] se hagan en el apoapsis).
- Test RV-1.08: los radios de apsides de la trayectoria propagada coinciden con los analíticos (`test_known_periapsis_and_apsides`).

## Familia de cónicas
| $e$ | Curva |
|---|---|
| 0 | Circunferencia |
| (0, 1) | Elipse |
| 1 | Parábola |
| > 1 | Hipérbola |

Relacionado: [[Energía específica]], [[Elementos orbitales clásicos]], [[Marco perifocal]].
