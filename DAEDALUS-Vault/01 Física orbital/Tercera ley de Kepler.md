---
tags:
  - fisica/geometria
aliases:
  - Periodo orbital
  - orbital_period
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/mechanics.py
---

# Tercera ley de Kepler

$$T = 2\pi\sqrt{\frac{a^{3}}{\mu}}$$

El periodo depende **solo del semieje mayor**, no de la excentricidad. Ejemplo validado: $a = 7\,000\,000$ m → $T \approx 5828.5$ s (~97 min).

## El test de cierre orbital
Propagar exactamente un periodo debe devolver el estado inicial. Es un test **cerrado y exacto** que no necesita ningún valor de referencia externo. Criterio Fase 1: < 1 m y < 0.001 m/s tras un periodo (RV-1.04).

## Usos futuros
- Tiempo de la [[Transferencia de Hohmann]]: medio periodo de la elipse de transferencia, $t = \pi\sqrt{a_t^3/\mu}$.
- Movimiento medio $n = 2\pi/T$ en las tasas seculares de la [[Perturbación J2]].

Relacionado: [[Anomalías orbitales]], [[Batería de validación de la Fase 1]].
