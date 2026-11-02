---
tags:
  - fisica/geometria
aliases:
  - Anomalía verdadera
  - Anomalía excéntrica
  - Anomalía media
  - ν
  - E
  - M
estado: implementado
fase: 1
---

# Anomalías orbitales

Tres ángulos que describen dónde está la nave en su órbita, medidos desde el periapsis:

| Anomalía | Símbolo | Naturaleza | Evolución temporal |
|---|---|---|---|
| Verdadera | $\nu$ | Ángulo geométrico real desde el foco | No uniforme (rápida en periapsis) |
| Excéntrica | $E$ | Construcción auxiliar sobre el círculo circunscrito | No uniforme |
| Media | $M$ | Ángulo ficticio | **Lineal en el tiempo** |

$$M = n\,(t - t_p), \qquad n = \sqrt{\mu/a^{3}} \;\text{(movimiento medio)}$$

$$\tan\frac{\nu}{2} = \sqrt{\frac{1+e}{1-e}}\,\tan\frac{E}{2}$$

En implementación se usa la forma con `atan2` para conservar el cuadrante.

## Puente tiempo ↔ posición
Pasar de un tiempo a una posición exige: $t \to M$ (trivial) $\to E$ (**[[Ecuación de Kepler]]**, trascendente) $\to \nu$ (fórmula de arriba).

> [!info] En DAEDALUS
> `OrbitalElements` usa la **anomalía verdadera** (NAIF usa anomalía media; DAEDALUS no). La cadena $t \to M \to E \to \nu$ solo existe en el oráculo de tests (`kepler_reference`).

Relacionado: [[Elementos orbitales clásicos]], [[Tercera ley de Kepler]].
