---
tags:
  - fisica/geometria
  - algoritmia/numerico
aliases:
  - Singularidades de los elementos
  - Órbita circular
  - Órbita ecuatorial
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/elements.py
---

# Casos degenerados de los elementos orbitales

> [!danger] Donde más implementaciones fallan en silencio
> Las fórmulas clásicas dividen por $e$ o por $n = |\mathbf{n}|$. En los casos degenerados producen `NaN` en escenarios **totalmente ordinarios**: una órbita geoestacionaria es circular y ecuatorial; una órbita de aparcamiento es casi circular.

| Caso | Qué se indefine | Sustituto | Convención DAEDALUS |
|---|---|---|---|
| **Circular** ($e \to 0$) | Periapsis → $\omega$ y $\nu$ | Argumento de latitud $u$ | $\omega = 0$, anomalía = $u$ |
| **Ecuatorial** ($i \to 0$ o $\pi$) | Nodo → $\Omega$ y $\omega$ | Longitud del periapsis $\tilde\omega$ | $\Omega = 0$, periapsis medido desde +x en el sentido del movimiento |
| **Circular y ecuatorial** | Periapsis y nodo | Longitud verdadera $\ell$ | $\Omega = \omega = 0$, anomalía = longitud; cubre también $i = \pi$ retrógrado |

## Detección
Umbral `SINGULAR_TOL = 32·eps_float64` (adimensional) aplicado a $e$, $|\mathbf{n}|/|\mathbf{h}|$ y $|\mathbf{h}|/(|\mathbf{r}||\mathbf{v}|)$. Solo absorbe redondeo; ver [[Umbral de degeneración numérica]].

## Validación
RV-1.12: los tests cubren órbita circular inclinada, elíptica ecuatorial y circular ecuatorial (`test_roundtrip_general_and_degenerate_geometry`).

Relacionado: [[Elementos orbitales clásicos]], [[Marco perifocal]].
