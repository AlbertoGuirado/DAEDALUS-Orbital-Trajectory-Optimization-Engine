---
tags:
  - fisica/maniobras
aliases:
  - Bielíptica
estado: planificado
fase: 2
---

# Transferencia bielíptica

Tres impulsos y un radio intermedio $r_b > r_2$: se sube primero muy por encima del objetivo y luego se baja.

> [!example] Contraintuitivo
> Para relaciones de radios suficientemente grandes **consume menos Δv que Hohmann**, a costa de un tiempo de transferencia mucho mayor.

| $r_2/r_1$ | Resultado |
|---|---|
| < ≈ 11.94 | Hohmann siempre mejor |
| ≈ 11.94 – 15.58 | Depende del radio intermedio |
| > ≈ 15.58 | Bielíptica mejor con $r_b$ suficientemente grande |

## Por qué importa al proyecto
1. Demuestra que **dos impulsos no siempre es óptimo** → el número de maniobras debe ser un parámetro del problema (RF-3.07), no una constante.
2. Es el caso canónico del **compromiso Δv frente a tiempo** → motiva la [[Dominancia de Pareto|optimización multiobjetivo]].

Requisito RF-2.05 (recomendado).

Relacionado: [[Transferencia de Hohmann]], [[Variables de decisión y normalización]].
