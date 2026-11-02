---
tags:
  - fisica/maniobras
aliases:
  - Cambio de inclinación
  - Maniobra combinada
estado: planificado
fase: 2
---

# Cambio de plano

Cambio puro de inclinación $\Delta i$ a rapidez $v$:

$$\Delta v = 2v\sin\frac{\Delta i}{2}$$

> [!warning] Maniobra cara
> 30° de cambio de plano en órbita baja cuesta aproximadamente la mitad que ir a GEO. Regla operativa: hacerlo donde la velocidad es mínima, **en el apoapsis**.

## Maniobra combinada (plano + rapidez)
$$\Delta v = \sqrt{v_1^{2} + v_2^{2} - 2v_1v_2\cos\Delta i}$$

Por la desigualdad triangular es **menor que la suma** de las dos maniobras por separado.

> [!tip] Caso de prueba para el optimizador
> Es una de las primeras no linealidades que un optimizador global puede "descubrir" por sí solo: si combina el cambio de plano con un impulso de altitud, la búsqueda global funciona.

Requisito RF-2.06 (recomendado).

Relacionado: [[Maniobra impulsiva]], [[Ecuación orbital de las cónicas]], [[Paisaje de optimización orbital]].
