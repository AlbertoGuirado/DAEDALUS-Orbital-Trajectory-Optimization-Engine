---
tags:
  - algoritmia/numerico
  - algoritmia/raices
aliases:
  - Bisección
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/kepler.py
---

# Método de bisección

Si $g$ es continua y $g(a)\,g(b) < 0$, hay una raíz en $[a, b]$. Se evalúa el punto medio y se descarta la mitad sin cambio de signo.

$$|x_k - x^*| \le \frac{b - a}{2^{k}}$$

| Propiedad | |
|---|---|
| Convergencia | **Lineal** (1 bit por iteración) |
| Requiere | Solo continuidad y un intervalo con cambio de signo |
| Robustez | **Garantizada**: no puede divergir |

## Uso en DAEDALUS

El solver integrado de producción utiliza bisección como salvaguarda de
[[Newton-Raphson]] cuando su propuesta sale del intervalo de raíz.
Los tests mantienen oráculos independientes por bisección pura.

`kepler_reference` resuelve la [[Ecuación de Kepler]] con **60 bisecciones**: el intervalo de $[0, 2\pi]$ se reduce a $2\pi / 2^{60} \approx 5\times10^{-18}$, por debajo del épsilon de float64 → precisión de máquina.

> [!tip] Por qué bisección y no Newton en un oráculo
> En un test, la referencia debe ser **trivialmente correcta** e independiente del código bajo prueba. La lentitud no importa; la garantía de convergencia sí. Ver [[Invariantes como oráculo de validación]].

Relacionado: [[Newton-Raphson]], [[Umbral de degeneración numérica]].
