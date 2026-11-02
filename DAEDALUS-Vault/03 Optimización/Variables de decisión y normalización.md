---
tags:
  - algoritmia/optimizacion
aliases:
  - Variables de decisión
  - Normalización
  - Hipercubo unidad
  - Parametrización
estado: planificado
fase: 3
---

# Variables de decisión y normalización

## Parametrización directa de $N$ impulsos
$$\mathbf{x} = \left[t_1, \Delta\mathbf{v}_1, \; t_2, \Delta\mathbf{v}_2, \;\ldots,\; t_N, \Delta\mathbf{v}_N\right] \qquad \dim = 4N$$

## Dirección del impulso: dos opciones
| Parametrización | Ventaja | Inconveniente |
|---|---|---|
| Cartesiana $(\Delta v_x, \Delta v_y, \Delta v_z)$ | Sin singularidades, límites simples | La magnitud no es variable directa → difícil imponer Δv máximo por impulso |
| **Esférica** (magnitud + 2 ángulos) | Magnitud directa, límites físicos naturales | Singularidad en los polos; ángulos periódicos |

> [!success] Recomendación
> Esférica, porque las restricciones de la Fase 4 se expresan sobre la magnitud. La periodicidad angular se trata **explícitamente** (envoltura del dominio o en el operador de mutación).

## Normalización (RF-3.08)
Los instantes van en miles de segundos; las magnitudes, en cientos de m/s. Un operador de mutación con un único factor de escala las trataría de forma incoherente. Todas las variables se mapean internamente al **hipercubo unidad** $[0,1]^n$:

$$x_{\text{norm}} = \frac{x - x_L}{x_U - x_L}$$

## Reducción de dimensión
- [[Problema de Lambert]]: optimizar solo instantes de salida/llegada → de 6 a 2 variables para dos impulsos.
- Número de maniobras como **parámetro** (RF-3.07), motivado por la [[Transferencia bielíptica]]. Si se vuelve variable entera, favorece los [[Algoritmos genéticos]].

Relacionado: [[Formulación del problema de optimización]], [[Paisaje de optimización orbital]].
