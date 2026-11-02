---
tags:
  - algoritmia/optimizacion
aliases:
  - J
  - Objetivo
  - Suma ponderada
estado: planificado
fase: 3
---

# Función objetivo

## Fase inicial: objetivo único
$$J(\mathbf{x}) = \Delta v_{\text{total}} = \sum_{k=1}^{N}\lVert\Delta\mathbf{v}_k\rVert$$

Suma de **módulos** (ver [[Maniobra impulsiva]]). Equivale a minimizar propelente **solo si todos los impulsos tienen el mismo $I_{sp}$** ([[Ecuación de Tsiolkovsky]]).

## Fase posterior: suma ponderada
$$J = w_1\,\frac{\Delta v}{\Delta v_{\text{ref}}} + w_2\,\frac{t}{t_{\text{ref}}} + w_3\,\frac{E}{E_{\text{ref}}} + w_4\,\frac{R}{R_{\text{ref}}}$$

> [!danger] Adimensionalizar antes de sumar (RF-3.21)
> Sumar m/s con segundos no tiene sentido dimensional; los pesos acabarían absorbiendo una conversión de unidades no documentada. Escalas de referencia naturales: los valores del baseline de [[Transferencia de Hohmann]].

> [!warning] Ninguna constante arbitraria (RF-3.20)
> Todo peso documenta su significado, unidad y la razón de su valor.

## Límite de la suma ponderada
**Nunca** alcanza las regiones no convexas del frente de Pareto, por mucho que se ajusten los pesos (resultado teórico). Por eso se recomienda multiobjetivo genuino: [[Dominancia de Pareto]], [[NSGA-II]].

Relacionado: [[Formulación del problema de optimización]], [[Manejo de restricciones]].
