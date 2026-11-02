---
tags:
  - algoritmia/optimizacion
  - algoritmia/multiobjetivo
aliases:
  - Pareto
  - Frente de Pareto
  - Multiobjetivo
estado: planificado
fase: 3
---

# Dominancia de Pareto

$\mathbf{x}_1$ **domina** a $\mathbf{x}_2$ si es al menos igual de buena en todos los objetivos y estrictamente mejor en alguno:

$$\mathbf{x}_1 \prec \mathbf{x}_2 \iff \forall i:\; J_i(\mathbf{x}_1) \le J_i(\mathbf{x}_2) \;\wedge\; \exists j:\; J_j(\mathbf{x}_1) < J_j(\mathbf{x}_2)$$

- **Conjunto de Pareto**: soluciones no dominadas (espacio de decisión).
- **Frente de Pareto**: su imagen en el espacio de objetivos.

## Por qué importa en DAEDALUS
El compromiso **Δv frente a tiempo** es real e irreducible; la [[Transferencia bielíptica]] lo ilustra: ahorra propelente a costa de mucho más tiempo. No hay una respuesta única: depende de prioridades de misión que el algoritmo no conoce.

> [!quote] Honestidad científica
> Presentar el frente traslada la decisión a quien corresponde — el diseñador de la misión — en lugar de esconderla en unos pesos arbitrarios.

## Suma ponderada frente a Pareto
| | Suma ponderada | Métodos de Pareto |
|---|---|---|
| Implementación | Trivial | Ordenación por dominancia + diversidad |
| Salida | Una solución por vector de pesos | Frente completo en una ejecución |
| Regiones **no convexas** del frente | **Nunca** las alcanza | Sí |
| Parámetros | Pesos difíciles de fijar | Sin pesos |

Requisito RF-3.22 (recomendado). Algoritmo: [[NSGA-II]]. Métrica: [[Hipervolumen]].

Relacionado: [[Función objetivo]].
