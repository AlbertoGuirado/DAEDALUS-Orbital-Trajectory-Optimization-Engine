---
tags:
  - algoritmia/optimizacion
  - algoritmia/metaheuristica
aliases:
  - GA
  - Algoritmo genético
estado: planificado
fase: 3
codigo: optimization/genetic_algorithm.py (vacío)
---

# Algoritmos genéticos

Población que evoluciona por **selección, cruce y mutación**, inspirada en la evolución biológica.

| Componente | Opciones habituales |
|---|---|
| Selección | Torneo, ruleta proporcional, por ranking |
| Cruce | Aritmético, **SBX** (simulated binary crossover), un punto, multipunto |
| Mutación | Gaussiana, polinómica, uniforme |
| Elitismo | Preservación explícita de los mejores |

## Frente a [[Differential Evolution]]
| | GA | DE |
|---|---|---|
| Flexibilidad de representación | **Alta**: enteros, permutaciones, mixtos | Vectores reales |
| Nº de parámetros | Muchos | Pocos ($NP$, $F$, $CR$) |
| Sensibilidad al ajuste | Alta | Baja |

> [!tip] Cuándo preferirlo en DAEDALUS
> Si el **número de maniobras** pasa a ser una variable entera del problema (RF-3.07), la representación mixta favorece al GA.

Requisito RF-3.11 (recomendado). Base de [[NSGA-II]] para multiobjetivo.

Relacionado: [[Optimización por enjambre de partículas]], [[Teorema No Free Lunch]].
