---
tags:
  - algoritmia/optimizacion
  - algoritmia/multiobjetivo
  - algoritmia/metaheuristica
aliases:
  - NSGA2
estado: planificado
fase: 3
---

# NSGA-II

Algoritmo de referencia en optimización multiobjetivo evolutiva (Deb et al., 2002). Es un [[Algoritmos genéticos|algoritmo genético]] con selección basada en [[Dominancia de Pareto]].

## Mecanismo 1: ordenación rápida por no dominancia
Clasifica la población en **frentes sucesivos**:
- Frente 1: las no dominadas.
- Frente 2: las no dominadas tras quitar el frente 1.
- …y así sucesivamente. Menor frente = mejor rango.

## Mecanismo 2: distancia de apiñamiento
Dentro de un mismo frente, prefiere las soluciones en regiones **menos pobladas** del espacio de objetivos. Mantiene el frente bien repartido en lugar de amontonado en una zona.

## Selección
Comparación por (rango, −apiñamiento): primero mejor frente; a igualdad, más aislada. Con **elitismo**: padres e hijos compiten juntos por la siguiente generación.

Relacionado: [[Hipervolumen]], [[Reglas de factibilidad de Deb]] (mismo autor; se combinan de forma natural).
