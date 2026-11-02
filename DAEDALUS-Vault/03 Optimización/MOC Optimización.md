---
tags:
  - moc
  - algoritmia/optimizacion
---

# MOC Optimización

Volver a [[00 Inicio]].

> [!warning] Todo planificado
> A fecha de la Fase 1 **no existe ningún optimizador** en el código. Los ficheros de `optimization/` están vacíos. Estas notas documentan la teoría que guiará la Fase 3 en adelante.

## Formulación
- [[Formulación del problema de optimización]] — variables, objetivo, restricciones, evaluación
- [[Variables de decisión y normalización]]
- [[Función objetivo]] — Δv total y suma ponderada adimensional
- [[Paisaje de optimización orbital]] — por qué hace falta búsqueda global

## Optimización global (metaheurísticas)
- [[Differential Evolution]] — primer algoritmo del proyecto
- [[Algoritmos genéticos]]
- [[Optimización por enjambre de partículas]]
- [[Teorema No Free Lunch]]

## Optimización local e híbridos
- [[Métodos de optimización local]] — Nelder-Mead, Powell, L-BFGS-B, SLSQP
- [[Algoritmos meméticos]] — global + local, lamarckiano frente a baldwiniano

## Control óptimo
- [[Control óptimo]] — formulación continua
- [[Principio del máximo de Pontryagin]] — y el vector primer de Lawden
- [[Métodos directos de transcripción]] — disparo, colocación, pseudoespectral

## Multiobjetivo
- [[Dominancia de Pareto]]
- [[NSGA-II]]
- [[Hipervolumen]]

## Restricciones
- [[Manejo de restricciones]] — penalizaciones, muerte, reparación
- [[Reglas de factibilidad de Deb]] — punto de partida recomendado
