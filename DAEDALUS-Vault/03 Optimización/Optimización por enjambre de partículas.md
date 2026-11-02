---
tags:
  - algoritmia/optimizacion
  - algoritmia/metaheuristica
aliases:
  - PSO
  - Particle Swarm
estado: planificado
fase: 3
---

# Optimización por enjambre de partículas

Cada partícula tiene posición $\mathbf{x}_i$ y velocidad $\mathbf{v}_i$, y recuerda su mejor posición $\mathbf{p}_i$; el enjambre comparte la mejor global $\mathbf{g}$.

$$\mathbf{v}_i^{(k+1)} = \underbrace{w\,\mathbf{v}_i^{(k)}}_{\text{inercia}} + \underbrace{c_1 r_1\left(\mathbf{p}_i - \mathbf{x}_i^{(k)}\right)}_{\text{cognitivo}} + \underbrace{c_2 r_2\left(\mathbf{g} - \mathbf{x}_i^{(k)}\right)}_{\text{social}}$$

$$\mathbf{x}_i^{(k+1)} = \mathbf{x}_i^{(k)} + \mathbf{v}_i^{(k+1)}$$

con $r_1, r_2 \sim U(0,1)$.

## Carácter
- Converge típicamente **rápido**…
- …con riesgo de **convergencia prematura** (todo el enjambre colapsa sobre un óptimo local).
- Mecanismo de equilibrio: **inercia $w$ decreciente** — exploración al inicio, explotación al final.

> [!note] "Velocidad" aquí no es física
> Es la velocidad de la partícula en el espacio de búsqueda normalizado, no una velocidad orbital.

Requisito RF-3.12 (recomendado).

Relacionado: [[Differential Evolution]], [[Algoritmos genéticos]], [[Teorema No Free Lunch]].
