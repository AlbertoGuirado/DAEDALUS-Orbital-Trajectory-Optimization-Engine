---
tags:
  - algoritmia/multiobjetivo
  - validacion
aliases:
  - Hypervolume
  - Métricas del frente de Pareto
estado: planificado
fase: 3
---

# Hipervolumen

Volumen del espacio de objetivos **dominado** por el frente aproximado y acotado por un **punto de referencia**.

| Métrica | Qué mide |
|---|---|
| **Hipervolumen** | Convergencia **y** diversidad a la vez |
| Distancia generacional | Distancia media al frente verdadero (si se conoce) |
| Espaciado | Uniformidad del reparto |

> [!important] Por qué es la métrica preferida
> Es la única de uso común **estrictamente monótona respecto a la dominancia de Pareto**: si un frente domina a otro, su hipervolumen es mayor.

> [!warning] El punto de referencia
> Su elección cambia el valor obtenido → debe **fijarse y documentarse** (coherente con RF-3.20). Candidato natural: derivado del baseline de [[Transferencia de Hohmann]].

Requisito RF-3.23 (recomendado).

Relacionado: [[Dominancia de Pareto]], [[NSGA-II]], [[Metodología de benchmark]].
