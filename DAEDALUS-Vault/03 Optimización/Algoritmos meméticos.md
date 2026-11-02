---
tags:
  - algoritmia/optimizacion
  - algoritmia/metaheuristica
aliases:
  - Híbrido
  - Optimizador híbrido
  - Lamarckiano
  - Baldwiniano
estado: planificado
fase: 3
codigo: optimization/hybrid_optimizer.py (vacío)
---

# Algoritmos meméticos

> [!abstract] En una frase
> La metaheurística localiza la **cuenca** del óptimo global; un método local **refina** dentro de ella.

```mermaid
graph LR
    A[Global: DE / GA / PSO] -->|cuenca prometedora| B[Local: Nelder-Mead / SLSQP]
    B --> C[Solución refinada]
```

## Tres estrategias de acoplamiento
| Estrategia | Qué hace | Compromiso |
|---|---|---|
| **Refinamiento final** | Búsqueda global completa; luego un refinamiento local del mejor | Simple; no aprovecha el refinamiento durante la búsqueda |
| **Lamarckiana** | El individuo refinado **sustituye** al original en la población | Convergencia rápida; riesgo de perder diversidad |
| **Baldwiniana** | El refinamiento mejora la **aptitud evaluada**, pero el genotipo no cambia | Preserva diversidad; aprovecha menos el refinamiento |

> [!success] Recomendación
> Empezar por el **refinamiento final**: es el más simple y el que menos interfiere con la comparación limpia entre metaheurísticas.

Requisito RF-3.14 (recomendado).

Relacionado: [[Differential Evolution]], [[Métodos de optimización local]].
