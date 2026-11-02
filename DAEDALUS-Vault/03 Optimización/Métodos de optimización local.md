---
tags:
  - algoritmia/optimizacion
  - algoritmia/local
aliases:
  - Nelder-Mead
  - Powell
  - L-BFGS-B
  - SLSQP
  - Trust-region
  - Optimización local
estado: planificado
fase: 3
---

# Métodos de optimización local

Refinan una solución dentro de su cuenca de atracción. Rápidos y precisos, pero **quedan atrapados** en óptimos locales del [[Paisaje de optimización orbital]].

| Método | Gradiente | Restricciones | Observación |
|---|---|---|---|
| **Nelder-Mead** | No | No directamente | Símplex de $n+1$ puntos que se refleja/expande/contrae. Robusto ante ruido, lento |
| **Powell** | No | No directamente | Direcciones conjugadas sin derivadas |
| **L-BFGS-B** | Sí | Solo límites | Cuasi-Newton con memoria limitada; eficiente con gradiente |
| **SLSQP** | Sí | Igualdades y desigualdades | Programación cuadrática secuencial; **apropiado al problema con restricciones** |
| **Trust-region** | Sí | Sí | Robusto, más caro por iteración |

Todos disponibles en `scipy.optimize.minimize`, ya dependencia del núcleo.

## Obtener gradientes
| Vía | Coste | Estado |
|---|---|---|
| Diferencias finitas | $n+1$ evaluaciones (= propagaciones) por gradiente | Viable |
| Diferenciación automática | Casi gratis… si el propagador es diferenciable | Exige reescribirlo; fuera de alcance |

Requisito RF-3.13 (recomendado). Uso principal: refinamiento en [[Algoritmos meméticos]].

Relacionado: [[Newton-Raphson]], [[Métodos directos de transcripción]].
