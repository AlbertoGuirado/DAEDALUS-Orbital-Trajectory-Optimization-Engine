---
tags:
  - fisica/invariantes
aliases:
  - h
  - specific_angular_momentum
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/mechanics.py
---

# Momento angular específico

$$\mathbf{h} = \mathbf{r} \times \mathbf{v} \qquad [\text{m}^2/\text{s}]$$

Se conserva **como vector completo**, no solo en módulo:
- **Dirección constante** → el movimiento es plano; el plano orbital es perpendicular a $\mathbf{h}$.
- **Módulo** ligado a la geometría:

$$h = \sqrt{\mu p}, \qquad p = a(1-e^{2})$$

donde $p$ es el *semilatus rectum*.

## Usos en el código
- Inclinación: $i = \arccos(h_z / h)$.
- Vector de nodo: $\mathbf{n} = \hat{\mathbf{z}} \times \mathbf{h}$ (`node_vector`).
- Detección de órbitas radiales: se rechazan si $|\mathbf{h}|/(|\mathbf{r}||\mathbf{v}|)$ está por debajo del [[Umbral de degeneración numérica]].

## Como herramienta de validación
| Componente | Detecta |
|---|---|
| Dirección | Errores que sacan la trayectoria de su plano |
| Módulo | Errores en la componente transversal de la velocidad |

Criterio Fase 1: error relativo del vector $\mathbf{h}$ < 10⁻⁷ en 10 periodos.

Relacionado: [[Elementos orbitales clásicos]], [[Invariantes como oráculo de validación]].
