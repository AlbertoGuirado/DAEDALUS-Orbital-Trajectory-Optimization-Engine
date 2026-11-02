---
tags:
  - algoritmia/numerico
  - fisica/fundamentos
aliases:
  - Estado
  - validate_state
  - Reducción a primer orden
estado: implementado
fase: 1
codigo: src/orbital_optimizer/dynamics/two_body.py
---

# Vector de estado

$$\mathbf{y} = [x,\; y,\; z,\; v_x,\; v_y,\; v_z]^{\mathsf{T}}$$

Posición en m, velocidad en m/s, en el [[Marco de referencia inercial]]. **Este orden es contractual** (RNF-07) y no debe alterarse.

## Reducción de orden
Los integradores de EDOs trabajan con sistemas de **primer orden** $\dot{\mathbf{y}} = \mathbf{f}(t, \mathbf{y})$. La ecuación de segundo orden del [[Problema de dos cuerpos]] se convierte en seis de primer orden tratando la velocidad como variable independiente:

$$\dot{\mathbf{r}} = \mathbf{v}, \qquad \dot{\mathbf{v}} = -\mu\mathbf{r}/r^3$$

Esta técnica es general: cualquier EDO de orden $n$ se reduce a $n$ ecuaciones de primer orden.

## Validación de entrada (`validate_state`)
- Dimensión exactamente 6, tipo real, todas las componentes finitas.
- Radio no nulo (singularidad del campo).
- Error → `ValueError`. Se trabaja sobre una **copia** del array del usuario.

## Representaciones alternativas
| Representación | Ventaja | Inconveniente |
|---|---|---|
| Cartesiana (esta) | Sin singularidades, ideal para integrar | Oculta la geometría |
| [[Elementos orbitales clásicos]] | Interpretable, 5 constantes | Singular en circular/ecuatorial |
| Equinocciales modificados | Sin singularidades e interpretable a medias | Menos intuitivos |

En la Fase 5, con [[Encendido de duración finita]], el estado se amplía a 7 componentes con la masa.

Relacionado: [[Runge-Kutta 4]], [[Dormand-Prince RK45]], [[Convención de unidades]].
