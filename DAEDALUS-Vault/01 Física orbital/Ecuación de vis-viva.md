---
tags:
  - fisica/invariantes
aliases:
  - Vis-viva
estado: implementado
fase: 1
---

# Ecuación de vis-viva

$$v^{2} = \mu\left(\frac{2}{r} - \frac{1}{a}\right)$$

Se obtiene igualando las dos expresiones de la [[Energía específica]]. Da la rapidez en cualquier punto de la órbita conociendo solo el radio y el semieje mayor.

## Casos particulares
- Órbita circular ($r = a$): $v_c = \sqrt{\mu/r}$
- Velocidad de escape ($a \to \infty$): $v_{esc} = \sqrt{2\mu/r}$

## Papel en DAEDALUS
- Es la herramienta con la que se derivan los impulsos de la [[Transferencia de Hohmann]]: se aplica en cada punto de tangencia.
- Es un test de validación de la Fase 1 (RV-1.07): la velocidad calculada a un radio dado coincide con la predicha.

Relacionado: [[Ecuación orbital de las cónicas]], [[Maniobra impulsiva]].
