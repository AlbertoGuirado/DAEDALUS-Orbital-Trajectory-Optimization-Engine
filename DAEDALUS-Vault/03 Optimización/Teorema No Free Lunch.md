---
tags:
  - algoritmia/optimizacion
  - metodologia
aliases:
  - No Free Lunch
  - NFL
estado: teoria
---

# Teorema No Free Lunch

> [!quote] Wolpert y Macready, 1997
> Promediado sobre **todos** los problemas posibles, ningún algoritmo de optimización supera a otro.

## Consecuencia práctica
La pregunta **"¿cuál es el mejor algoritmo?" no tiene respuesta en abstracto**. Solo tiene sentido preguntar cuál funciona mejor:
- en **esta clase concreta** de problemas,
- con **este presupuesto** de evaluaciones,
- según **estos criterios**.

## Por qué está en el corazón de DAEDALUS
- Justifica el principio 2: **no afirmar superioridad sin benchmark** ([[Principios de ingeniería]]).
- Justifica la [[Metodología de benchmark]]: múltiples semillas, escenarios fijados de antemano, presupuesto idéntico.
- Obliga a **declarar la clase de problemas** sobre la que se concluye (amenaza de validez: conclusión generalizada en exceso).

Relacionado: [[Differential Evolution]], [[Tests estadísticos no paramétricos]].
