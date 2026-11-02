---
tags:
  - algoritmia/optimizacion
  - algoritmia/control-optimo
aliases:
  - Transcripción directa
  - Disparo múltiple
  - Colocación
  - Pseudoespectral
  - Shooting
estado: planificado
fase: 3
---

# Métodos directos de transcripción

Discretizan el problema de [[Control óptimo]] continuo en un problema de **programación no lineal** (NLP) de dimensión finita, resoluble con [[Métodos de optimización local|SLSQP]] u otro NLP solver.

| Método | Idea |
|---|---|
| **Disparo simple** | Parametrizar el control e integrar toda la trayectoria de una vez. Muy sensible: un cambio al inicio se amplifica al final |
| **Disparo múltiple** | Partir el intervalo en tramos integrados por separado + restricciones de continuidad en las uniones. Mucho mejor condicionado |
| **Colocación** | Discretizar estado y control; la dinámica se impone como **restricción algebraica** en puntos de colocación (sin integrador) |
| **Pseudoespectral** | Colocación con polinomios globales en nodos de cuadratura de Gauss; convergencia espectral para soluciones suaves |

## Ventajas e inconvenientes
- ✅ Radio de convergencia amplio, sin costados, restricciones naturales.
- ❌ Menor precisión que los indirectos; NLP potencialmente de dimensión muy alta.

> [!tip] Encaje con la arquitectura
> Se reduce a un NLP → se resuelve con las mismas herramientas que el resto del proyecto. El disparo simple/múltiple reutiliza directamente el propagador existente.

Referencia: Betts, *Practical Methods for Optimal Control Using Nonlinear Programming*.

Relacionado: [[Principio del máximo de Pontryagin]], [[Problema de Lambert]].
