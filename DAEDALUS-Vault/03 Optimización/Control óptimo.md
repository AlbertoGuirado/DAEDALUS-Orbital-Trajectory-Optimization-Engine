---
tags:
  - algoritmia/optimizacion
  - algoritmia/control-optimo
aliases:
  - Optimal control
estado: planificado
fase: 3
codigo: optimization/optimal_control.py (vacío)
---

# Control óptimo

> [!abstract] En una frase
> Aborda el problema en forma **continua**: hallar la ley de empuje $\mathbf{u}(t)$ completa que minimiza un funcional, sujeta a la dinámica.

$$\min_{\mathbf{u}(\cdot)} \; J = \underbrace{\phi(\mathbf{x}(t_f), t_f)}_{\text{coste terminal}} + \int_{t_0}^{t_f} \underbrace{L(\mathbf{x},\mathbf{u},t)}_{\text{coste corriente}}\,dt \quad \text{s.a.} \quad \dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u},t)$$

Frente a la optimización paramétrica (unos pocos impulsos), aquí la incógnita es una **función**: dimensión infinita.

## Dos familias
| | [[Principio del máximo de Pontryagin\|Indirectos]] | [[Métodos directos de transcripción\|Directos]] |
|---|---|---|
| Idea | "Optimizar y luego discretizar": condiciones necesarias → problema de contorno | "Discretizar y luego optimizar": → programación no lineal (NLP) |
| Precisión | Alta | Menor |
| Radio de convergencia | **Muy estrecho** | Amplio |
| Necesita | Derivar costados y estimarlos inicialmente | Solo un NLP solver |

> [!success] Recomendación del proyecto
> Si se aborda (RF-3.15, opcional), por vía **directa**. La teoría indirecta se usa para **verificar** soluciones, mediante el vector primer.

Relacionado: [[Encendido de duración finita]], [[Formulación del problema de optimización]].
