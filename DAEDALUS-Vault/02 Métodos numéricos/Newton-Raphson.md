---
tags:
  - algoritmia/numerico
  - algoritmia/raices
aliases:
  - Método de Newton
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/kepler.py
---

# Newton-Raphson

Método iterativo para resolver $g(x) = 0$ linealizando en cada iteración:

$$x_{k+1} = x_k - \frac{g(x_k)}{g'(x_k)}$$

| Propiedad | |
|---|---|
| Convergencia | **Cuadrática** cerca de la raíz (los dígitos correctos se duplican por iteración) |
| Requiere | Derivada $g'$ y una semilla razonable |
| Riesgo | Divergir u oscilar si $g' \approx 0$ o la semilla es mala |

## Aplicación canónica: [[Ecuación de Kepler]]
$g(E) = E - e\sin E - M$, $g'(E) = 1 - e\cos E$:

$$E_{k+1} = E_k - \frac{E_k - e\sin E_k - M}{1 - e\cos E_k}$$

Semilla $E_0 = M$ para $e$ moderada; $E_0 = \pi$ cuando $e \to 1$ (ahí $g'$ se acerca a 0 en el periapsis).

## Conexión con optimización
Aplicado a $\nabla J = 0$ da el método de Newton de optimización; sus variantes cuasi-Newton (BFGS, L-BFGS-B) están en [[Métodos de optimización local]].

> [!info] En DAEDALUS
> El solver integrado de Kepler usa Newton salvaguardado: cuando la propuesta sale del intervalo de raíz se sustituye por bisección. Los oráculos independientes de tests usan solo [[Método de bisección]]. Ver [[Cierre documental de la Fase 1]].
