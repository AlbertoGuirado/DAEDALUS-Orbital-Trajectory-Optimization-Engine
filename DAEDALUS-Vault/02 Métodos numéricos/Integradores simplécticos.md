---
tags:
  - algoritmia/numerico
  - algoritmia/integracion
  - decision-de-diseno
aliases:
  - Verlet
  - Yoshida
  - Integración geométrica
estado: teoria
fase: 1
---

# Integradores simplécticos

> [!abstract] En una frase
> Integradores que preservan la estructura hamiltoniana del sistema, de modo que el error de energía queda **acotado y oscilante** en lugar de derivar secularmente.

Ejemplos: Verlet de velocidad (orden 2), Yoshida (órdenes superiores por composición).

## Por qué NO se adoptan (por ahora)
| Motivo | Detalle |
|---|---|
| Paso fijo | Ineficientes en órbitas excéntricas (mismo problema que RK4) |
| Fuerzas no conservativas | Empuje y arrastre — el objeto de las fases 4-5 — rompen la estructura hamiltoniana que justifica su ventaja |

> [!tip] Decisión con datos, no a priori
> Medir la deriva de energía de [[Dormand-Prince RK45]] en los regímenes de interés y reconsiderar solo si resulta limitante. En la Fase 1, la deriva relativa medida en un periodo fue ~9.4×10⁻¹⁰ ([[Resultados de validación de la Fase 1]]).

Referencia: Hairer, Lubich y Wanner, *Geometric Numerical Integration*.

Relacionado: [[Comparativa RK4 frente a RK45]], [[Energía específica]].
