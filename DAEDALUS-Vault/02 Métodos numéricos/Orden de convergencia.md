---
tags:
  - algoritmia/numerico
  - validacion
aliases:
  - Test del factor 16
  - Orden observado
estado: implementado
fase: 1
codigo: tests/test_propagator.py
---

# Orden de convergencia

Un método de orden $p$ tiene error global $E(h) \approx C\,h^{p}$. Con dos pasos $h$ y $h/2$:

$$p_{\text{obs}} = \log_2 \frac{E(h)}{E(h/2)}$$

Para [[Runge-Kutta 4]], $p = 4$ → halvar el paso divide el error por $2^4 = \mathbf{16}$.

## Por qué es un test tan exigente
Un error sutil en los coeficientes del esquema (un 2 donde va un 1, un $h/2$ mal puesto) suele dar un método que **funciona** pero de orden menor. Un test de "error pequeño" no lo detecta; el test de orden sí.

## En DAEDALUS (RV-1.09)
- Pasos $T/100$, $T/200$, $T/400$.
- Criterio: orden observado entre **3.7 y 4.5** — entorno asintótico de orden 4, admitiendo términos de orden superior que aún no son despreciables.
- Sin `sample_times`, porque el muestreo altera la malla de RK4.
- Test: `test_rk4_fourth_order_convergence`.

> [!note] Límite del test
> Si se reduce demasiado el paso se alcanza el suelo de redondeo y el orden observado se desploma. Ver [[Análisis de error numérico]].

Relacionado: [[Batería de validación de la Fase 1]], [[Tolerancias justificadas]].
