---
tags:
  - algoritmia/numerico
  - algoritmia/integracion
aliases:
  - RK4
  - Runge-Kutta clásico
estado: implementado
fase: 1
codigo: src/orbital_optimizer/simulation/propagator.py
---

# Runge-Kutta 4

> [!abstract] En una frase
> Integrador explícito de paso fijo que promedia cuatro pendientes por paso para lograr error global de orden 4.

## Esquema
$$\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(t_n,\; \mathbf{y}_n)\\
\mathbf{k}_2 &= \mathbf{f}\!\left(t_n+\tfrac{h}{2},\; \mathbf{y}_n + \tfrac{h}{2}\mathbf{k}_1\right)\\
\mathbf{k}_3 &= \mathbf{f}\!\left(t_n+\tfrac{h}{2},\; \mathbf{y}_n + \tfrac{h}{2}\mathbf{k}_2\right)\\
\mathbf{k}_4 &= \mathbf{f}(t_n+h,\; \mathbf{y}_n + h\,\mathbf{k}_3)\\
\mathbf{y}_{n+1} &= \mathbf{y}_n + \frac{h}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)
\end{aligned}$$

Intuición: $\mathbf{k}_1$ pendiente al inicio, $\mathbf{k}_2$ y $\mathbf{k}_3$ estimaciones en el punto medio, $\mathbf{k}_4$ al final; la media ponderada (1-2-2-1, como Simpson) cancela los términos de Taylor hasta $h^4$.

## Propiedades
| | |
|---|---|
| Error local | $\mathcal{O}(h^5)$ |
| Error global | $\mathcal{O}(h^4)$ → al **halvar** el paso, el error cae ~**16×** |
| Coste | 4 evaluaciones de $\mathbf{f}$ por paso |
| Control de error | Ninguno |
| Simpléctico | No → la energía deriva secularmente en propagaciones muy largas |

## En DAEDALUS (`_rk4`)
- Implementado **a mano**, paso = `max_step`.
- El paso se **acorta** para aterrizar exactamente en el instante final y en cada instante de muestreo pedido.

> [!warning] Muestrear altera la malla (RF-1.20)
> Pedir `sample_times` cambia los pasos y, por tanto, el error. Los estudios de convergencia deben hacerse con `sample_times=None`.

- `rtol`/`atol` se validan pero **RK4 no los usa**.
- Uso recomendado: **estudios de [[Orden de convergencia]]**, donde el paso prescrito es justo lo que se quiere controlar. Test: `test_rk4_fourth_order_convergence` (orden observado entre 3.7 y 4.5 con $T/100$, $T/200$, $T/400$).

Relacionado: [[Dormand-Prince RK45]], [[Comparativa RK4 frente a RK45]], [[Análisis de error numérico]].
