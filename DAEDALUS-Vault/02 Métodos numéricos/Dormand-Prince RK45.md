---
tags:
  - algoritmia/numerico
  - algoritmia/integracion
aliases:
  - RK45
  - Dormand-Prince
  - DOPRI5
  - solve_ivp
  - Paso adaptativo
estado: implementado
fase: 1
codigo: src/orbital_optimizer/simulation/propagator.py
---

# Dormand-Prince RK45

> [!abstract] En una frase
> Runge-Kutta **embebido** de 7 etapas que calcula a la vez una solución de orden 5 y otra de orden 4; su diferencia estima el error local y gobierna el tamaño del paso.

Es el método por defecto de `scipy.integrate.solve_ivp` (`method="RK45"`), que DAEDALUS usa como **integrador por defecto** para simulación.

## Control adaptativo del paso
$$h_{\text{nuevo}} = h \cdot \sigma \left(\frac{\text{tol}}{\lVert \mathbf{err}\rVert}\right)^{1/5}, \qquad \text{tol}_i = \text{atol}_i + \text{rtol}\cdot\lvert y_i \rvert$$

con $\sigma < 1$ factor de seguridad. Si el error supera la tolerancia, el paso se rechaza y se repite más pequeño.

**FSAL** (*First Same As Last*): la última evaluación de un paso es la primera del siguiente → coste efectivo de **6 evaluaciones por paso aceptado**.

## Configuración en DAEDALUS
| Parámetro | Valor | Por qué |
|---|---|---|
| `rtol` | 1e-10 | Se rechaza < 100·eps en vez de dejar que SciPy la recorte en silencio |
| `atol` | `[1e-3]*3 + [1e-6]*3` | mm en posición, µm/s en velocidad — ver [[Tolerancias rtol y atol]] |
| `max_step` | 60 s | Techo explícito para ejemplos terrestres, configurable |

Las muestras pedidas se **interpolan** sobre la malla adaptativa (salida densa); no alteran los pasos, a diferencia de [[Runge-Kutta 4]].

> [!danger] Tolerancia local ≠ error global
> El control de error es **local** (por paso). Una tolerancia estricta no garantiza un error global pequeño: los errores se acumulan. Por eso la validación mide el error global contra soluciones analíticas. Ver [[Análisis de error numérico]].

Fallo del solver → `RuntimeError`, nunca éxito parcial (`test_solver_failure_is_explicit`).

Relacionado: [[Comparativa RK4 frente a RK45]], [[Resultados de validación de la Fase 1]] (788 evaluaciones para un periodo).
