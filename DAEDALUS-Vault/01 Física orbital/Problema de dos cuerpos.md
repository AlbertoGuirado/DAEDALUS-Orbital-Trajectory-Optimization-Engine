---
tags:
  - fisica/dinamica
aliases:
  - Dos cuerpos
  - Two-body
  - Ecuación de Newton orbital
estado: implementado
fase: 1
codigo: src/orbital_optimizer/dynamics/two_body.py
---

# Problema de dos cuerpos

> [!abstract] En una frase
> La ecuación fundamental del proyecto: el cuerpo menor se acelera hacia el central con intensidad inversamente proporcional al cuadrado de la distancia.

$$\ddot{\mathbf{r}} = -\frac{\mu\,\mathbf{r}}{\lVert\mathbf{r}\rVert^{3}}$$

Para la posición relativa entre dos masas puntuales, restar sus ecuaciones de Newton da $\ddot{\mathbf r}=-G(M+m)\mathbf r/r^3$. En la aproximación de primario fijo y $m \ll M$, se usa $\mu\simeq GM$. DAEDALUS recibe μ constante y expresa el estado respecto al primario, sin integrar su movimiento baricéntrico. Es una EDO vectorial de **segundo orden**; para integrarla se reduce a primer orden en el [[Vector de estado]]:

$$\dot{\mathbf{y}} = \left[v_x,\; v_y,\; v_z,\; -\frac{\mu x}{r^{3}},\; -\frac{\mu y}{r^{3}},\; -\frac{\mu z}{r^{3}}\right]$$

En el código: `two_body_dynamics(t, state, mu)`, firma compatible con `scipy.integrate.solve_ivp`.

## Propiedades del campo vectorial
| Propiedad | Consecuencia |
|---|---|
| **Autónomo** ($f$ no depende de $t$) | El argumento `t` solo existe por compatibilidad con SciPy |
| **Singular en el origen** ($a \sim r^{-2}$) | `validate_state` rechaza radio nulo; no existe un detector de colisiones o de contacto con una superficie |
| **Conservativo** ($U = -\mu/r$) | Existen [[Energía específica]], [[Momento angular específico]] y [[Vector de excentricidad]] conservados |
| **Solución analítica cerrada** (cónicas) | Se puede validar el integrador contra una verdad exacta — algo rarísimo en simulación |

## Condicionamiento numérico
Cerca del periapsis de órbitas muy excéntricas, velocidad y aceleración crecen bruscamente: un paso fijo adecuado en el apoapsis es insuficiente en el periapsis. Motiva preferir [[Dormand-Prince RK45]] y limitar el paso máximo.

## Dominio de la Fase 1
Elipses no radiales y ligadas: $a > 0$, $0 \le e < 1$. **No** hay superficie del cuerpo central, colisiones, empuje, masa variable, arrastre ni J2. Una trayectoria que atraviese el centro no es físicamente continuable con este modelo.

Relacionado: [[Ecuación orbital de las cónicas]], [[Perturbaciones orbitales]], [[Marco de referencia inercial]].

## Derivación de invariantes y unidades

En vuelo libre, $d(v^2/2-\mu/r)/dt=\mathbf v\cdot\ddot{\mathbf r}
+\mu(\mathbf r\cdot\mathbf v)/r^3=0$, y
$d(\mathbf r\times\mathbf v)/dt=\mathbf r\times\ddot{\mathbf r}=0$.
Así se obtienen energía específica (m²/s²) y momento angular específico
(m²/s). Posición en m, velocidad en m/s, tiempo en s y μ en m³/s².
Son invariantes del modelo continuo entre impulsos, no una garantía de
conservación exacta del integrador ni a través de un salto de velocidad.

**Implemented:** campo newtoniano, invariantes y propagación elíptica.
**Known Limitations:** la validación del campo solo exige radio no nulo y
estado real finito; el dominio elíptico adicional pertenece al propagador.
**Planned:** perturbaciones y modelos de propulsión extendidos.
Fuente: [Relative Motion in the Two-Body Problem](https://orbital-mechanics.space/the-n-body-problem/two-body-relative-motion.html),
consultada el 2026-11-02. Ver [[Validación de la Fase 2 e integración]].
