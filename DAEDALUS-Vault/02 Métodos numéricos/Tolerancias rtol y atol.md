---
tags:
  - algoritmia/numerico
aliases:
  - rtol
  - atol
  - Tolerancia vectorial
estado: implementado
fase: 1
---

# Tolerancias rtol y atol

Los integradores adaptativos aceptan un paso si el error estimado de cada componente cumple:

$$|\text{err}_i| \le \text{atol}_i + \text{rtol}\cdot|y_i|$$

- **rtol** (relativa) domina cuando $|y_i|$ es grande.
- **atol** (absoluta) domina cuando $|y_i| \to 0$ y evita exigir precisión infinita cerca de cero.

## Por qué `atol` es vectorial en DAEDALUS
En LEO, $r \sim 7\times10^{6}$ m y $v \sim 7.5\times10^{3}$ m/s: **tres órdenes de magnitud** de diferencia y unidades distintas. Un `atol` escalar sería absurdamente laxo para la velocidad o absurdamente estricto para la posición.

| Componentes | `atol` | Significado |
|---|---|---|
| $x, y, z$ | 1e-3 m | Milímetro |
| $v_x, v_y, v_z$ | 1e-6 m/s | Micrómetro por segundo |

> [!note] Mismo principio en los criterios de aceptación
> Los errores de posición y velocidad se evalúan **por separado** (p. ej. < 1 m y < 0.001 m/s), nunca combinados en una norma única. Ver [[Convención de unidades]].

> [!warning] Suelo de `rtol`
> El propagador rechaza `rtol < 100·eps_float64`: por debajo, el redondeo domina y SciPy la recortaría en silencio. Ver [[Análisis de error numérico]].

Relacionado: [[Dormand-Prince RK45]], [[Tolerancias justificadas]].
