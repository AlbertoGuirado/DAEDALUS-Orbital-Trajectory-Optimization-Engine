---
tags:
  - fisica/invariantes
aliases:
  - Vector de Laplace-Runge-Lenz
  - LRL
  - eccentricity_vector
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/mechanics.py
---

# Vector de excentricidad

Vector de Laplace-Runge-Lenz normalizado:

$$\mathbf{e} = \frac{\mathbf{v} \times \mathbf{h}}{\mu} - \frac{\mathbf{r}}{r} = \frac{1}{\mu}\left[\left(v^{2}-\frac{\mu}{r}\right)\mathbf{r} - (\mathbf{r}\cdot\mathbf{v})\,\mathbf{v}\right]$$

- **Apunta del foco al periapsis.**
- Su **módulo es la excentricidad** $e$.

> [!important] El detector más sensible
> Su conservación es **específica del potencial $1/r$**. Cualquier perturbación (o cualquier error del integrador) la rompe antes que a los demás invariantes. Es el primero que se degrada cuando algo va mal: detecta la rotación espuria de la línea de apsides.

## Detalle de validación
El criterio de la Fase 1 usa **error absoluto** (< 10⁻⁷), no relativo, porque en una órbita circular $e = 0$ y el error relativo dividiría por cero.

Relacionado: [[Energía específica]], [[Momento angular específico]], [[Elementos orbitales clásicos]], [[Perturbación J2]] (rotación de apsides).
