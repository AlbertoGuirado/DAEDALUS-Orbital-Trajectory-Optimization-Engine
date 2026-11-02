---
tags:
  - fisica/perturbaciones
aliases:
  - J2
  - Achatamiento terrestre
  - Regresión nodal
  - Órbita heliosíncrona
estado: planificado
fase: 5
---

# Perturbación J2

Efecto del **achatamiento terrestre**, la perturbación dominante en órbita baja.

$$\mathbf{a}_{J_2} = -\frac{3J_2\mu R_\oplus^{2}}{2r^{5}}\begin{bmatrix}x\left(1-\dfrac{5z^{2}}{r^{2}}\right)\\[2mm] y\left(1-\dfrac{5z^{2}}{r^{2}}\right)\\[2mm] z\left(3-\dfrac{5z^{2}}{r^{2}}\right)\end{bmatrix}$$

## Efectos seculares (no promedian a cero)
$$\dot\Omega = -\frac{3}{2}\,n\,J_2\left(\frac{R_\oplus}{p}\right)^{2}\cos i \qquad \text{(regresión nodal)}$$

$$\dot\omega = \frac{3}{4}\,n\,J_2\left(\frac{R_\oplus}{p}\right)^{2}\left(5\cos^{2} i - 1\right) \qquad \text{(rotación de apsides)}$$

## Consecuencias de ingeniería
> [!example] Órbitas heliosíncronas
> Eligiendo la inclinación, el plano precesa exactamente una vuelta al año → geometría constante respecto al Sol.

> [!example] Inclinación crítica ≈ 63.4°
> Con $5\cos^2 i = 1$ la rotación de apsides se anula. Lo explotan las órbitas **Molniya** para mantener el apoapsis sobre latitudes altas.

## Validación
Las tasas seculares son una referencia **analítica independiente** para verificar la implementación numérica: RV-5.02 (reproducir la tasa de regresión nodal).

Relacionado: [[Perturbaciones orbitales]], [[Constantes físicas]], [[Elementos orbitales clásicos]].
