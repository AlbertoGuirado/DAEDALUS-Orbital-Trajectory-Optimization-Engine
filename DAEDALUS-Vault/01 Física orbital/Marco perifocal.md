---
tags:
  - fisica/geometria
aliases:
  - PQW
  - Plano perifocal
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/elements.py
---

# Marco perifocal

Marco ligado a la órbita: eje $\hat{p}$ hacia el periapsis, $\hat{q}$ a 90° en el sentido del movimiento, $\hat{w}$ paralelo a $\mathbf{h}$. En él la órbita es plana y el estado es trivial:

$$\mathbf{r}_{pqw} = \frac{p}{1+e\cos\nu}\begin{bmatrix}\cos\nu\\ \sin\nu\\ 0\end{bmatrix}, \qquad \mathbf{v}_{pqw} = \sqrt{\frac{\mu}{p}}\begin{bmatrix}-\sin\nu\\ e+\cos\nu\\ 0\end{bmatrix}$$

## Rotación al marco inercial
Tres rotaciones elementales sucesivas (Euler 3-1-3):

$$\mathbf{r} = R_z(\Omega)\,R_x(i)\,R_z(\omega)\;\mathbf{r}_{pqw}$$

y lo mismo para $\mathbf{v}$. Es el núcleo de `elements_to_state`.

> [!note] Convención de signos
> El documento 03 escribe la matriz como $R_z(-\Omega)R_x(-i)R_z(-\omega)$ con matrices de rotación *de ejes*; `docs/phase1.md` la escribe como $R_z(\Omega)R_x(i)R_z(\omega)$ con matrices de rotación *de vectores*. Son la misma transformación.

Relacionado: [[Elementos orbitales clásicos]], [[Ecuación orbital de las cónicas]], [[Marco de referencia inercial]].
