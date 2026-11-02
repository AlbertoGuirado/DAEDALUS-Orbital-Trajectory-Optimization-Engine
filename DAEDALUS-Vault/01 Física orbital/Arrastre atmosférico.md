---
tags:
  - fisica/perturbaciones
aliases:
  - Drag
  - Coeficiente balístico
estado: planificado
fase: 5
---

# Arrastre atmosférico

$$\mathbf{a}_{\text{arrastre}} = -\frac{1}{2}\rho\,\frac{C_D A}{m}\,\lVert\mathbf{v}_{\text{rel}}\rVert\,\mathbf{v}_{\text{rel}}$$

- $\mathbf{v}_{\text{rel}}$: velocidad relativa a la atmósfera, que **corrota con la Tierra**.
- $C_D A/m$: inverso del coeficiente balístico.

## Naturaleza
**No conservativa.** Disipa energía de forma monótona, reduce el semieje mayor y circulariza la órbita. No admite oráculo por invariante conservado: se valida por **comportamiento cualitativo** (decaimiento monótono de la [[Energía específica]]).

> [!warning] La dificultad real es la densidad
> $\rho(h)$ varía **órdenes de magnitud** con la actividad solar; es la mayor fuente de incertidumbre en predicción LEO. Un modelo exponencial basta para estudios comparativos, pero esta limitación debe declararse junto a **cualquier** resultado con arrastre (RF-5.05).

Relevante por debajo de ~600 km; despreciable por encima de ~1000 km.

Relacionado: [[Perturbaciones orbitales]], [[Perturbación J2]].
