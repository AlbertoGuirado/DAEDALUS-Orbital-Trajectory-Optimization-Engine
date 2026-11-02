---
tags:
  - algoritmia/numerico
aliases:
  - SINGULAR_TOL
  - Épsilon de máquina
  - eps
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/elements.py
---

# Umbral de degeneración numérica

$$\texttt{SINGULAR\_TOL} = 32 \cdot \varepsilon_{\text{float64}} \approx 7.1\times10^{-15}$$

donde $\varepsilon_{\text{float64}} = 2^{-52} \approx 2.22\times10^{-16}$ es el **épsilon de máquina**: el menor número tal que $1 + \varepsilon \ne 1$ en doble precisión.

## Dónde se aplica
Cocientes **adimensionales**, para que el umbral no dependa de la escala:
| Cantidad | Si está por debajo… |
|---|---|
| $e$ | Órbita circular |
| $\lvert\mathbf{n}\rvert/\lvert\mathbf{h}\rvert$ | Órbita ecuatorial |
| $\lvert\mathbf{h}\rvert/(\lvert\mathbf{r}\rvert\lvert\mathbf{v}\rvert)$ | Órbita radial → se rechaza |

> [!important] Filosofía
> El umbral **solo absorbe redondeo**; no es una tolerancia física. Una órbita con $e = 10^{-6}$ es elíptica de verdad y se trata como tal. Cerca de la singularidad, los ángulos individuales están mal condicionados, así que los tests comparan el **estado reconstruido**, no los ángulos.

## Otros usos del épsilon en el código
- `rtol` mínima = $100\varepsilon$ en el propagador.
- Patrón general: elegir umbrales como múltiplos pequeños de $\varepsilon$ cuando el objetivo es distinguir "cero exacto" de "cero por redondeo".

Relacionado: [[Casos degenerados de los elementos orbitales]], [[Análisis de error numérico]].
