---
tags:
  - fisica/fundamentos
aliases:
  - μ
  - mu
  - EARTH_MU
  - Parámetro gravitacional
estado: implementado
fase: 1
codigo: src/orbital_optimizer/dynamics/two_body.py
---

# Constantes físicas

| Constante | Valor | Unidad | Fuente | En el código |
|---|---|---|---|---|
| $\mu_\oplus$ | 3.986004418 × 10¹⁴ | m³/s² | WGS-84 / EGM-96 | `EARTH_MU` ✅ |
| $R_\oplus$ | 6 378 137 | m | Radio ecuatorial WGS-84 | Fase 4/5 |
| $J_2$ | 1.08262668 × 10⁻³ | — | EGM-96 | Fase 5 |
| $g_0$ | 9.80665 | m/s² | Definición CGPM | Fase 2 |

## Parámetro gravitacional estándar
$$\mu = G(M + m) \approx GM$$
Se usa μ directamente y no $G$ y $M$ por separado porque μ se conoce con muchísima más precisión que cualquiera de los dos factores. `validate_mu` exige un escalar real, finito y estrictamente positivo; todas las funciones aceptan un `mu` personalizado.

> [!note] $g_0$ no es la gravedad local
> Es una constante de definición que solo sirve para convertir impulso específico en velocidad de escape efectiva en la [[Ecuación de Tsiolkovsky]].

Relacionado: [[Convención de unidades]], [[Problema de dos cuerpos]], [[Perturbación J2]].
