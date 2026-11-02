---
tags:
  - algoritmia/optimizacion
aliases:
  - Landscape
  - Espacio de búsqueda
estado: teoria
fase: 3
---

# Paisaje de optimización orbital

Las características del problema de transferencias determinan qué algoritmos tienen sentido:

| Característica | Implicación algorítmica |
|---|---|
| **No convexo**, múltiples óptimos locales | Un método puramente local queda atrapado → búsqueda global |
| **Objetivo caro** (cada evaluación = una propagación) | El presupuesto de evaluaciones es la métrica justa |
| **Sin gradiente analítico** | Métodos sin derivadas, o diferencias finitas ($n+1$ evaluaciones por gradiente) |
| Posibles **discontinuidades** por restricciones | Complica los métodos de gradiente |
| Variables de **escalas dispares** | [[Variables de decisión y normalización\|Normalización]] obligatoria |
| **Topología periódica** en los ángulos | Operadores que respeten la periodicidad |

> [!success] Conclusión
> Estas características justifican **metaheurísticas poblacionales** y, en particular, [[Differential Evolution]] como primer algoritmo. Después, refinamiento local en esquema [[Algoritmos meméticos|memético]].

## Alternativa: diferenciación automática
Si el propagador fuese diferenciable, los gradientes serían exactos y baratos. Exige reescribirlo → fuera del alcance actual.

Relacionado: [[Teorema No Free Lunch]], [[Métodos de optimización local]].
