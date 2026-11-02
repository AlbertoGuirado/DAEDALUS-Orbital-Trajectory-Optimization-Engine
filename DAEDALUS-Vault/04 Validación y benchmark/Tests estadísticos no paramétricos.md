---
tags:
  - validacion
  - estadistica
aliases:
  - Wilcoxon
  - Mann-Whitney
  - Friedman
  - Holm
  - Nemenyi
  - Comparaciones múltiples
estado: planificado
fase: 3
---

# Tests estadísticos no paramétricos

Los resultados de metaheurísticas **no siguen una distribución normal** (asimétricos, con colas) → los tests paramétricos (t de Student, ANOVA) no son apropiados.

| Test | Uso |
|---|---|
| **Wilcoxon** de rangos con signo | Dos algoritmos sobre los **mismos** problemas (datos emparejados) |
| **Mann-Whitney U** | Dos muestras **independientes** |
| **Friedman** | **Varios** algoritmos sobre varios problemas |
| Post-hoc de **Holm** o **Nemenyi** | Comparaciones por pares tras un Friedman significativo |

Todos en `scipy.stats` (`wilcoxon`, `mannwhitneyu`, `friedmanchisquare`).

## Corrección por comparaciones múltiples
5 algoritmos → 10 comparaciones por pares. Sin corregir, a α = 0.05:

$$P(\ge 1 \text{ falso positivo}) = 1 - 0.95^{10} \approx 40\%$$

> [!danger] La corrección no es opcional

## Tamaño del efecto
Significación ≠ relevancia. Con suficientes ejecuciones, **cualquier** diferencia mínima sale significativa. Reportar siempre la **magnitud** de la diferencia y valorar si es operativamente relevante (¿importan 0.3 m/s de Δv?).

Referencia: Derrac et al., 2011, *A Practical Tutorial on the Use of Nonparametric Statistical Tests…*

Relacionado: [[Metodología de benchmark]], [[Teorema No Free Lunch]].
