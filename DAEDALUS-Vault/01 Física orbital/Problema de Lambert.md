---
tags:
  - fisica/maniobras
  - algoritmia/numerico
aliases:
  - Lambert
estado: planificado
fase: 2
---

# Problema de Lambert

> [!abstract] En una frase
> Dados dos vectores de posición $\mathbf{r}_1$, $\mathbf{r}_2$ y un tiempo de vuelo $\Delta t$, hallar la órbita que los une.

Es un **problema de contorno de dos puntos**. Su solución da las velocidades necesarias en ambos extremos, y de ellas los dos impulsos.

## Por qué es estructuralmente importante
**Reduce la dimensión del espacio de búsqueda.** En vez de optimizar magnitudes y direcciones de impulso, se optimizan solo los instantes de salida y llegada; Lambert resuelve el resto. Una transferencia de dos impulsos pasa típicamente de **seis a dos** variables de decisión.

## Coste
Exige un solucionador robusto: transferencias de múltiples revoluciones y configuraciones casi colineales. Trabajo no trivial → requisito **opcional** RF-2.07.

Referencia: Battin, *An Introduction to the Mathematics and Methods of Astrodynamics*.

Relacionado: [[Variables de decisión y normalización]], [[Transferencia de Hohmann]], [[Métodos directos de transcripción]].
