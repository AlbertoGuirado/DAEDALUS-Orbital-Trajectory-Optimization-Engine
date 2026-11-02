---
tags:
  - validacion
  - metodologia
aliases:
  - Oráculo de validación
  - Invariantes
  - Cantidades conservadas
estado: implementado
fase: 1
---

# Invariantes como oráculo de validación

> [!abstract] En una frase
> Una cantidad que la física conserva **exactamente** y que el integrador solo conserva **aproximadamente** es una medida directa y honesta del error acumulado.

## Por qué el problema de dos cuerpos es un regalo
En simulación científica rara vez existe la solución exacta. El [[Problema de dos cuerpos]] la tiene (cónicas) y además ofrece **tres** cantidades conservadas. Se puede construir una batería de validación que no depende de valores tabulados externos **ni de resultados previos del propio código**.

## Cada invariante detecta una familia de errores
| Invariante | Detecta especialmente bien |
|---|---|
| [[Energía específica]] | Deriva secular, errores en la magnitud de la aceleración |
| [[Momento angular específico]] — dirección | Errores que sacan la trayectoria de su plano |
| [[Momento angular específico]] — módulo | Errores en la velocidad transversal |
| [[Vector de excentricidad]] | Rotación espuria de la línea de apsides |

> [!warning] Por qué los tres y no uno
> Un integrador puede conservar bien la energía y a la vez hacer precesar la línea de apsides de forma no física.

## Jerarquía de oráculos en DAEDALUS
1. **Solución analítica exacta**: órbita circular; elipse por [[Ecuación de Kepler]] resuelta con [[Método de bisección]].
2. **Invariantes conservados**: energía, $\mathbf{h}$, $\mathbf{e}$.
3. **Relaciones físicas independientes**: [[Ecuación de vis-viva]], radios de apsides, [[Tercera ley de Kepler|cierre del periodo]].
4. **Consistencia interna**: ida y vuelta de elementos, acuerdo RK4 ↔ RK45, [[Orden de convergencia]].

> [!danger] Principio vinculante (RNF-14)
> Un test **nunca** compara contra valores producidos antes por el propio código. Eso es un test de regresión, no de corrección.

## Cuando dejen de servir
Con perturbaciones, $\mathbf{e}$ deja de conservarse; con arrastre o empuje, tampoco la energía. Nuevos oráculos: tasas seculares de [[Perturbación J2]], límite impulsivo del [[Encendido de duración finita]], [[Ecuación de Tsiolkovsky]].

Relacionado: [[Batería de validación de la Fase 1]], [[Niveles de evidencia]].
